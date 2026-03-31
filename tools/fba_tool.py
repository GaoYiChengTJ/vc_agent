"""FBA simulation tools for metabolic engineering.

Design: stateful working model.
- Modification actions (reset, add_pathway, knockout, overexpress, media)
  change the working model in place.
- Analysis actions (maximize, envelope) read the working model via
  ``with model:`` context so they do NOT corrupt state.
- GEM query functions (search_model, query_gpr, get_metabolite_reactions)
  load independent copies and are unaffected by the working model.
"""

from pathlib import Path

import cobra

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

AVAILABLE_MODELS = {
    "e_coli_core": DATA_DIR / "e_coli_core.xml",
    "iMM904": DATA_DIR / "iMM904.xml",
}

# Inorganic exchange reactions that should stay open when switching carbon source
_INORGANIC_EXCHANGES = frozenset({
    "EX_o2_e", "EX_h2o_e", "EX_h_e",
    "EX_nh4_e", "EX_pi_e", "EX_co2_e",
    "EX_so4_e", "EX_fe2_e", "EX_k_e",
    "EX_na1_e", "EX_cl_e",
})

# ── Model loading (cached) ──────────────────────────────────────────────

_MODEL_CACHE: dict = {}


def _load_model(model_name: str):
    """Return a fresh copy of the specified metabolic model.

    The parsed model is cached after the first load so that subsequent
    calls only pay the cost of ``model.copy()`` (~0.1 s) instead of
    re-parsing the SBML XML (~1.5 s).
    """
    if model_name not in _MODEL_CACHE:
        path = AVAILABLE_MODELS.get(model_name)
        if path is None:
            raise ValueError(
                f"Unknown model '{model_name}'. "
                f"Available: {', '.join(sorted(AVAILABLE_MODELS))}"
            )
        _MODEL_CACHE[model_name] = cobra.io.read_sbml_model(str(path))
    return _MODEL_CACHE[model_name].copy()


def _get_biomass_reaction_id(model) -> str:
    """Get the biomass reaction ID from the model's current objective."""
    for rxn in model.reactions:
        if rxn.objective_coefficient != 0:
            return rxn.id
    raise ValueError("No biomass/objective reaction found in model.")


# ═══════════════════════════════════════════════════════════════════════════
# Stateful working model
# ═══════════════════════════════════════════════════════════════════════════

_WORKING_MODEL = None
_WORKING_MODEL_NAME: str | None = None


def _require_model():
    """Return the working model or raise if none loaded."""
    if _WORKING_MODEL is None:
        raise RuntimeError("No model loaded. Call fba(action='reset') first.")
    return _WORKING_MODEL


# ── 1. reset ─────────────────────────────────────────────────────────────

def fba_reset(model_name: str) -> str:
    """Load a fresh copy of a GEM as the working model."""
    global _WORKING_MODEL, _WORKING_MODEL_NAME
    try:
        _WORKING_MODEL = _load_model(model_name)
        _WORKING_MODEL_NAME = model_name
    except Exception as exc:
        return f"Error loading model: {exc}"

    n_mets = len(_WORKING_MODEL.metabolites)
    n_rxns = len(_WORKING_MODEL.reactions)
    n_genes = len(_WORKING_MODEL.genes)
    try:
        biomass_id = _get_biomass_reaction_id(_WORKING_MODEL)
    except ValueError:
        biomass_id = "(none)"

    return (
        f"Model '{model_name}' loaded as working model.\n"
        f"Metabolites: {n_mets}  Reactions: {n_rxns}  Genes: {n_genes}\n"
        f"Biomass reaction: {biomass_id}"
    )


# ── 2. add_pathway ──────────────────────────────────────────────────────

def fba_add_pathway(reactions: list[dict]) -> str:
    """Add heterologous reactions to the working model (persistent)."""
    try:
        model = _require_model()
    except RuntimeError as exc:
        return str(exc)

    if not reactions:
        return "Error: reactions list must be non-empty."

    for i, rxn_dict in enumerate(reactions):
        if not isinstance(rxn_dict, dict):
            return f"Error: reactions[{i}] is not a dict."
        if "id" not in rxn_dict or "reaction_string" not in rxn_dict:
            return (
                f"Error: reactions[{i}] must have keys 'id' and 'reaction_string'. "
                f"Got keys: {list(rxn_dict.keys())}"
            )

    added_ids = []
    try:
        for rxn_dict in reactions:
            rid = rxn_dict["id"]
            if rid in model.reactions:
                added_ids.append(f"{rid} (already exists, skipped)")
                continue
            new_rxn = cobra.Reaction(rid)
            model.add_reactions([new_rxn])
            new_rxn.build_reaction_from_string(rxn_dict["reaction_string"])
            added_ids.append(rid)
    except Exception as exc:
        return (
            f"Error adding reaction '{rxn_dict['id']}': {exc}\n"
            f"Successfully added before failure: {added_ids}"
        )

    lines = [f"Added {len(added_ids)} reaction(s) to working model:"]
    for rd in reactions:
        status = "skipped" if rd["id"] + " (already exists, skipped)" in added_ids else "added"
        lines.append(f"  {rd['id']:25s}  {rd['reaction_string']}")
    lines.append(f"\nModel now has {len(model.reactions)} reactions.")
    return "\n".join(lines)


# ── 3. knockout ─────────────────────────────────────────────────────────

def fba_knockout(knockouts: list[str]) -> str:
    """Knock out genes in the working model (persistent)."""
    try:
        model = _require_model()
    except RuntimeError as exc:
        return str(exc)

    model_gene_ids = {g.id for g in model.genes}
    unknown = [g for g in knockouts if g not in model_gene_ids]
    if unknown:
        return (
            f"Error: gene(s) {unknown} not found in model. "
            f"Available genes (first 20): "
            f"{sorted(model_gene_ids)[:20]} ..."
        )

    try:
        cobra.manipulation.knock_out_model_genes(model, knockouts)
    except Exception as exc:
        return f"Error during knockout: {exc}"

    # Report which reactions were disabled
    disabled = []
    for g_id in knockouts:
        gene = model.genes.get_by_id(g_id)
        for rxn in gene.reactions:
            if rxn.lower_bound == 0 and rxn.upper_bound == 0:
                disabled.append(f"  {rxn.id:25s}  (bounds forced to [0, 0])")

    lines = [f"Knocked out {len(knockouts)} gene(s): {knockouts}"]
    if disabled:
        lines.append(f"Disabled reactions ({len(disabled)}):")
        lines.extend(disabled)
    else:
        lines.append("No reactions fully disabled (genes may have OR-logic redundancy in GPR).")
    return "\n".join(lines)


# ── 4. overexpress ──────────────────────────────────────────────────────

def fba_overexpress(gene_id: str, min_flux: float) -> str:
    """Force minimum flux on gene-associated reactions (persistent)."""
    try:
        model = _require_model()
    except RuntimeError as exc:
        return str(exc)

    model_gene_ids = {g.id for g in model.genes}
    if gene_id not in model_gene_ids:
        return (
            f"Error: gene '{gene_id}' not found in model. "
            f"Available genes (first 20): "
            f"{sorted(model_gene_ids)[:20]} ..."
        )

    gene = model.genes.get_by_id(gene_id)
    associated_rxns = list(gene.reactions)
    if not associated_rxns:
        return f"Gene '{gene_id}' has no associated reactions in the model."

    changes = []
    for rxn in associated_rxns:
        old_lb = rxn.lower_bound
        rxn.lower_bound = max(rxn.lower_bound, min_flux)
        changes.append(
            f"  {rxn.id:25s}  lb: {old_lb:.2f} → {rxn.lower_bound:.2f}  "
            f"GPR: {rxn.gene_reaction_rule}"
        )

    lines = [
        f"Overexpressed gene '{gene_id}' ({gene.name or 'N/A'}):",
        f"Forced min flux ≥ {min_flux} on {len(associated_rxns)} reaction(s):",
    ]
    lines.extend(changes)
    return "\n".join(lines)


# ── 5. media ────────────────────────────────────────────────────────────

def fba_media(carbon_source: str, oxygen_lower_bound: float) -> str:
    """Change growth medium on the working model (persistent)."""
    try:
        model = _require_model()
    except RuntimeError as exc:
        return str(exc)

    ex_ids = sorted(r.id for r in model.reactions if r.id.startswith("EX_"))
    if carbon_source not in ex_ids:
        return (
            f"Error: '{carbon_source}' is not a valid exchange reaction. "
            f"Available exchange reactions (first 20): {', '.join(ex_ids[:20])}"
        )

    # Close all organic carbon exchanges except the specified one
    closed = []
    for rxn in model.reactions:
        if rxn.id.startswith("EX_") and rxn.id != carbon_source:
            if rxn.lower_bound < 0 and rxn.id not in _INORGANIC_EXCHANGES:
                rxn.lower_bound = 0.0
                closed.append(rxn.id)

    model.reactions.get_by_id(carbon_source).lower_bound = -10.0
    model.reactions.get_by_id("EX_o2_e").lower_bound = oxygen_lower_bound

    if oxygen_lower_bound == 0:
        condition = "anaerobic"
    elif oxygen_lower_bound > -1000:
        condition = f"micro-aerobic (O2 lb={oxygen_lower_bound})"
    else:
        condition = "aerobic"

    return (
        f"Media changed on working model:\n"
        f"  Carbon source: {carbon_source} (lb=-10.0)\n"
        f"  Oxygen: {condition}\n"
        f"  Closed {len(closed)} other organic exchanges."
    )


# ── 6. maximize ─────────────────────────────────────────────────────────

def fba_maximize(target_reaction: str, min_biomass_fraction: float = 0.0) -> str:
    """Maximize flux through target reaction on the working model.

    Uses ``with model:`` context — does NOT modify the working model.
    """
    try:
        model = _require_model()
    except RuntimeError as exc:
        return str(exc)

    if target_reaction not in model.reactions:
        return (
            f"Error: reaction '{target_reaction}' not found in model. "
            f"Did you forget to add it via fba(action='add_pathway')?"
        )

    try:
        biomass_rxn_id = _get_biomass_reaction_id(model)
    except ValueError as exc:
        return f"Error: {exc}"

    try:
        with model:
            # If min_biomass_fraction > 0, constrain biomass
            if min_biomass_fraction > 0:
                wt_solution = model.optimize()
                if wt_solution.status != "optimal":
                    return "Error: model is infeasible; cannot compute biomass constraint."
                wt_biomass = wt_solution.objective_value
                biomass_rxn = model.reactions.get_by_id(biomass_rxn_id)
                biomass_rxn.lower_bound = min_biomass_fraction * wt_biomass

            # Set objective to target reaction
            model.objective = target_reaction
            solution = model.optimize()

            if solution.status != "optimal":
                return (
                    f"Model is infeasible when maximising '{target_reaction}' "
                    f"(min_biomass_fraction={min_biomass_fraction})."
                )

            max_product_flux = solution.objective_value
            biomass_flux = solution.fluxes[biomass_rxn_id]

            # Top 10 flux-carrying reactions by absolute flux
            flux_series = solution.fluxes
            top_rxns = flux_series.abs().sort_values(ascending=False).head(10)

            top_lines = []
            for rid in top_rxns.index:
                top_lines.append(f"  {rid:30s}  flux={flux_series[rid]:.6f}")

            return (
                f"Model: {_WORKING_MODEL_NAME}\n"
                f"Target reaction: {target_reaction}\n"
                f"Min biomass fraction: {min_biomass_fraction}\n"
                f"Max product flux: {max_product_flux:.6f}\n"
                f"Biomass flux: {biomass_flux:.6f}\n"
                f"\nTop 10 flux-carrying reactions:\n"
                + "\n".join(top_lines)
            )
    except Exception as exc:
        return f"Error during product maximisation: {exc}"


# ── 7. envelope ─────────────────────────────────────────────────────────

def fba_envelope(target_reaction: str, steps: int = 10) -> str:
    """Production envelope on the working model.

    Uses ``with model:`` context — does NOT modify the working model.
    """
    try:
        model = _require_model()
    except RuntimeError as exc:
        return str(exc)

    if target_reaction not in model.reactions:
        return (
            f"Error: reaction '{target_reaction}' not found in model. "
            f"Did you forget to add it via fba(action='add_pathway')?"
        )

    if steps < 1:
        return "Error: steps must be >= 1."

    try:
        biomass_rxn_id = _get_biomass_reaction_id(model)
    except ValueError as exc:
        return f"Error: {exc}"

    try:
        with model:
            # Get wild-type biomass
            wt_solution = model.optimize()
            if wt_solution.status != "optimal":
                return "Error: model is infeasible; cannot compute production envelope."
            wt_biomass = wt_solution.objective_value

            biomass_rxn = model.reactions.get_by_id(biomass_rxn_id)

            # Switch objective to target
            model.objective = target_reaction

            # Sweep
            rows = []
            for i in range(steps + 1):
                fraction = i / steps
                biomass_min = fraction * wt_biomass
                biomass_rxn.lower_bound = biomass_min

                sol = model.optimize()
                if sol.status == "optimal":
                    max_prod = sol.objective_value
                else:
                    max_prod = float("nan")

                rows.append((biomass_min, max_prod))

            # Format as text table
            lines = [
                f"Model: {_WORKING_MODEL_NAME}",
                f"Target reaction: {target_reaction}",
                f"Wild-type biomass: {wt_biomass:.6f}",
                f"Steps: {steps}",
                "",
                f"{'biomass_min':>15s} | {'max_product_flux':>18s}",
                f"{'-' * 15}-+-{'-' * 18}",
            ]
            for bm, prod in rows:
                prod_str = f"{prod:.6f}" if prod == prod else "infeasible"
                lines.append(f"{bm:>15.6f} | {prod_str:>18s}")

            return "\n".join(lines)
    except Exception as exc:
        return f"Error during production envelope computation: {exc}"


# ═══════════════════════════════════════════════════════════════════════════
# GEM query functions (independent of working model)
# ═══════════════════════════════════════════════════════════════════════════

# ── Query GPR ────────────────────────────────────────────────────────────

def query_gpr(model_name: str, gene_id: str) -> str:
    """Look up the Gene-Protein-Reaction associations for a gene."""
    try:
        model = _load_model(model_name)
    except Exception as exc:
        return f"Error loading model: {exc}"

    model_gene_ids = {g.id for g in model.genes}
    if gene_id not in model_gene_ids:
        return (
            f"Error: gene '{gene_id}' not found in model. "
            f"Available genes (first 20): "
            f"{sorted(model_gene_ids)[:20]} ..."
        )

    gene = model.genes.get_by_id(gene_id)
    associated_rxns = list(gene.reactions)

    if not associated_rxns:
        return f"Gene '{gene_id}' ({gene.name}) has no associated reactions."

    lines = [
        f"Gene: {gene_id}  Name: {gene.name or 'N/A'}",
        f"Associated reactions ({len(associated_rxns)}):",
    ]
    for r in sorted(associated_rxns, key=lambda x: x.id):
        lines.append(f"  Reaction: {r.id}")
        lines.append(f"    Name: {r.name}")
        lines.append(f"    GPR:  {r.gene_reaction_rule}")
        lines.append(f"    Stoichiometry: {r.reaction}")
        lines.append(f"    Bounds: [{r.lower_bound}, {r.upper_bound}]")
    return "\n".join(lines)


# ── Search metabolites and reactions ─────────────────────────────────────

def search_model(model_name: str, keyword: str) -> str:
    """Search metabolites and reactions by keyword (case-insensitive substring match)."""
    try:
        model = _load_model(model_name)
    except Exception as exc:
        return f"Error loading model: {exc}"

    if not keyword or not keyword.strip():
        return "Error: keyword must be a non-empty string."

    keyword_lower = keyword.strip().lower()
    max_results = 30

    try:
        # Search metabolites
        matched_mets = []
        for met in model.metabolites:
            if keyword_lower in met.id.lower() or keyword_lower in (met.name or "").lower():
                matched_mets.append(met)
                if len(matched_mets) >= max_results:
                    break

        # Search reactions
        matched_rxns = []
        for rxn in model.reactions:
            if keyword_lower in rxn.id.lower() or keyword_lower in (rxn.name or "").lower():
                matched_rxns.append(rxn)
                if len(matched_rxns) >= max_results:
                    break

        lines = [f"Model: {model_name}", f"Search keyword: '{keyword.strip()}'", ""]

        lines.append(f"Matching metabolites ({len(matched_mets)}"
                     f"{'+' if len(matched_mets) >= max_results else ''}):")
        if matched_mets:
            for met in matched_mets:
                lines.append(
                    f"  {met.id:30s}  name={met.name or 'N/A':30s}  "
                    f"formula={met.formula or 'N/A':15s}  compartment={met.compartment or 'N/A'}"
                )
        else:
            lines.append("  (none)")

        lines.append("")
        lines.append(f"Matching reactions ({len(matched_rxns)}"
                     f"{'+' if len(matched_rxns) >= max_results else ''}):")
        if matched_rxns:
            for rxn in matched_rxns:
                lines.append(
                    f"  {rxn.id:30s}  name={rxn.name or 'N/A':30s}  "
                    f"reaction={rxn.reaction}"
                )
        else:
            lines.append("  (none)")

        return "\n".join(lines)
    except Exception as exc:
        return f"Error during search: {exc}"


# ── Metabolite → reactions (with coefficients, flux, GPR) ─────────────

def get_metabolite_reactions(model_name: str, metabolite_id: str) -> str:
    """Return every reaction that consumes or produces a specific metabolite.

    For each reaction the output includes the stoichiometric coefficient
    (negative = consumed, positive = produced), the wild-type baseline
    flux, and the GPR rule.  Transport reactions, other compartment
    variants, and exchange reactions are reported as well.
    """
    try:
        model = _load_model(model_name)
    except Exception as exc:
        return f"Error loading model: {exc}"

    if metabolite_id not in model.metabolites:
        return (
            f"Error: metabolite '{metabolite_id}' not found in model. "
            f"Use gem(action='search') first to find the correct ID."
        )

    met = model.metabolites.get_by_id(metabolite_id)

    # ── Run baseline FBA once ────────────────────────────────────────
    try:
        solution = model.optimize()
        has_flux = solution.status == "optimal"
    except Exception:
        has_flux = False

    # ── Classify reactions ───────────────────────────────────────────
    consuming: list[tuple] = []   # (rxn, coeff, flux)
    producing: list[tuple] = []

    for rxn in sorted(met.reactions, key=lambda r: r.id):
        coeff = rxn.metabolites[met]
        flux = solution.fluxes[rxn.id] if has_flux else float("nan")
        entry = (rxn, coeff, flux)
        if coeff < 0:
            consuming.append(entry)
        else:
            producing.append(entry)

    # Sort each group by absolute flux descending (highest competition first)
    if has_flux:
        consuming.sort(key=lambda e: abs(e[2]), reverse=True)
        producing.sort(key=lambda e: abs(e[2]), reverse=True)

    # ── Other compartments & transport & exchange ────────────────────
    base_id = met.id.rsplit("_", 1)[0]  # e.g. "pyr" from "pyr_c"
    other_compartments = []
    for m in model.metabolites:
        if m.id != met.id and m.id.rsplit("_", 1)[0] == base_id:
            other_compartments.append(m.id)

    exchange_rxns = [
        r for r in model.reactions
        if r.id.startswith("EX_") and any(
            m.id.rsplit("_", 1)[0] == base_id for m in r.metabolites
        )
    ]

    # ── Format output ────────────────────────────────────────────────
    def _fmt_rxn(rxn, coeff, flux):
        gpr = rxn.gene_reaction_rule or "(no GPR)"
        flux_str = f"{flux:+.4f}" if flux == flux else "N/A"
        return (
            f"  {rxn.id:25s}  coeff={coeff:+.1f}  "
            f"flux={flux_str:>10s}  GPR: {gpr}\n"
            f"  {'':25s}  {rxn.reaction}"
        )

    lines = [
        f"Metabolite: {met.id}  ({met.name or 'N/A'})  "
        f"[compartment: {met.compartment}]",
        "",
    ]

    lines.append(f"Consuming reactions ({len(consuming)}):")
    if consuming:
        for rxn, coeff, flux in consuming:
            lines.append(_fmt_rxn(rxn, coeff, flux))
    else:
        lines.append("  (none)")

    lines.append("")
    lines.append(f"Producing reactions ({len(producing)}):")
    if producing:
        for rxn, coeff, flux in producing:
            lines.append(_fmt_rxn(rxn, coeff, flux))
    else:
        lines.append("  (none)")

    lines.append("")
    lines.append(
        f"Other compartments: "
        f"{', '.join(sorted(other_compartments)) if other_compartments else '(none)'}"
    )
    if exchange_rxns:
        for er in exchange_rxns:
            lines.append(
                f"Exchange: {er.id}  bounds=[{er.lower_bound}, {er.upper_bound}]"
            )
    else:
        lines.append("Exchange: (none)")

    return "\n".join(lines)
