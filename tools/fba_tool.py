"""FBA simulation tools for metabolic engineering."""

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


def _load_model(model_name: str):
    """Load a fresh copy of the specified metabolic model."""
    path = AVAILABLE_MODELS.get(model_name)
    if path is None:
        raise ValueError(
            f"Unknown model '{model_name}'. "
            f"Available: {', '.join(sorted(AVAILABLE_MODELS))}"
        )
    return cobra.io.read_sbml_model(str(path))


def _get_biomass_reaction_id(model) -> str:
    """Get the biomass reaction ID from the model's current objective."""
    for rxn in model.reactions:
        if rxn.objective_coefficient != 0:
            return rxn.id
    raise ValueError("No biomass/objective reaction found in model.")


# ── 1. Gene knockout ────────────────────────────────────────────────────────

def simulate_knockout(
    model_name: str,
    target_reaction: str,
    knockouts: list[str],
) -> str:
    """Knock out specified genes and report flux results."""
    try:
        model = _load_model(model_name)
    except Exception as exc:
        return f"Error loading model: {exc}"

    if target_reaction not in model.reactions:
        available = ", ".join(sorted(r.id for r in model.reactions))
        return (
            f"Error: reaction '{target_reaction}' not found in model. "
            f"Available reactions: {available}"
        )

    model_gene_ids = {g.id for g in model.genes}
    unknown = [g for g in knockouts if g not in model_gene_ids]
    if unknown:
        return (
            f"Error: gene(s) {unknown} not found in model. "
            f"Available genes (first 20): "
            f"{sorted(model_gene_ids)[:20]} ..."
        )

    try:
        with model:
            cobra.manipulation.knock_out_model_genes(model, knockouts)
            solution = model.optimize()

            if solution.status != "optimal":
                return (
                    f"Model status after knockout: {solution.status}. "
                    f"The model is infeasible — the organism cannot grow "
                    f"with gene(s) {knockouts} knocked out."
                )

            biomass_flux = solution.objective_value
            target_flux = solution.fluxes[target_reaction]

            return (
                f"Model: {model_name}\n"
                f"Knockout: {knockouts}\n"
                f"Biomass flux: {biomass_flux:.6f}\n"
                f"Target reaction '{target_reaction}' flux: {target_flux:.6f}"
            )
    except Exception as exc:
        return f"Error during FBA simulation: {exc}"


# ── 2. Change media / environment ───────────────────────────────────────────

def change_media_and_simulate(
    model_name: str,
    carbon_source: str,
    oxygen_lower_bound: float,
    target_reaction: str,
) -> str:
    """Modify exchange reaction bounds and run FBA."""
    try:
        model = _load_model(model_name)
    except Exception as exc:
        return f"Error loading model: {exc}"

    ex_ids = sorted(r.id for r in model.reactions if r.id.startswith("EX_"))
    if carbon_source not in ex_ids:
        return (
            f"Error: '{carbon_source}' is not a valid exchange reaction. "
            f"Available exchange reactions: {', '.join(ex_ids)}"
        )

    if target_reaction not in model.reactions:
        available = ", ".join(sorted(r.id for r in model.reactions))
        return (
            f"Error: reaction '{target_reaction}' not found. "
            f"Available reactions: {available}"
        )

    try:
        with model:
            for rxn in model.reactions:
                if rxn.id.startswith("EX_") and rxn.id != carbon_source:
                    if rxn.lower_bound < 0 and rxn.id not in _INORGANIC_EXCHANGES:
                        rxn.lower_bound = 0.0

            model.reactions.get_by_id(carbon_source).lower_bound = -10.0
            model.reactions.get_by_id("EX_o2_e").lower_bound = oxygen_lower_bound

            solution = model.optimize()

            if solution.status != "optimal":
                return (
                    f"Model is infeasible under the specified media "
                    f"(carbon={carbon_source}, O2 lb={oxygen_lower_bound})."
                )

            biomass_flux = solution.objective_value
            target_flux = solution.fluxes[target_reaction]
            o2_flux = solution.fluxes["EX_o2_e"]

            if oxygen_lower_bound == 0:
                condition = "anaerobic"
            elif oxygen_lower_bound > -1000:
                condition = f"micro-aerobic (O2 lb={oxygen_lower_bound})"
            else:
                condition = "aerobic"
            return (
                f"Model: {model_name}\n"
                f"Media condition: {condition}, carbon source={carbon_source}\n"
                f"O2 exchange flux: {o2_flux:.6f}\n"
                f"Biomass flux: {biomass_flux:.6f}\n"
                f"Target reaction '{target_reaction}' flux: {target_flux:.6f}"
            )
    except Exception as exc:
        return f"Error during media simulation: {exc}"


# ── 3. Add heterologous reaction ────────────────────────────────────────────

def add_heterologous_reaction(
    model_name: str,
    reaction_id: str,
    reaction_string: str,
) -> str:
    """Add a new reaction to the model and maximise its flux."""
    try:
        model = _load_model(model_name)
    except Exception as exc:
        return f"Error loading model: {exc}"

    if reaction_id in model.reactions:
        return f"Error: reaction '{reaction_id}' already exists in the model."

    try:
        biomass_rxn_id = _get_biomass_reaction_id(model)

        new_rxn = cobra.Reaction(reaction_id)
        model.add_reactions([new_rxn])
        new_rxn.build_reaction_from_string(reaction_string)
    except Exception as exc:
        return f"Error building reaction from string: {exc}"

    try:
        baseline_sol = model.optimize()
        baseline_biomass = (
            baseline_sol.objective_value
            if baseline_sol.status == "optimal"
            else 0.0
        )

        model.objective = reaction_id
        solution = model.optimize()

        if solution.status != "optimal":
            return (
                f"Reaction '{reaction_id}' ({reaction_string}) was added, "
                f"but FBA is infeasible when maximising it. "
                f"The substrates may not be producible by the model."
            )

        max_flux = solution.objective_value
        biomass_under_max = solution.fluxes[biomass_rxn_id]

        return (
            f"Model: {model_name}\n"
            f"Heterologous reaction added: {reaction_id}\n"
            f"Reaction: {reaction_string}\n"
            f"Max theoretical flux: {max_flux:.6f}\n"
            f"Baseline biomass (before objective switch): {baseline_biomass:.6f}\n"
            f"Biomass when maximising product: {biomass_under_max:.6f}"
        )
    except Exception as exc:
        return f"Error during heterologous reaction FBA: {exc}"


# ── 4. Simulate overexpression (gene-level, GPR-aware) ──────────────────────

def simulate_overexpression(
    model_name: str,
    gene_id: str,
    forced_lower_bound: float,
) -> str:
    """Look up reactions via GPR, force minimum flux, observe biomass impact."""
    try:
        model = _load_model(model_name)
    except Exception as exc:
        return f"Error loading model: {exc}"

    # Validate gene exists
    model_gene_ids = {g.id for g in model.genes}
    if gene_id not in model_gene_ids:
        return (
            f"Error: gene '{gene_id}' not found in model. "
            f"Available genes (first 20): "
            f"{sorted(model_gene_ids)[:20]} ..."
        )

    # GPR lookup: find reactions associated with this gene
    gene = model.genes.get_by_id(gene_id)
    associated_rxns = list(gene.reactions)
    if not associated_rxns:
        return f"Gene '{gene_id}' has no associated reactions in the model."

    gpr_info = "\n".join(
        f"  {r.id:20s}  GPR: {r.gene_reaction_rule}"
        for r in associated_rxns
    )

    try:
        # Baseline
        baseline_sol = model.optimize()
        if baseline_sol.status != "optimal":
            return "Error: baseline model is infeasible."
        baseline_biomass = baseline_sol.objective_value

        # Overexpression: force minimum flux on all associated reactions
        with model:
            per_rxn_results = []
            for rxn in associated_rxns:
                old_lb = rxn.lower_bound
                rxn.lower_bound = max(rxn.lower_bound, forced_lower_bound)
                per_rxn_results.append(
                    (rxn.id, old_lb, rxn.lower_bound)
                )

            solution = model.optimize()

            if solution.status != "optimal":
                bounds_desc = ", ".join(
                    f"{rid} lb: {old}->{new}"
                    for rid, old, new in per_rxn_results
                )
                return (
                    f"Model is infeasible when overexpressing gene '{gene_id}'.\n"
                    f"Forced bounds: {bounds_desc}\n"
                    f"The forced flux exceeds the metabolic network capacity."
                )

            new_biomass = solution.objective_value
            pct_change = (
                (new_biomass - baseline_biomass) / baseline_biomass * 100
                if baseline_biomass != 0
                else 0.0
            )

            flux_details = "\n".join(
                f"  {rxn.id:20s}  baseline={baseline_sol.fluxes[rxn.id]:.6f}"
                f"  new={solution.fluxes[rxn.id]:.6f}"
                for rxn in associated_rxns
            )

            return (
                f"Model: {model_name}\n"
                f"Overexpression of gene '{gene_id}' "
                f"(forced lower_bound={forced_lower_bound}):\n"
                f"\nGPR associations:\n{gpr_info}\n"
                f"\nBiomass impact:\n"
                f"  Baseline: {baseline_biomass:.6f}\n"
                f"  New:      {new_biomass:.6f}  ({pct_change:+.2f}%)\n"
                f"\nFlux changes per reaction:\n{flux_details}"
            )
    except Exception as exc:
        return f"Error during overexpression simulation: {exc}"


# ── 5. Query GPR ────────────────────────────────────────────────────────────

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


# ── 6. Search metabolites and reactions ────────────────────────────────────

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


# ── 7. Add pathway (multiple reactions) ───────────────────────────────────

def add_pathway(model_name: str, reactions: list[dict]) -> str:
    """Add multiple reactions to the model and verify feasibility with FBA."""
    try:
        model = _load_model(model_name)
    except Exception as exc:
        return f"Error loading model: {exc}"

    if not reactions:
        return "Error: reactions list must be non-empty."

    # Validate reaction dicts
    for i, rxn_dict in enumerate(reactions):
        if not isinstance(rxn_dict, dict):
            return f"Error: reactions[{i}] is not a dict."
        if "id" not in rxn_dict or "reaction_string" not in rxn_dict:
            return (
                f"Error: reactions[{i}] must have keys 'id' and 'reaction_string'. "
                f"Got keys: {list(rxn_dict.keys())}"
            )
        if rxn_dict["id"] in model.reactions:
            return f"Error: reaction '{rxn_dict['id']}' already exists in the model."

    try:
        biomass_rxn_id = _get_biomass_reaction_id(model)
    except Exception as exc:
        return f"Error identifying biomass reaction: {exc}"

    # Add all reactions
    added_ids = []
    try:
        for rxn_dict in reactions:
            new_rxn = cobra.Reaction(rxn_dict["id"])
            model.add_reactions([new_rxn])
            new_rxn.build_reaction_from_string(rxn_dict["reaction_string"])
            added_ids.append(rxn_dict["id"])
    except Exception as exc:
        return (
            f"Error adding reaction '{rxn_dict['id']}': {exc}\n"
            f"Successfully added before failure: {added_ids}"
        )

    # Run FBA with default objective (biomass) to verify feasibility
    try:
        solution = model.optimize()

        if solution.status != "optimal":
            return (
                f"Reactions added: {added_ids}\n"
                f"WARNING: model is infeasible after adding the pathway "
                f"(status: {solution.status})."
            )

        biomass_flux = solution.objective_value
        flux_lines = []
        for rid in added_ids:
            flux_lines.append(f"  {rid:30s}  flux={solution.fluxes[rid]:.6f}")

        return (
            f"Model: {model_name}\n"
            f"Added reactions ({len(added_ids)}):\n"
            + "\n".join(
                f"  {rd['id']:30s}  {rd['reaction_string']}"
                for rd in reactions
            )
            + f"\n\nBiomass flux after addition: {biomass_flux:.6f}\n"
            f"\nFlux through new reactions (biomass-optimal):\n"
            + "\n".join(flux_lines)
        )
    except Exception as exc:
        return f"Error during FBA after adding pathway: {exc}"


# ── 8. Maximize product flux ──────────────────────────────────────────────

def maximize_product(
    model_name: str,
    target_reaction: str,
    min_biomass_fraction: float = 0.0,
    extra_reactions: list[dict] = None,
) -> str:
    """Maximize flux through a target reaction, optionally constraining minimum biomass."""
    try:
        model = _load_model(model_name)
    except Exception as exc:
        return f"Error loading model: {exc}"

    # Add extra reactions if provided
    if extra_reactions:
        for i, rxn_dict in enumerate(extra_reactions):
            if not isinstance(rxn_dict, dict):
                return f"Error: extra_reactions[{i}] is not a dict."
            if "id" not in rxn_dict or "reaction_string" not in rxn_dict:
                return (
                    f"Error: extra_reactions[{i}] must have keys 'id' and 'reaction_string'. "
                    f"Got keys: {list(rxn_dict.keys())}"
                )
        try:
            for rxn_dict in extra_reactions:
                if rxn_dict["id"] not in model.reactions:
                    new_rxn = cobra.Reaction(rxn_dict["id"])
                    model.add_reactions([new_rxn])
                    new_rxn.build_reaction_from_string(rxn_dict["reaction_string"])
        except Exception as exc:
            return f"Error adding extra reaction '{rxn_dict['id']}': {exc}"

    # Validate target reaction
    if target_reaction not in model.reactions:
        available = ", ".join(sorted(r.id for r in model.reactions))
        return (
            f"Error: reaction '{target_reaction}' not found in model. "
            f"Available reactions: {available}"
        )

    try:
        biomass_rxn_id = _get_biomass_reaction_id(model)
    except Exception as exc:
        return f"Error identifying biomass reaction: {exc}"

    try:
        # If min_biomass_fraction > 0, compute wild-type biomass and constrain
        if min_biomass_fraction > 0:
            wt_solution = model.optimize()
            if wt_solution.status != "optimal":
                return "Error: wild-type model is infeasible; cannot compute biomass constraint."
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
            f"Model: {model_name}\n"
            f"Target reaction: {target_reaction}\n"
            f"Min biomass fraction: {min_biomass_fraction}\n"
            f"Max product flux: {max_product_flux:.6f}\n"
            f"Biomass flux: {biomass_flux:.6f}\n"
            f"\nTop 10 flux-carrying reactions:\n"
            + "\n".join(top_lines)
        )
    except Exception as exc:
        return f"Error during product maximisation: {exc}"


# ── 9. Production envelope ────────────────────────────────────────────────

def production_envelope(
    model_name: str,
    target_reaction: str,
    steps: int = 10,
    extra_reactions: list[dict] = None,
) -> str:
    """Sweep biomass constraint and compute max product flux at each level."""
    try:
        model = _load_model(model_name)
    except Exception as exc:
        return f"Error loading model: {exc}"

    # Add extra reactions if provided
    if extra_reactions:
        for i, rxn_dict in enumerate(extra_reactions):
            if not isinstance(rxn_dict, dict):
                return f"Error: extra_reactions[{i}] is not a dict."
            if "id" not in rxn_dict or "reaction_string" not in rxn_dict:
                return (
                    f"Error: extra_reactions[{i}] must have keys 'id' and 'reaction_string'. "
                    f"Got keys: {list(rxn_dict.keys())}"
                )
        try:
            for rxn_dict in extra_reactions:
                if rxn_dict["id"] not in model.reactions:
                    new_rxn = cobra.Reaction(rxn_dict["id"])
                    model.add_reactions([new_rxn])
                    new_rxn.build_reaction_from_string(rxn_dict["reaction_string"])
        except Exception as exc:
            return f"Error adding extra reaction '{rxn_dict['id']}': {exc}"

    # Validate target reaction
    if target_reaction not in model.reactions:
        available = ", ".join(sorted(r.id for r in model.reactions))
        return (
            f"Error: reaction '{target_reaction}' not found in model. "
            f"Available reactions: {available}"
        )

    if steps < 1:
        return "Error: steps must be >= 1."

    try:
        biomass_rxn_id = _get_biomass_reaction_id(model)
    except Exception as exc:
        return f"Error identifying biomass reaction: {exc}"

    try:
        # Get wild-type biomass
        wt_solution = model.optimize()
        if wt_solution.status != "optimal":
            return "Error: wild-type model is infeasible; cannot compute production envelope."
        wt_biomass = wt_solution.objective_value

        biomass_rxn = model.reactions.get_by_id(biomass_rxn_id)
        original_lb = biomass_rxn.lower_bound

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

        # Restore original lower bound
        biomass_rxn.lower_bound = original_lb

        # Format as text table
        lines = [
            f"Model: {model_name}",
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
