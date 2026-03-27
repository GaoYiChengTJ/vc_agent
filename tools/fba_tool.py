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


# ── 4. Simulate overexpression ──────────────────────────────────────────────

def simulate_overexpression(
    model_name: str,
    target_reaction: str,
    forced_lower_bound: float,
) -> str:
    """Force a minimum flux on a reaction and observe biomass impact."""
    try:
        model = _load_model(model_name)
    except Exception as exc:
        return f"Error loading model: {exc}"

    if target_reaction not in model.reactions:
        available = ", ".join(sorted(r.id for r in model.reactions))
        return (
            f"Error: reaction '{target_reaction}' not found. "
            f"Available reactions: {available}"
        )

    try:
        baseline_sol = model.optimize()
        if baseline_sol.status != "optimal":
            return "Error: baseline model is infeasible."
        baseline_biomass = baseline_sol.objective_value
        baseline_target_flux = baseline_sol.fluxes[target_reaction]

        with model:
            rxn = model.reactions.get_by_id(target_reaction)
            rxn.lower_bound = forced_lower_bound

            solution = model.optimize()

            if solution.status != "optimal":
                return (
                    f"Model is infeasible when forcing '{target_reaction}' "
                    f"lower_bound={forced_lower_bound}. The forced flux "
                    f"exceeds the metabolic network capacity."
                )

            new_biomass = solution.objective_value
            new_target_flux = solution.fluxes[target_reaction]
            biomass_change = new_biomass - baseline_biomass
            pct_change = (
                (biomass_change / baseline_biomass * 100)
                if baseline_biomass != 0
                else 0.0
            )

            return (
                f"Model: {model_name}\n"
                f"Overexpression of '{target_reaction}' "
                f"(forced lower_bound={forced_lower_bound}):\n"
                f"Baseline biomass: {baseline_biomass:.6f}\n"
                f"New biomass:      {new_biomass:.6f}  "
                f"({pct_change:+.2f}%)\n"
                f"Baseline '{target_reaction}' flux: "
                f"{baseline_target_flux:.6f}\n"
                f"New '{target_reaction}' flux:      "
                f"{new_target_flux:.6f}"
            )
    except Exception as exc:
        return f"Error during overexpression simulation: {exc}"
