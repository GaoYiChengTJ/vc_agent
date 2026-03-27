"""DNA codon-optimization tool."""

from dnachisel import (
    AvoidPattern,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceTranslation,
)


def optimize_sequence(sequence: str) -> str:
    """Codon-optimize a DNA coding sequence for the target host organism.

    Args:
        sequence: A DNA coding sequence (multiple of 3, A/T/G/C only).

    Returns:
        A human-readable report with the optimized sequence, or an error
        message if something goes wrong.
    """
    seq = sequence.strip().upper()

    # --- basic validation ---
    invalid_chars = set(seq) - {"A", "T", "G", "C"}
    if invalid_chars:
        return (
            f"Error: sequence contains invalid characters: {invalid_chars}. "
            f"Only A, T, G, C are allowed."
        )

    if len(seq) == 0:
        return "Error: empty sequence provided."

    if len(seq) % 3 != 0:
        return (
            f"Error: sequence length ({len(seq)}) is not a multiple of 3. "
            f"Please provide a complete coding sequence."
        )

    try:
        problem = DnaOptimizationProblem(
            sequence=seq,
            constraints=[
                EnforceTranslation(),          # keep the protein identical
                AvoidPattern("GGTCTC"),         # BsaI forward site
                AvoidPattern("GAGACC"),         # BsaI reverse complement
            ],
            objectives=[
                CodonOptimize(species="s_cerevisiae"),
            ],
        )

        problem.resolve_constraints()
        problem.optimize()

        if not problem.all_constraints_pass():
            failed = [
                str(c) for c in problem.constraints if not c.evaluate(problem).passes
            ]
            return (
                f"Warning: optimisation finished but some constraints failed: "
                f"{failed}.\nOptimized sequence: {problem.sequence}"
            )

        return (
            f"Original sequence:  {seq}\n"
            f"Optimized sequence: {problem.sequence}\n"
            f"Length: {len(problem.sequence)} bp\n"
            f"All constraints passed: True\n"
            f"Constraints: translation preserved, BsaI sites (GGTCTC/GAGACC) avoided"
        )
    except Exception as exc:
        return f"Error during sequence optimisation: {exc}"
