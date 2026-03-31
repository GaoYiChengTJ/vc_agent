"""
MCP Server for Metabolic Engineering Tools

7 coarse-grained tools:

  kegg          — KEGG search/compound/pathway/reaction/orthology (5 actions)
  gem           — GEM model search, GPR lookup, metabolite reactions (3 actions)
  uniprot       — UniProt/InterPro search and entry retrieval (3 actions)
  protein       — Protein analysis: annotate, interactions, kinetics, structure (4 actions)
  fba           — Stateful FBA: reset, add_pathway, knockout, overexpress,
                  media, maximize, envelope (7 actions)
  pubmed_search — PubMed literature search
  dna_optimize  — Codon optimization for S. cerevisiae
"""

import sys
from pathlib import Path
from typing import Literal, Optional

# Ensure project root is importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from mcp.server.fastmcp import FastMCP

# ── Import tool functions ─────────────────────────────────────────────

from tools.db_tools import (
    search_kegg,
    query_kegg_compound,
    query_kegg_pathway,
    query_kegg_reaction,
    query_kegg_orthology,
    search_pubmed,
    uniprot_search,
    uniprot_entry,
    interpro_entry,
)
from tools.fba_tool import (
    # GEM queries (stateless)
    search_model,
    query_gpr,
    get_metabolite_reactions,
    # FBA simulation (stateful)
    fba_reset,
    fba_add_pathway,
    fba_knockout,
    fba_overexpress,
    fba_media,
    fba_maximize,
    fba_envelope,
)
from tools.dna_tool import optimize_sequence
from tools.protein_tools import (
    protein_annotation,
    protein_interactions,
    enzyme_params,
    protein_structure,
)

# ── Create MCP server ─────────────────────────────────────────────────────

mcp = FastMCP("synbio-tools")


# ══════════════════════════════════════════════════════════════════════════
# Tool 1: KEGG (5 actions)
# ══════════════════════════════════════════════════════════════════════════

@mcp.tool()
def kegg(
    action: Literal["search", "compound", "pathway", "reaction", "orthology"],
    database: Optional[str] = None,
    keyword: Optional[str] = None,
    compound_id: Optional[str] = None,
    pathway_id: Optional[str] = None,
    reaction_id: Optional[str] = None,
    ko_id: Optional[str] = None,
) -> str:
    """Query KEGG databases for compounds, pathways, reactions, and orthology.

    Actions:
      action="search"     — Search a KEGG database by keyword.
                             Required: database (one of 'compound','reaction','enzyme','module','pathway'), keyword.
      action="compound"   — Get detailed info for a KEGG compound (e.g. 'C00186').
                             Required: compound_id.
      action="pathway"    — Get detailed info for a KEGG pathway (e.g. 'map00943') or module (e.g. 'M00941').
                             Required: pathway_id.
      action="reaction"   — Get detailed info for a KEGG reaction (e.g. 'R00703').
                             Required: reaction_id.
      action="orthology"  — Get a KEGG orthology entry (e.g. 'K00016').
                             Required: ko_id.
    """
    if action == "search":
        if not database or not keyword:
            return "Error: action='search' requires 'database' and 'keyword'."
        return search_kegg(database, keyword)
    elif action == "compound":
        if not compound_id:
            return "Error: action='compound' requires 'compound_id'."
        return query_kegg_compound(compound_id)
    elif action == "pathway":
        if not pathway_id:
            return "Error: action='pathway' requires 'pathway_id'."
        return query_kegg_pathway(pathway_id)
    elif action == "reaction":
        if not reaction_id:
            return "Error: action='reaction' requires 'reaction_id'."
        return query_kegg_reaction(reaction_id)
    elif action == "orthology":
        if not ko_id:
            return "Error: action='orthology' requires 'ko_id'."
        return query_kegg_orthology(ko_id)
    else:
        return f"Error: unknown action '{action}'."


# ══════════════════════════════════════════════════════════════════════════
# Tool 2: GEM (3 actions)
# ══════════════════════════════════════════════════════════════════════════

@mcp.tool()
def gem(
    action: Literal["search", "gpr", "reactions"],
    model_name: str = "iMM904",
    keyword: Optional[str] = None,
    gene_id: Optional[str] = None,
    metabolite_id: Optional[str] = None,
) -> str:
    """Query a genome-scale metabolic model (GEM).

    Available models: 'iMM904' (S. cerevisiae), 'e_coli_core' (E. coli).

    Actions:
      action="search"    — Search for metabolites and reactions by keyword (substring match).
                            Required: keyword.
      action="gpr"       — Look up Gene-Protein-Reaction associations for a gene.
                            Required: gene_id.
      action="reactions" — List ALL reactions that consume or produce a specific metabolite.
                            Returns stoichiometric coefficients, baseline flux, GPR,
                            other compartment variants, and exchange reactions.
                            Required: metabolite_id (exact model ID, e.g. 'pyr_c').
    """
    if action == "search":
        if not keyword:
            return "Error: action='search' requires 'keyword'."
        return search_model(model_name, keyword)
    elif action == "gpr":
        if not gene_id:
            return "Error: action='gpr' requires 'gene_id'."
        return query_gpr(model_name, gene_id)
    elif action == "reactions":
        if not metabolite_id:
            return "Error: action='reactions' requires 'metabolite_id' (exact model ID, e.g. 'pyr_c')."
        return get_metabolite_reactions(model_name, metabolite_id)
    else:
        return f"Error: unknown action '{action}'."


# ══════════════════════════════════════════════════════════════════════════
# Tool 3: UniProt / InterPro (3 actions)
# ══════════════════════════════════════════════════════════════════════════

@mcp.tool()
def uniprot(
    action: Literal["search", "entry", "interpro"],
    query: Optional[str] = None,
    max_results: int = 15,
    accession: Optional[str] = None,
    entry_id: Optional[str] = None,
) -> str:
    """Search UniProt for enzymes/transporters and retrieve protein or domain details.

    IMPORTANT: transporters MUST be searched here — KEGG does not index them.

    Actions:
      action="search"   — Search UniProt. Uses UniProt query syntax,
                           e.g. 'ec:1.1.1.27 AND reviewed:true'.
                           Required: query. Optional: max_results (default 15).
      action="entry"    — Get detailed info for a UniProt accession (e.g. 'Q9JXQ6').
                           Returns function, kinetics, cofactors, PDB structures.
                           Required: accession.
      action="interpro" — Get protein family/domain info from InterPro (e.g. 'IPR000292')
                           or Pfam (e.g. 'PF01226').
                           Required: entry_id.
    """
    if action == "search":
        if not query:
            return "Error: action='search' requires 'query'."
        return uniprot_search(query, max_results)
    elif action == "entry":
        if not accession:
            return "Error: action='entry' requires 'accession'."
        return uniprot_entry(accession)
    elif action == "interpro":
        if not entry_id:
            return "Error: action='interpro' requires 'entry_id'."
        return interpro_entry(entry_id)
    else:
        return f"Error: unknown action '{action}'."


# ══════════════════════════════════════════════════════════════════════════
# Tool 4: Protein analysis (4 actions)
# ══════════════════════════════════════════════════════════════════════════

@mcp.tool()
def protein(
    action: Literal["annotate", "interactions", "kinetics", "structure"],
    accession: Optional[str] = None,
    protein_id: Optional[str] = None,
    organism_taxid: Optional[int] = None,
    ec_number: Optional[str] = None,
    organism: Optional[str] = None,
) -> str:
    """Deep protein-level analysis for enzyme optimisation.

    Actions:
      action="annotate"      — Structural/functional features: domains, TM helices,
                                signal peptide, localization, cofactors, MW, mutagenesis data.
                                Required: accession (UniProt ID, e.g. 'Q9SWR5').
      action="interactions"  — Protein-protein interactions from STRING DB.
                                Required: protein_id (gene name or UniProt ID).
                                Optional: organism_taxid (default 4932 for S. cerevisiae).
      action="kinetics"      — Enzyme kinetic parameters (Km, kcat, Vmax, Ki) and known
                                mutations with measured effects from SABIO-RK / UniProt.
                                Required: ec_number (e.g. '1.14.14.87').
                                Optional: organism (e.g. 'Glycine max').
      action="structure"     — 3D structure from PDB (experimental) and AlphaFold (predicted).
                                Active site, binding site, and metal binding residue annotation.
                                Required: accession (UniProt ID).
    """
    if action == "annotate":
        if not accession:
            return "Error: action='annotate' requires 'accession' (UniProt ID)."
        return protein_annotation(accession)
    elif action == "interactions":
        pid = protein_id or accession
        if not pid:
            return "Error: action='interactions' requires 'protein_id' or 'accession'."
        return protein_interactions(pid, organism_taxid or 4932)
    elif action == "kinetics":
        if not ec_number:
            return "Error: action='kinetics' requires 'ec_number' (e.g. '1.1.1.27')."
        return enzyme_params(ec_number, organism)
    elif action == "structure":
        if not accession:
            return "Error: action='structure' requires 'accession' (UniProt ID)."
        return protein_structure(accession)
    else:
        return f"Error: unknown action '{action}'."


# ══════════════════════════════════════════════════════════════════════════
# Tool 5: FBA — stateful (7 actions)
# ══════════════════════════════════════════════════════════════════════════

@mcp.tool()
def fba(
    action: Literal[
        "reset", "add_pathway", "knockout", "overexpress",
        "media", "maximize", "envelope",
    ],
    model_name: str = "iMM904",
    reactions: Optional[list[dict]] = None,
    knockouts: Optional[list[str]] = None,
    gene_id: Optional[str] = None,
    min_flux: Optional[float] = None,
    carbon_source: Optional[str] = None,
    oxygen_lower_bound: Optional[float] = None,
    target_reaction: Optional[str] = None,
    min_biomass_fraction: float = 0.0,
    steps: int = 10,
) -> str:
    """Stateful Flux Balance Analysis on a genome-scale metabolic model.

    The working model persists across calls. Modification actions change it
    in place; analysis actions read it without modifying.

    WORKFLOW: reset → add_pathway → knockout → overexpress → media → maximize/envelope.
    Call reset to start a new scenario.

    Modification actions (change working model):
      action="reset"        — Load a fresh model. Required: model_name.
      action="add_pathway"  — Add heterologous reactions. Required: reactions (list of {id, reaction_string}).
      action="knockout"     — Knock out genes. Required: knockouts (list of gene IDs).
      action="overexpress"  — Force min flux on gene-associated reactions.
                               Required: gene_id, min_flux.
      action="media"        — Change carbon source and oxygen level.
                               Required: carbon_source, oxygen_lower_bound.
                               oxygen_lower_bound: 0=anaerobic, -0.5=microaerobic, -1000=aerobic.

    Analysis actions (read-only, do not modify working model):
      action="maximize"     — Maximize flux through a target reaction.
                               Required: target_reaction. Optional: min_biomass_fraction (default 0.0).
      action="envelope"     — Production envelope: sweep biomass 0-100%, compute max product flux.
                               Required: target_reaction. Optional: steps (default 10).
    """
    if action == "reset":
        return fba_reset(model_name)

    elif action == "add_pathway":
        if not reactions:
            return "Error: action='add_pathway' requires 'reactions' (list of {id, reaction_string})."
        return fba_add_pathway(reactions)

    elif action == "knockout":
        if not knockouts:
            return "Error: action='knockout' requires 'knockouts' (list of gene IDs)."
        return fba_knockout(knockouts)

    elif action == "overexpress":
        if not gene_id or min_flux is None:
            return "Error: action='overexpress' requires 'gene_id' and 'min_flux'."
        return fba_overexpress(gene_id, min_flux)

    elif action == "media":
        if not carbon_source or oxygen_lower_bound is None:
            return "Error: action='media' requires 'carbon_source' and 'oxygen_lower_bound'."
        return fba_media(carbon_source, oxygen_lower_bound)

    elif action == "maximize":
        if not target_reaction:
            return "Error: action='maximize' requires 'target_reaction'."
        return fba_maximize(target_reaction, min_biomass_fraction)

    elif action == "envelope":
        if not target_reaction:
            return "Error: action='envelope' requires 'target_reaction'."
        return fba_envelope(target_reaction, steps)

    else:
        return f"Error: unknown action '{action}'."


# ══════════════════════════════════════════════════════════════════════════
# Tool 6 & 7: Auxiliary (kept as-is)
# ══════════════════════════════════════════════════════════════════════════

@mcp.tool()
def pubmed_search(query: str, max_results: int = 10) -> str:
    """Search PubMed for relevant literature.
    Returns titles, authors, and abstracts."""
    return search_pubmed(query, max_results)


@mcp.tool()
def dna_optimize(sequence: str) -> str:
    """Codon-optimize a DNA coding sequence for S. cerevisiae.
    Maximizes CAI and removes BsaI sites while preserving protein."""
    return optimize_sequence(sequence)


# ── Run ───────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    mcp.run(transport="stdio")
