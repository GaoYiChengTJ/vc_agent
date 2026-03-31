"""
Metabolic Engineering Agent — 4-Phase Framework

Built on Claude Agent SDK.  Registers all tools organized by the 4-phase
pipeline described in Framework Design.md:

  Phase 1  Target Identification   (KEGG)
  Phase 2  Host Assessment         (GEM / FBA)
  Phase 3  Enzyme Sourcing         (UniProt / InterPro)
  Phase 4  Simulation & Optim.     (FBA)
  Optional Sequence Engineering    (DNAchisel)
"""

import sys
from pathlib import Path

import os
import anyio
from dotenv import load_dotenv

load_dotenv(Path(__file__).resolve().parent.parent / ".env")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

API_ENV = {
    "ANTHROPIC_API_KEY": os.environ.get("ANTHROPIC_API_KEY", ""),
    "ANTHROPIC_BASE_URL": os.environ.get("ANTHROPIC_BASE_URL", ""),
}

from claude_agent_sdk import (
    AssistantMessage,
    ClaudeAgentOptions,
    ClaudeSDKClient,
    ResultMessage,
    TextBlock,
    create_sdk_mcp_server,
    tool,
)

# ── Import all tool functions ─────────────────────────────────────────────

# Phase 1: KEGG
from tools.db_tools import (
    search_kegg,
    query_kegg_compound,
    query_kegg_pathway,
    query_kegg_reaction,
    query_kegg_orthology,
    search_pubmed,
)

# Phase 2 & 4: FBA
from tools.fba_tool import (
    search_model,
    simulate_knockout,
    change_media_and_simulate,
    add_heterologous_reaction,
    simulate_overexpression,
    query_gpr,
    add_pathway,
    maximize_product,
    production_envelope,
)

# Phase 3: UniProt / InterPro
from tools.db_tools import (
    uniprot_search,
    uniprot_entry,
    interpro_entry,
)

# Optional: DNA
from tools.dna_tool import optimize_sequence

MODEL_NAME_DESC = (
    "Name of the metabolic model. "
    "Available: 'e_coli_core' (E. coli), 'iMM904' (S. cerevisiae)."
)

# ══════════════════════════════════════════════════════════════════════════
# Phase 1 Tools — Target Identification (KEGG)
# ══════════════════════════════════════════════════════════════════════════

@tool(
    "search_kegg",
    "[Phase 1] Search a KEGG database (compound / reaction / enzyme / "
    "module / pathway) by keyword.  Start here to find the target compound.",
    {"database": str, "keyword": str},
)
async def _search_kegg(args: dict):
    return {"content": [{"type": "text", "text": search_kegg(args["database"], args["keyword"])}]}


@tool(
    "query_kegg_compound",
    "[Phase 1] Get detailed info for a KEGG compound ID (e.g. 'C00186'). "
    "Returns pathways, reactions, enzymes, and database links.",
    {"compound_id": str},
)
async def _query_kegg_compound(args: dict):
    return {"content": [{"type": "text", "text": query_kegg_compound(args["compound_id"])}]}


@tool(
    "query_kegg_pathway",
    "[Phase 1] Get detailed info for a KEGG pathway (e.g. 'map00943') or "
    "module (e.g. 'M00941'). Returns reactions, orthology, and compounds.",
    {"pathway_id": str},
)
async def _query_kegg_pathway(args: dict):
    return {"content": [{"type": "text", "text": query_kegg_pathway(args["pathway_id"])}]}


@tool(
    "query_kegg_reaction",
    "[Phase 1] Get detailed info for a KEGG reaction (e.g. 'R00703'). "
    "Returns equation, EC number, and orthology (KO ID).",
    {"reaction_id": str},
)
async def _query_kegg_reaction(args: dict):
    return {"content": [{"type": "text", "text": query_kegg_reaction(args["reaction_id"])}]}


@tool(
    "query_kegg_orthology",
    "[Phase 1] Get a KEGG orthology entry (e.g. 'K00016'). Shows gene "
    "orthologs across organisms. Use after finding a KO ID from a reaction.",
    {"ko_id": str},
)
async def _query_kegg_orthology(args: dict):
    return {"content": [{"type": "text", "text": query_kegg_orthology(args["ko_id"])}]}


# ══════════════════════════════════════════════════════════════════════════
# Phase 2 Tools — Host Assessment (GEM)
# ══════════════════════════════════════════════════════════════════════════

@tool(
    "search_model",
    "[Phase 2] Search the host GEM for metabolites and reactions matching a "
    "keyword. Use to check which Phase 1 metabolites/reactions already exist.",
    {"model_name": str, "keyword": str},
)
async def _search_model(args: dict):
    return {"content": [{"type": "text", "text": search_model(args["model_name"], args["keyword"])}]}


@tool(
    "query_gpr",
    "[Phase 2/4] Look up Gene-Protein-Reaction associations for a gene in "
    "the host GEM. Use before knockout/overexpression to understand the "
    "gene's metabolic role.",
    {"model_name": str, "gene_id": str},
)
async def _query_gpr(args: dict):
    return {"content": [{"type": "text", "text": query_gpr(args["model_name"], args["gene_id"])}]}


# ══════════════════════════════════════════════════════════════════════════
# Phase 3 Tools — Enzyme Sourcing (UniProt / InterPro)
# ══════════════════════════════════════════════════════════════════════════

@tool(
    "uniprot_search",
    "[Phase 3] Search UniProt for enzymes/transporters. Uses UniProt query "
    "syntax, e.g. 'ec:1.1.1.27 AND reviewed:true' or "
    "'gene:lgtA AND organism_name:neisseria'. "
    "IMPORTANT: transporters MUST be searched here — KEGG does not index them.",
    {"query": str, "max_results": int},
)
async def _uniprot_search(args: dict):
    return {"content": [{"type": "text", "text": uniprot_search(args["query"], args.get("max_results", 15))}]}


@tool(
    "uniprot_entry",
    "[Phase 3] Get detailed info for a UniProt accession (e.g. 'Q9JXQ6'). "
    "Returns function, kinetics (Km/Vmax), cofactors, PDB structures, "
    "and cross-references.",
    {"accession": str},
)
async def _uniprot_entry(args: dict):
    return {"content": [{"type": "text", "text": uniprot_entry(args["accession"])}]}


@tool(
    "interpro_entry",
    "[Phase 3] Get protein family/domain info from InterPro (e.g. "
    "'IPR000292') or Pfam (e.g. 'PF01226'). Returns domain architecture, "
    "species distribution, and structure counts.",
    {"entry_id": str},
)
async def _interpro_entry(args: dict):
    return {"content": [{"type": "text", "text": interpro_entry(args["entry_id"])}]}


# ══════════════════════════════════════════════════════════════════════════
# Phase 4 Tools — Simulation & Optimization (FBA)
# ══════════════════════════════════════════════════════════════════════════

@tool(
    "add_heterologous_reaction",
    "[Phase 4] Add a single heterologous reaction to the GEM and compute "
    "its max theoretical flux. reaction_string uses COBRApy format, "
    "e.g. 'pyr_c + nadh_c + h_c --> lac__L_c + nad_c'.",
    {"model_name": str, "reaction_id": str, "reaction_string": str},
)
async def _add_heterologous_reaction(args: dict):
    return {"content": [{"type": "text", "text": add_heterologous_reaction(
        args["model_name"], args["reaction_id"], args["reaction_string"]
    )}]}


@tool(
    "add_pathway",
    "[Phase 4] Add multiple heterologous reactions at once and verify "
    "feasibility. Input: list of dicts with 'id' and 'reaction_string' keys.",
    {"model_name": str, "reactions": list},
)
async def _add_pathway(args: dict):
    return {"content": [{"type": "text", "text": add_pathway(args["model_name"], args["reactions"])}]}


@tool(
    "maximize_product",
    "[Phase 4] Maximize flux through a target reaction (e.g. an exchange "
    "reaction) with optional minimum biomass constraint. Use to find the "
    "theoretical production ceiling.",
    {"model_name": str, "target_reaction": str, "min_biomass_fraction": float},
)
async def _maximize_product(args: dict):
    extra = args.get("extra_reactions")
    return {"content": [{"type": "text", "text": maximize_product(
        args["model_name"], args["target_reaction"],
        args.get("min_biomass_fraction", 0.0), extra
    )}]}


@tool(
    "production_envelope",
    "[Phase 4] Sweep biomass constraint from 0% to 100% and compute max "
    "product flux at each level. Reveals the growth-production tradeoff.",
    {"model_name": str, "target_reaction": str, "steps": int},
)
async def _production_envelope(args: dict):
    extra = args.get("extra_reactions")
    return {"content": [{"type": "text", "text": production_envelope(
        args["model_name"], args["target_reaction"],
        args.get("steps", 10), extra
    )}]}


@tool(
    "simulate_knockout",
    "[Phase 4] Knock out one or more genes and report biomass and target "
    "reaction flux. Use to test elimination of competing pathways. "
    "Pass extra_reactions (list of {id, reaction_string} dicts) to test "
    "on a model with heterologous pathway added.",
    {"model_name": str, "target_reaction": str, "knockouts": list[str]},
)
async def _simulate_knockout(args: dict):
    return {"content": [{"type": "text", "text": simulate_knockout(
        args["model_name"], args["target_reaction"], args["knockouts"],
        args.get("extra_reactions"),
    )}]}


@tool(
    "simulate_overexpression",
    "[Phase 4] Force a minimum flux on all reactions associated with a gene "
    "(via GPR). Use to test precursor supply enhancement. "
    "Pass extra_reactions to test on a model with heterologous pathway. "
    "For endogenous products, include an export reaction in extra_reactions.",
    {"model_name": str, "gene_id": str, "forced_lower_bound": float},
)
async def _simulate_overexpression(args: dict):
    return {"content": [{"type": "text", "text": simulate_overexpression(
        args["model_name"], args["gene_id"], args["forced_lower_bound"],
        args.get("extra_reactions"),
    )}]}


@tool(
    "change_media_and_simulate",
    "[Phase 4] Change growth medium (carbon source, O₂ level) and run FBA. "
    "Use oxygen_lower_bound=0 for anaerobic, -1000 for fully aerobic. "
    "Pass extra_reactions to test on a model with heterologous pathway.",
    {"model_name": str, "carbon_source": str, "oxygen_lower_bound": float, "target_reaction": str},
)
async def _change_media(args: dict):
    return {"content": [{"type": "text", "text": change_media_and_simulate(
        args["model_name"], args["carbon_source"],
        args["oxygen_lower_bound"], args["target_reaction"],
        args.get("extra_reactions"),
    )}]}


# ══════════════════════════════════════════════════════════════════════════
# Optional — Sequence Engineering
# ══════════════════════════════════════════════════════════════════════════

@tool(
    "optimize_sequence",
    "[Optional] Codon-optimize a DNA coding sequence for S. cerevisiae. "
    "Maximizes CAI and removes BsaI sites while preserving protein.",
    {"sequence": str},
)
async def _optimize_sequence(args: dict):
    return {"content": [{"type": "text", "text": optimize_sequence(args["sequence"])}]}


# ══════════════════════════════════════════════════════════════════════════
# Auxiliary
# ══════════════════════════════════════════════════════════════════════════

@tool(
    "search_pubmed",
    "[Auxiliary] Search PubMed for relevant literature. Returns titles, "
    "authors, and abstracts.",
    {"query": str, "max_results": int},
)
async def _search_pubmed(args: dict):
    return {"content": [{"type": "text", "text": search_pubmed(args["query"], args.get("max_results", 10))}]}


# ══════════════════════════════════════════════════════════════════════════
# MCP Server — register all tools
# ══════════════════════════════════════════════════════════════════════════

mcp_server = create_sdk_mcp_server(
    "synbio-tools",
    tools=[
        # Phase 1: KEGG
        _search_kegg,
        _query_kegg_compound,
        _query_kegg_pathway,
        _query_kegg_reaction,
        _query_kegg_orthology,
        # Phase 2: Host assessment
        _search_model,
        _query_gpr,
        # Phase 3: Enzyme sourcing
        _uniprot_search,
        _uniprot_entry,
        _interpro_entry,
        # Phase 4: Simulation & optimization
        _add_heterologous_reaction,
        _add_pathway,
        _maximize_product,
        _production_envelope,
        _simulate_knockout,
        _simulate_overexpression,
        _change_media,
        # Optional
        _optimize_sequence,
        # Auxiliary
        _search_pubmed,
    ],
)

# ══════════════════════════════════════════════════════════════════════════
# System Prompt — 4-Phase Pipeline
# ══════════════════════════════════════════════════════════════════════════

SYSTEM_PROMPT = """\
You are a metabolic engineering design agent.  When the user provides a
(target product, host organism) pair, follow this 4-phase workflow to
systematically design a complete engineering strategy.

═══ Phase 1 · Target Identification (KEGG) ═══
Use search_kegg to find the target compound.  Then query_kegg_compound to
get pathways, reactions, and enzyme EC numbers.  Drill into each pathway/
module with query_kegg_pathway and each reaction with query_kegg_reaction
to extract the full biosynthetic route and KO IDs.
If KEGG has no direct compound entry, fall back to searching by pathway
or module keywords and infer the route from structural analogs.

═══ Phase 2 · Host Assessment (GEM) ═══
Use search_model to check which Phase 1 metabolites and reactions already
exist in the host genome-scale model (GEM).  Flag missing ones as needing
heterologous introduction.  Run baseline FBA to understand current flux
distribution.  Use query_gpr to identify genes behind competing pathways
and product-consuming reactions.

Available models: 'e_coli_core' (E. coli), 'iMM904' (S. cerevisiae).
Select based on the host organism.

═══ Phase 3 · Enzyme Sourcing (UniProt / InterPro) ═══
For each missing reaction, use uniprot_search with the EC number to find
candidate enzymes across species.  Use uniprot_entry to get kinetics,
cofactors, and structure info for top candidates.  Use interpro_entry to
check protein family and domain architecture.
CRITICAL: Transporters are NOT indexed in KEGG.  You MUST use
uniprot_search with substrate keywords to find them.

═══ Phase 4 · Simulation & Optimization (FBA) ═══
This is a continuous iterative workflow:
 4.1  Add heterologous reactions → verify feasibility
 4.2  maximize_product → theoretical ceiling
 4.3  production_envelope → growth-production tradeoff
 4.4  Metabolite balance analysis → identify bottlenecks
 4.5  simulate_knockout → eliminate competing pathways
 4.6  simulate_overexpression → enhance precursor supply
 4.7  change_media_and_simulate → O₂ / carbon source sweep
 4.8  Transporter engineering → product export

Iterate: each intervention may reveal new bottlenecks.

═══ Output ═══
After completing all phases, output:
1. The complete engineering strategy (gene targets, enzyme sources, conditions)
2. The reasoning chain (which tool calls led to which conclusions)
3. Quantitative FBA predictions (max flux, yield, production envelope)

Summarize findings after each phase before proceeding to the next.
"""


# ══════════════════════════════════════════════════════════════════════════
# Agent runner
# ══════════════════════════════════════════════════════════════════════════

async def run_agent(user_input: str) -> str:
    """Send a prompt to the agent and return the final text response."""
    options = ClaudeAgentOptions(
        model="claude-sonnet-4-6",
        system_prompt=SYSTEM_PROMPT,
        mcp_servers={"synbio": mcp_server},
        max_turns=30,
        env=API_ENV,
    )

    result_text = ""
    async with ClaudeSDKClient(options=options) as client:
        await client.query(user_input)
        async for message in client.receive_response():
            if isinstance(message, ResultMessage):
                result_text = message.result
            elif isinstance(message, AssistantMessage):
                for block in message.content:
                    if isinstance(block, TextBlock):
                        print(f"[Agent] {block.text}")

    return result_text


async def main():
    print("Metabolic Engineering Agent — 4-Phase Framework")
    print("Available models: e_coli_core, iMM904")
    print("Type 'quit' to exit.")
    print("-" * 50)
    while True:
        try:
            user_msg = input("\nYou: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nBye!")
            break
        if not user_msg or user_msg.lower() in ("quit", "exit"):
            print("Bye!")
            break
        answer = await run_agent(user_msg)
        print(f"\nAgent:\n{answer}")


if __name__ == "__main__":
    anyio.run(main)
