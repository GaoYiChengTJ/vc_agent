# Diagram Prompt for Metabolic Engineering Agentic System

Draw a system architecture diagram for an LLM-powered metabolic engineering design agent. The system takes a target compound name (e.g., "Hydroxytyrosol") and a host organism model as input, and outputs a complete engineering strategy report. Use a clean, professional style suitable for a research paper.

## Overall structure

The diagram should show three layers from top to bottom:

**Top layer: User Input** — a single box labeled "User Input: Product Name + Host Model".

**Middle layer: LLM Orchestrator** — a large central box labeled "LLM Orchestrator (Claude)". This is the brain of the system. It calls tools, reasons about results, and drives a 4-phase sequential pipeline. All tool calls originate from this box.

**Bottom layer: 6 MCP Tools** — six tool boxes arranged in a row beneath the LLM. Each tool connects to its external/local data source below it:

| Tool | Actions | Data Source | Source Type |
|------|---------|-------------|-------------|
| kegg | search, compound, pathway, reaction, orthology | KEGG REST API | External |
| gem | search, gpr, reactions | GEM Models (SBML) | Local |
| uniprot | search, entry, interpro | UniProt REST API | External |
| fba | reset, add_pathway, knockout, overexpress, media, maximize, envelope | COBRApy Solver (GLPK) | Local |
| pubmed_search | (single action) | PubMed API | External |
| dna_optimize | (single action) | Codon Tables | Local |

The fba tool should be visually distinct (highlighted border or different color) because it is **stateful** — modifications persist across calls. All other tools are stateless.

## 4-Phase Pipeline (inside or alongside the LLM box)

Show 4 phases flowing left to right (or top to bottom), connected by arrows carrying their data outputs:

**Phase 1: Pathway Discovery**
- Tools used: kegg, pubmed_search (fallback)
- Input: product name
- Output: **Reaction Chain Table** (substrate→product, EC number, KO ID for each step) + cofactor tally
- Data arrow: "Reaction Chain + Cofactor Demands" → Phase 2

**Phase 2: Host Assessment**
- Tools used: gem (search → reactions → gpr)
- Input: Reaction Chain from Phase 1
- Internal workflow (6 steps):
  1. Metabolite existence check (gem search) → classify existing vs missing
  2. For every existing intermediate: gem reactions → get ALL consuming/producing reactions with flux + GPR
  3. Cofactor supply assessment (gem reactions on nadph_c etc.)
  4. Compartmentalization check (read compartment info from step 2)
  5. Product export check (read exchange info from step 2)
  6. GPR lookup for key targets
- Output: **Gap List** — a structured table with 6 gap types: missing reactions, competing pathways (with flux), product degradation, cofactor shortages, compartment issues, export gaps
- Data arrow: Gap List splits into two:
  - "Missing ECs, Transport/Export gaps" → Phase 3
  - "KO targets (ranked by flux), OE targets, Media hints" → Phase 4

**Phase 3: Enzyme Sourcing**
- Tools used: uniprot, kegg (orthology)
- Input: Gap List from Phase 2 (missing ECs, transporter needs)
- Output: **Enzyme List** — recommended enzymes with UniProt ID, organism, kinetics, cofactors
- Data arrow: "Enzyme List + Reaction Equations (COBRApy format)" → Phase 4

**Phase 4: FBA Simulation & Iterative Optimization**
- Tools used: fba (stateful)
- Input: Gap List (KO/OE targets) from Phase 2 + Enzyme List (reaction equations) from Phase 3
- Internal workflow — show this as an iterative loop:
  - **4.1 Baseline**: reset → add_pathway → maximize → get theoretical ceiling
  - **4.2 Iterative loop** (show with a circular arrow):
    - PROPOSE: LLM selects KOs + OEs from gap list
    - TEST: reset → add_pathway → knockout → overexpress → maximize
    - INTERPRET: LLM reads result (infeasible? improved? stagnant?)
    - STOP CHECK: flux ≥ 95% max? improvement < 5%? 3 rounds done?
    - If no: ADJUST strategy → loop back to PROPOSE
    - If yes: exit loop
  - **4.3 Final Envelope**: envelope on best strategy
- Output: **Best Strategy + Scenario Comparison Table + Production Envelope**

## Final Output (rightmost or bottom box)

A report box with 5 sections:
1. **Biosynthetic Pathway** — complete reaction chain with enzyme sources (from Phase 1+2+3)
2. **Host Modifications** — three sub-tables: (a) heterologous genes to express, (b) genes to knock out, (c) genes to overexpress (from Phase 2+3+4)
3. **Quantitative Predictions** — optimization rounds table + production envelope (from Phase 4)
4. **Recommended Conditions** — carbon source, oxygen level (from Phase 4)
5. **Known Limitations** — what FBA cannot predict

## Key visual elements to emphasize

1. **Data flow arrows** between phases should be labeled with the data type being passed (Reaction Chain, Gap List, Enzyme List, Best Strategy).

2. **Phase 4 iterative loop** should be visually distinct — use a circular arrow or loop indicator to show the propose→test→interpret→adjust cycle, with stopping criteria noted.

3. **Phase 2's gem(reactions) tool** is the key innovation — it systematically queries the metabolic model for all consuming/producing reactions of each metabolite (by stoichiometry, not keyword guessing), returning flux values and GPR in one call.

4. **FBA stateful nature** — show that modifications (add_pathway, knockout, overexpress, media) accumulate on a working model, while analysis actions (maximize, envelope) are read-only.

5. **LLM reasoning** happens between tool calls in Phase 4: the LLM interprets FBA results and decides whether to continue iterating or stop. This is not a fixed pipeline — the LLM drives the optimization loop.
