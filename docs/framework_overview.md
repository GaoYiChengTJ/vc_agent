# Metabolic Engineering Agentic System — Framework Overview

## System Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                          User Input                                 │
│                  Product Name + Host Model                          │
│              e.g. "Hydroxytyrosol", "iMM904"                       │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│                     LLM Orchestrator                                │
│                      (Claude)                                       │
│                                                                     │
│  Responsibilities:                                                  │
│  · Interpret user intent                                            │
│  · Sequence tool calls across 5 phases                              │
│  · Reason about intermediate results                                │
│  · Propose and iteratively refine engineering strategies             │
│  · Synthesize final report                                          │
└────┬──────────┬──────────┬──────────┬──────────┬──────────┬─────────┘
     │          │          │          │          │          │
     ▼          ▼          ▼          ▼          ▼          ▼
  ┌──────┐  ┌──────┐  ┌─────────┐ ┌──────┐  ┌────────┐ ┌────────┐
  │ kegg │  │ gem  │  │ uniprot │ │ fba  │  │pubmed  │ │dna_opt │
  │      │  │      │  │         │ │      │  │_search │ │_imize  │
  └──┬───┘  └──┬───┘  └────┬────┘ └──┬───┘  └───┬────┘ └───┬────┘
     │         │           │         │           │          │
     ▼         ▼           ▼         ▼           ▼          ▼
  ┌──────┐  ┌──────┐  ┌─────────┐ ┌──────┐  ┌────────┐ ┌────────┐
  │ KEGG │  │ GEM  │  │ UniProt │ │COBRApy│ │PubMed  │ │Codon   │
  │ REST │  │Models│  │ REST    │ │Solver │ │  API   │ │Tables  │
  │ API  │  │(SBML)│  │ API     │ │(GLPK) │ │        │ │        │
  └──────┘  └──────┘  └─────────┘ └──────┘  └────────┘ └────────┘
  External   Local      External   Local     External    Local
```

## Tool Inventory

| Tool | Actions | Type | Data Source | Purpose |
|------|---------|------|-------------|---------|
| **kegg** | search, compound, pathway, reaction, orthology | Stateless | KEGG REST API | Pathway discovery |
| **gem** | search, gpr, reactions | Stateless | Local SBML models | Host model queries |
| **uniprot** | search, entry, interpro | Stateless | UniProt/InterPro REST API | Enzyme sourcing |
| **protein** | annotate, interactions, kinetics, structure | Stateless | UniProt/STRING/BRENDA/PDB/AlphaFold | Protein-level analysis |
| **fba** | reset, add_pathway, knockout, overexpress, media, maximize, envelope | **Stateful** | Local COBRApy solver | Flux simulation |
| **pubmed_search** | — | Stateless | PubMed API | Literature evidence |
| **dna_optimize** | — | Stateless | Local codon tables | Codon optimization |

## 5-Phase Pipeline

```
Phase 1             Phase 2            Phase 3           Phase 4            Phase 5
Pathway Discovery   Host Assessment    Enzyme Sourcing   Protein Analysis   FBA Optimization
─────────────────   ───────────────    ───────────────   ────────────────   ────────────────

┌──────────────┐   ┌──────────────┐   ┌──────────────┐  ┌──────────────┐  ┌──────────────┐
│              │   │              │   │              │  │              │  │              │
│ KEGG search  │   │ Metabolite   │   │ UniProt      │  │ Annotate     │  │ FBA reset    │
│      │       │   │ existence    │   │ search by EC │  │ (TM, signal  │  │      │       │
│      ▼       │   │ check        │   │      │       │  │  peptide,    │  │      ▼       │
│ Compound     │   │      │       │   │      ▼       │  │  cofactors,  │  │ Add pathway  │
│ details      │   │      ▼       │   │ Entry detail │  │  mutagenesis)│  │      │       │
│      │       │   │ Metabolite   │   │ (Km, PDB)    │  │      │       │  │      ▼       │
│      ▼       │   │ reactions    │   │      │       │  │      ▼       │  │ ┌──────────┐ │
│ Pathway /    │   │ (competing,  │   │      ▼       │  │ Kinetics     │  │ │ Propose  │ │
│ Module       │   │  supply,     │   │ Transporter  │  │ (Km, kcat,   │  │ │ strategy │ │
│ details      │   │  cofactor,   │   │ search       │  │  Ki, known   │  │ │ from gap │ │
│      │       │   │  transport,  │   │ (if needed)  │  │  mutations)  │  │ │ list     │ │
│      ▼       │   │  exchange)   │   │              │  │      │       │  │ └─────┬────┘ │
│ Reaction     │   │      │       │   │              │  │      ▼       │  │       │      │
│ details      │   │      ▼       │   │              │  │ Interactions │  │       ▼      │
│ (EC, KO,     │   │ GPR lookup   │   │              │  │ (CYP→CPR,   │  │ ┌──────────┐ │
│  equation)   │   │ for targets  │   │              │  │  complexes)  │  │ │  Test    │ │
│      │       │   │              │   │              │  │      │       │  │ │(maximize)│ │
│      ▼       │   │              │   │              │  │      ▼       │  │ └─────┬────┘ │
│ [PubMed      │   │              │   │              │  │ Structure    │  │       │      │
│  fallback]   │   │              │   │              │  │ (optional,   │  │       ▼      │
│              │   │              │   │              │  │  bottleneck  │  │ ┌──────────┐ │
│              │   │              │   │              │  │  enzymes)    │  │ │Interpret │ │
└──────┬───────┘   └──────┬───────┘   └──────┬───────┘  └──────┬───────┘  │ │results  │─┼─┐
       │                  │                  │                 │          │ └─────┬────┘ │ │
       ▼                  ▼                  ▼                 ▼          │       │      │ │
  Reaction          Gap List           Enzyme List       Protein-Level   │       ▼      │ │
  Chain Table                                           Risk Table       │ ┌──────────┐ │ │
                                                                        │ │ Adjust   │ │ │
                                                                        │ │ strategy │◄┼─┘
                                                                        │ └─────┬────┘ │
                                                                        │       │      │
                                                                        │  Stop?──No──►│
                                                                        │   │Yes       │
                                                                        │   ▼          │
                                                                        │  Envelope    │
                                                                        └──────┬───────┘
                                                                               │
                                                                               ▼
                                                                         Best Strategy
                                                                         + Envelope
```

## Data Flow Between Phases

```
Phase 1 ──────────────────────────────────► Phase 2
  OUTPUT:                                    INPUT:
  · Reaction chain table                     · Intermediate metabolite names
    (substrate→product, EC, KO)              · Cofactor demands per molecule
  · Cofactor tally (N NADPH, M ATP)          · Target product name

Phase 2 ──────────────────────────────────► Phase 3
  OUTPUT:                                    INPUT (from gap list):
  · Gap list table:                          · Missing EC numbers → search enzymes
    - Missing reactions (EC)                 · "No export" flag → search transporter
    - Competing pathways (gene, flux)        · Cofactor shortage → search regenerating enzyme
    - Cofactor shortages
    - Compartment issues
    - Export gaps

Phase 3 ──────────────────────────────────► Phase 4
  OUTPUT:                                    INPUT:
  · Enzyme list (UniProt IDs,                · UniProt accessions → annotate
    organism, EC)                            · EC numbers → kinetics
                                             · Gene names → interactions

Phase 2 + 4 ─────────────────────────────► Phase 5
  Phase 2 OUTPUT:                            INPUT:
  · Competing pathway genes + flux           · KO targets (ranked by flux)
  · Supply bottleneck genes                  · OE targets
  · Cofactor gaps                            · Media hints (aerobic requirement?)
  Phase 4 OUTPUT:                            · Co-expression needs (CPR for CYP450)
  · Protein risk table (TM, Km, kcat,        · Known beneficial mutations
    partner needs, mutations)                · Multi-copy flags (low kcat enzymes)
                                             · Reaction equations in COBRApy format

Phase 1 + 2 + 3 + 4 + 5 ────────────────► Final Output
```

## Phase 2 Detail: Host Assessment

```
                    Phase 1 output
                   (metabolite names)
                         │
                         ▼
              ┌─────────────────────┐
              │ Step 1: Search      │  gem(action="search")
              │ Map names → IDs     │  × N intermediates
              └──────────┬──────────┘
                         │
                    existing IDs         missing names
                    [pyr_c, ...]         [naringenin, ...]
                         │                     │
                         ▼                     ▼
              ┌─────────────────────┐    (record as
              │ Step 2: Reactions   │    "missing reaction"
              │ For EVERY existing  │     in gap list)
              │ intermediate +      │
              │ target product      │  gem(action="reactions")
              │                     │  × M existing metabolites
              │ Returns:            │
              │ · Consuming rxns    │──► Competing pathways
              │   (with flux, GPR)  │
              │ · Producing rxns    │──► Supply sources
              │   (with flux)       │
              │ · Other compartments│──► Compartment info
              │ · Exchange rxns     │──► Export info
              └──────────┬──────────┘
                         │
          ┌──────────────┼──────────────┐
          ▼              ▼              ▼
   ┌────────────┐ ┌────────────┐ ┌────────────┐
   │ Step 3     │ │ Step 4     │ │ Step 5     │
   │ Cofactor   │ │ Compartment│ │ Product    │
   │ supply     │ │ check      │ │ export     │
   │ assessment │ │ (_c/_m/_e) │ │ check      │
   │            │ │            │ │            │
   │ reactions  │ │ (read from │ │ (read from │
   │ ("nadph_c")│ │  Step 2    │ │  Step 2    │
   │ supply vs  │ │  output)   │ │  output)   │
   │ demand     │ │            │ │            │
   └─────┬──────┘ └─────┬──────┘ └─────┬──────┘
         │              │              │
         └──────────────┼──────────────┘
                        ▼
              ┌─────────────────────┐
              │ Step 6: GPR lookup  │  gem(action="gpr")
              │ Confirm gene-rxn    │  × K target genes
              │ associations        │
              └──────────┬──────────┘
                         │
                         ▼
                ┌─────────────────┐
                │   GAP LIST      │
                │                 │
                │ · Missing rxns  │──► Phase 3 (search enzymes)
                │ · Competing     │──► Phase 5 (KO targets)
                │ · Degradation   │──► Phase 5 (KO targets)
                │ · Cofactor gaps │──► Phase 3 + Phase 5
                │ · Compartment   │──► Phase 3 (transporter) + Phase 5
                │ · Export gaps   │──► Phase 3 (exporter) + Phase 5
                └─────────────────┘
```

## Phase 5 Detail: Iterative FBA Optimization

```
          Phase 2 gap list          Phase 3+4 enzyme+protein data
          (KO/OE targets)           (reaction equations, risk table)
                │                          │
                └────────────┬─────────────┘
                             ▼
                  ┌─────────────────────┐
                  │ 5.1 BASELINE        │
                  │                     │
                  │ reset → add_pathway │
                  │ → maximize(0%)      │  Theoretical ceiling
                  │ → maximize(10%)     │  (pathway only, no host mods)
                  └──────────┬──────────┘
                             │
                      baseline_max = X
                             │
                             ▼
                  ┌─────────────────────┐
            ┌────►│ 5.2 PROPOSE         │
            │     │                     │
            │     │ LLM selects from    │
            │     │ gap list:           │
            │     │ · KO targets        │
            │     │ · OE targets        │
            │     │ · Media conditions  │
            │     └──────────┬──────────┘
            │                │
            │                ▼
            │     ┌─────────────────────┐
            │     │ TEST                │
            │     │                     │
            │     │ reset → add_pathway │
            │     │ → knockout(...)     │
            │     │ → overexpress(...)  │
            │     │ → maximize(10%)     │
            │     └──────────┬──────────┘
            │                │
            │          result = Y
            │                │
            │                ▼
            │     ┌─────────────────────┐
            │     │ INTERPRET           │
            │     │                     │
            │     │ Infeasible?         │──► Remove lethal KO
            │     │ Y ≈ baseline?       │──► Read top-10 flux,
            │     │                     │    find new targets
            │     │ Y >> baseline?      │──► Check stopping criteria
            │     └──────────┬──────────┘
            │                │
            │                ▼
            │     ┌─────────────────────┐
            │     │ STOP?               │
            │     │                     │
            │     │ · flux ≥ 95% max?   │──Yes──►┐
            │     │ · improvement < 5%? │──Yes──►│
            │     │ · all options tried?│──Yes──►│
            │     │ · 3 rounds done?    │──Yes──►│
            │     │                     │        │
            │     └──────────┬──────────┘        │
            │                │No                 │
            │                ▼                   │
            │     ┌─────────────────────┐        │
            │     │ ADJUST              │        │
            │     │                     │        │
            │     │ · Remove lethal KOs │        │
            │     │ · Adjust OE levels  │        │
            │     │ · Add new targets   │        │
            │     │   from top-10 flux  │        │
            └─────┤ · Try alt. media    │        │
                  └─────────────────────┘        │
                                                 │
                                                 ▼
                                      ┌─────────────────────┐
                                      │ 5.3 FINAL ENVELOPE  │
                                      │                     │
                                      │ envelope(steps=10)  │
                                      │ on best strategy    │
                                      └──────────┬──────────┘
                                                 │
                                                 ▼
                                          Best Strategy
                                        + Scenario Table
                                        + Production Envelope
```

## Final Output Structure

```
┌─────────────────────────────────────────────────────────────────┐
│                        FINAL REPORT                             │
│                                                                 │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │ 1. BIOSYNTHETIC PATHWAY                                   │  │
│  │    Complete reaction chain table                          │  │
│  │    Step | Substrate→Product | EC | Enzyme | Source |      │  │
│  │    Native? | UniProt                                      │  │
│  │    + Cofactor demands per molecule                        │  │
│  │                                      ◄── Phase 1+2+3     │  │
│  └───────────────────────────────────────────────────────────┘  │
│                                                                 │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │ 2. HOST MODIFICATIONS                                     │  │
│  │    a) Heterologous genes to express                       │  │
│  │       Gene | Organism | UniProt | MW | Km | kcat |        │  │
│  │       TM? | Mutations | Notes                             │  │
│  │                                      ◄── Phase 3+4       │  │
│  │    b) Genes to knock out                                  │  │
│  │       Gene | Yeast ID | Reaction | WT flux | Rationale    │  │
│  │                                      ◄── Phase 2+5       │  │
│  │    c) Genes to overexpress                                │  │
│  │       Gene | Yeast ID | Reaction | Baseline | Forced      │  │
│  │                                      ◄── Phase 2+5       │  │
│  └───────────────────────────────────────────────────────────┘  │
│                                                                 │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │ 3. QUANTITATIVE PREDICTIONS                               │  │
│  │    Optimization rounds table                              │  │
│  │    Round | Strategy | Biomass | Product | vs Baseline     │  │
│  │    + Production envelope                                  │  │
│  │                                      ◄── Phase 5         │  │
│  └───────────────────────────────────────────────────────────┘  │
│                                                                 │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │ 3a. PROTEIN ENGINEERING NOTES                             │  │
│  │    · Signal peptide / TM domain handling                  │  │
│  │    · Known beneficial mutations                           │  │
│  │    · Co-expression partners (CPR for CYP450)              │  │
│  │    · Multi-copy integration flags (low kcat)              │  │
│  │    · ER expansion / heme supply recommendations           │  │
│  │                                      ◄── Phase 4         │  │
│  └───────────────────────────────────────────────────────────┘  │
│                                                                 │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │ 4. RECOMMENDED CONDITIONS                                 │  │
│  │    Carbon source, oxygen, supplements                     │  │
│  │                                      ◄── Phase 5         │  │
│  └───────────────────────────────────────────────────────────┘  │
│                                                                 │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │ 5. KNOWN LIMITATIONS                                      │  │
│  │    What the system cannot predict                         │  │
│  └───────────────────────────────────────────────────────────┘  │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## FBA Tool State Machine

```
                    ┌──────────┐
                    │  No      │
          ┌────────►│  Model   │◄────────────────────┐
          │         │  Loaded  │                      │
          │         └────┬─────┘                      │
          │              │                            │
          │         reset(model)                      │
          │              │                            │
          │              ▼                            │
          │         ┌──────────┐                      │
          │         │  Clean   │                      │
          │    ┌───►│  Model   │◄──────┐              │
          │    │    │  Loaded  │       │              │
          │    │    └────┬─────┘       │              │
          │    │         │             │              │
          │    │    add_pathway   reset(model)        │
          │    │         │             │              │
          │    │         ▼             │              │
          │    │    ┌──────────┐       │              │
          │    │    │ Pathway  │       │              │
          │    │    │ Added    │───────┤              │
          │    │    └────┬─────┘       │              │
          │    │         │             │              │
          │    │    knockout /         │              │
          │    │    overexpress /      │              │
          │    │    media              │              │
          │    │         │             │              │
          │    │         ▼             │              │
          │    │    ┌──────────┐       │              │
          │    │    │ Modified │       │              │
          │    │    │ Model    │───────┘              │
          │    │    └────┬─────┘                      │
          │    │         │                            │
          │    │    knockout /                        │
          │    │    overexpress /                     │
          │    │    media (accumulate)                │
          │    │         │                            │
          │    │         ▼                            │
          │    │    ┌──────────────┐                  │
          │    │    │ maximize /   │ ◄── read-only    │
          │    │    │ envelope     │     (with model:)│
          │    │    │              │                  │
          │    │    │ Returns      │                  │
          │    │    │ results      │                  │
          │    │    └──────┬───────┘                  │
          │    │           │                          │
          │    │      LLM decides                    │
          │    │           │                          │
          │    │     ┌─────┴──────┐                   │
          │    │     ▼            ▼                   │
          │    │  Continue    New scenario             │
          │    │  modifying   (reset)                  │
          │    │     │            │                    │
          │    └─────┘            └────────────────────┘
          │
          └── reset(different_model)

  Modification actions: reset, add_pathway, knockout, overexpress, media
     → change _WORKING_MODEL in place, persist across calls

  Analysis actions: maximize, envelope
     → use `with model:` context, do NOT modify _WORKING_MODEL
```

## Key Design Principles

1. **Tool atomicity**: Each tool action does ONE thing. Modifications and analysis are separate.

2. **Stateful FBA**: Model modifications accumulate across calls. `reset` starts a new scenario. Analysis actions are read-only.

3. **Phase 2 is tool-driven**: `gem(action="reactions")` systematically finds ALL competing/producing reactions via stoichiometry, not LLM keyword guessing.

4. **Phase 5 is LLM-driven**: The LLM proposes strategies based on Phase 2 evidence, tests via FBA, interprets results, and iterates. Stopping criteria are outcome-based with a safety cap.

5. **Data flows forward**: Phase 1 → reaction chain → Phase 2 → gap list → Phase 3 → enzyme list → Phase 4 → protein risk table → Phase 5 → validated strategy → Final Output. Each phase's output is the next phase's input.

6. **Final output is user-actionable**: Contains everything a researcher needs to start lab work — specific genes to order, specific genes to knock out, expected yields, and known risks.
