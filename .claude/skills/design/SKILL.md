---
name: design
description: "Design a metabolic engineering strategy for producing a target compound in a microbial host. Executes a 4-phase pipeline: (1) KEGG pathway discovery, (2) host GEM assessment, (3) UniProt enzyme sourcing, (4) FBA simulation & optimization. Invoke with /design <product_name> [host_model]"
argument-hint: "<product_name> [host_model: iMM904 | e_coli_core]"
---

# Metabolic Engineering Design — 4-Phase Pipeline

Design a complete engineering strategy for producing **$ARGUMENTS** in a microbial host.

Parse arguments:
- First argument = target product name (required)
- Second argument = host GEM model name (optional, default `iMM904` for *S. cerevisiae*)
- Available models: `iMM904` (yeast), `e_coli_core` (E. coli)

Execute the following phases **sequentially** (1 → 2 → 3 → 3.5 → 4). After each phase, print a summary before moving to the next.

## Available MCP tools (7 total)

| Tool | Actions |
|------|---------|
| `mcp__synbio__kegg` | search, compound, pathway, reaction, orthology |
| `mcp__synbio__gem` | search, gpr, reactions |
| `mcp__synbio__uniprot` | search, entry, interpro |
| `mcp__synbio__protein` | annotate, interactions, kinetics, structure |
| `mcp__synbio__fba` | **Stateful**: reset, add_pathway, knockout, overexpress, media; **Analysis**: maximize, envelope |
| `mcp__synbio__pubmed_search` | *(no action param)* |
| `mcp__synbio__dna_optimize` | *(no action param)* |

---

## Phase 1 · Target Identification (KEGG)

**Goal**: Find the full biosynthetic pathway from product name to every reaction, enzyme (EC), and precursor.

**Steps**:

1. `mcp__synbio__kegg(action="search", database="compound", keyword=<product_name>)` → get compound ID
2. If found: `mcp__synbio__kegg(action="compound", compound_id=<ID>)` → extract PATHWAY, REACTION, ENZYME lists
3. If NOT found OR if the compound has no REACTION/PATHWAY entries:
   - Fall back to `mcp__synbio__kegg(action="search", database="pathway", keyword=...)` and `mcp__synbio__kegg(action="search", database="module", keyword=...)` to find the pathway indirectly
   - Search for the **precursor compound** instead (e.g., for 3'-hydroxygenistein, search genistein first, then find the hydroxylation reaction)
   - If still nothing, use `mcp__synbio__pubmed_search(query=...)` to find biosynthesis literature for pathway hints
4. For each relevant pathway/module: `mcp__synbio__kegg(action="pathway", pathway_id=<ID>)` → get reaction steps and orthology
5. For each key reaction: `mcp__synbio__kegg(action="reaction", reaction_id=<ID>)` → get EC number, KO ID, equation
6. If needed, trace **upstream modules** (e.g. phenylpropanoid pathway upstream of flavonoid pathway)

**Collect**: A complete reaction chain from host precursors to the target product, with EC numbers and KO IDs for every step.

**Print**: Summary table of all reaction steps (substrate → product, EC, KO).

---

## Phase 2 · Host Assessment (GEM)

**Goal**: Build a complete gap list — missing reactions, competing pathways, product degradation, cofactor shortages, compartment issues, and export gaps — that directly drives Phase 3 (enzyme sourcing) and Phase 4 (FBA).

**Steps**:

### Step 1 — Metabolite existence check
For the target product and each Phase 1 intermediate:
`mcp__synbio__gem(action="search", model_name=<model>, keyword=<name>)` → exists? get metabolite ID.

Classify each Phase 1 reaction as **existing** or **missing** (needs heterologous introduction).

### Step 2 — Competing pathways & product degradation (tool-driven, not guessing)
For **every existing Phase 1 intermediate and the target product** (if it exists in the model):
`mcp__synbio__gem(action="reactions", model_name=<model>, metabolite_id=<id>)`
→ Returns ALL consuming reactions (= competing pathways) with stoichiometric coefficients, baseline flux, and GPR.
→ High-flux consumers are the priority knockout/attenuation targets for Phase 4.

Consuming reactions on precursors = **competing pathways**. Consuming reactions on the product = **product degradation** (e.g. CYB2 oxidizes L-lactate back to pyruvate).

### Step 3 — Cofactor supply assessment
From Phase 1, count per-molecule cofactor demands (e.g. 3 NADPH, 1 ATP, 2 NADH per product).
For each heavily consumed cofactor:
`mcp__synbio__gem(action="reactions", model_name=<model>, metabolite_id="nadph_c")`
→ Compare total producing flux vs pathway demand.
→ If supply < demand, flag for Phase 3 (search cofactor-regenerating enzymes) and Phase 4 (overexpression).

### Step 4 — Compartmentalization check
From Step 1 results, note compartment suffixes (_c, _m, _x, _e) for each existing metabolite.
The heterologous pathway typically runs in the cytoplasm (_c). If a key precursor only exists in another compartment (e.g. accoa_m but not accoa_c):
→ Check `reactions` output for transport reactions between compartments.
→ If no transport exists, flag for Phase 3 (search transporter) or Phase 4 (add transport reaction).

### Step 5 — Product export check
Check if an exchange reaction exists for the product:
→ The `reactions` output already reports exchange reactions at the bottom.
→ If no exchange reaction and no native transporter, flag for Phase 3 (search heterologous exporter).

### Step 6 — GPR lookup for key targets
For competing pathway reactions identified in Step 2 and cofactor genes from Step 3:
`mcp__synbio__gem(action="gpr", model_name=<model>, gene_id=<gene>)` → confirm gene-reaction associations before recommending knockouts/overexpression.

**Print**: A structured gap list table:

| Gap Type | Detail | Metabolite/Reaction | Flux | → Phase 3 | → Phase 4 |
|----------|--------|---------------------|------|-----------|-----------|
| Missing reaction | ... | ... | — | Search EC X.X.X.X | Add to extra_reactions |
| Competing pathway | ... | ... | 8.32 | — | Knockout gene Y |
| Product degradation | ... | ... | 0.05 | — | Knockout gene Z |
| Cofactor shortage | Need N NADPH, supply = M | nadph_c | M | Search NADPH enzyme | Overexpress ZWF1 |
| Compartment issue | accoa only in _m | accoa_m | — | Search ACS | Add transport rxn |
| No product export | No EX_ reaction | product_c | — | Search transporter | Add EX reaction |

---

## Phase 3 · Enzyme Sourcing (UniProt / InterPro)

**Goal**: For each missing reaction, find the best candidate heterologous enzyme.

**Steps**:

1. For each missing EC number: `mcp__synbio__uniprot(action="search", query="ec:<EC> AND reviewed:true", max_results=10)` → candidate list
2. For the top 1-2 candidates: `mcp__synbio__uniprot(action="entry", accession=<ID>)` → kinetics (Km, Vmax), cofactors, PDB structures
3. **CRITICAL**: If a transporter is needed for product export, you MUST search UniProt — KEGG does not index transporters:
   `mcp__synbio__uniprot(action="search", query="<substrate> transporter AND reviewed:true")` or search by protein family name
4. Optionally: `mcp__synbio__uniprot(action="interpro", entry_id=<IPR/PF_ID>)` for protein family info
5. Optionally: `mcp__synbio__kegg(action="orthology", ko_id=<KO>)` for cross-species gene distribution

**Prioritize**: Enzymes from the same domain as the host (eukaryotic for yeast), with crystal structures, and favorable kinetics.

**Print**: Table of recommended enzymes (UniProt ID, organism, gene, length, EC, notes).

---

## Phase 3.5 · Protein-Level Analysis

**Goal**: For each heterologous enzyme selected in Phase 3, query protein-level data to identify expression risks, kinetic bottlenecks, known beneficial mutations, and interaction requirements. This data directly informs Phase 4 engineering strategy.

**Steps**:

### Step 1 — Annotation (all heterologous enzymes)
For each candidate enzyme from Phase 3:
`mcp__synbio__protein(action="annotate", accession=<UniProt_ID>)`
→ Check:
  - **Transmembrane domains?** → may need N-terminal truncation or ER membrane expansion
  - **Signal/transit peptide?** → may need removal for cytoplasmic expression in yeast
  - **Cofactor requirements** → cross-reference with Phase 2 cofactor gap list
  - **Mutagenesis data** → known mutations that improve activity or remove regulation

### Step 2 — Kinetics (all heterologous ECs)
For each EC number in the pathway:
`mcp__synbio__protein(action="kinetics", ec_number=<EC>, organism=<source_organism>)`
→ Check:
  - **Km**: lower = higher affinity = better at low intracellular substrate
  - **kcat**: if very low (<1 /s), enzyme may need multi-copy expression or engineering
  - **Ki**: inhibition by product or pathway metabolite → may need de-inhibition mutation
  - **Known mutations**: feedback-resistant variants (e.g. ARO4-K229L), activity-improved variants

### Step 3 — Interactions (for enzymes that need partners)
For enzymes flagged in Step 1 as membrane-bound or requiring cofactors/partners:
`mcp__synbio__protein(action="interactions", protein_id=<gene_name>, organism_taxid=<taxid>)`
→ Check:
  - **Electron transfer partners** (CYP450 → CPR/FDX): which reductase partner is needed?
  - **Protein complex members**: disrupting a complex via KO may have side effects
  - **Co-expression needs**: chaperones, assembly factors

### Step 4 — Structure (optional, for bottleneck enzymes only)
For enzymes where Steps 1-3 identified a kinetic or expression bottleneck:
`mcp__synbio__protein(action="structure", accession=<UniProt_ID>)`
→ Check:
  - **PDB structures available?** → enables rational design
  - **Active site residues** → targets for rational mutagenesis
  - **Cofactor binding sites** → assess electron/substrate access geometry
  - **AlphaFold model confidence** → pLDDT > 70 for reliable analysis

**Decision criteria** — use Phase 3.5 data to:
- Replace wild-type enzymes with known improved variants (e.g. ACC1-S659A/S1157A)
- Add co-expression requirements to Phase 4 (CPR for CYP450, chaperones for large proteins)
- Flag enzymes needing multi-copy integration (low kcat)
- Recommend ER membrane expansion if >2 membrane-bound CYP450 enzymes compete for ER space
- Recommend VHb / heme supply OE if ≥1 heme-dependent CYP450 is in the pathway

**Print**: Updated enzyme table with protein-level data:

| UniProt | Gene | EC | MW | TM? | Km | kcat | Known mutations | Expression notes |
|---------|------|----|----|-----|-----|------|-----------------|------------------|
| Q9SWR5 | IFS2 | 1.14.14.87 | 58 kDa | Yes (ER) | 8 µM | 0.5/s | multi-copy +3.2x | Needs CPR; consider ER expansion |
| P28012 | CHI1 | 5.5.1.6 | 25 kDa | No | 15 µM | — | none | Low risk |

---

## Phase 4 · Simulation & Optimization (FBA)

**Goal**: Validate the pathway in silico, propose a combined engineering strategy, and iteratively optimize it.

The FBA tool is **stateful**: call `reset` to load a fresh model, then apply modifications (add_pathway, knockout, overexpress, media) which persist across calls. Analysis actions (maximize, envelope) read the model without modifying it.

**Workflow**: `reset → add_pathway → [knockout/overexpress/media] → maximize/envelope`. Call `reset` again to start a new scenario.

### 4.1 Baseline

Establish the theoretical maximum with the heterologous pathway alone (no host modifications):

```
mcp__synbio__fba(action="reset", model_name=<model>)
mcp__synbio__fba(action="add_pathway", reactions=[
    {"id": "RXN1", "reaction_string": "sub_c + cofactor_c --> prod_c + ..."},
    {"id": "EX_product", "reaction_string": "product_c -->"},
    ...
])
mcp__synbio__fba(action="maximize", target_reaction="EX_product", min_biomass_fraction=0.0)
mcp__synbio__fba(action="maximize", target_reaction="EX_product", min_biomass_fraction=0.1)
```

The reaction_string must use model metabolite IDs (e.g. `pyr_c + nadh_c + h_c --> lac__L_c + nad_c`).
If the target product already exists but has no exchange reaction, include one in the pathway.

Record the **baseline theoretical max** — this is the ceiling. All subsequent optimization is measured against it.

### 4.2 Iterative Strategy Optimization

This is a **propose → test → interpret → adjust** loop driven by your reasoning.

**Round 1 — Propose initial strategy from Phase 2 gap list:**

Based on Phase 2's competing pathways (ranked by flux), cofactor bottlenecks, and supply gaps, propose a combined strategy and test it in one scenario:

```
mcp__synbio__fba(action="reset", model_name=<model>)
mcp__synbio__fba(action="add_pathway", reactions=[...])
mcp__synbio__fba(action="knockout", knockouts=["gene_A", "gene_B"])
mcp__synbio__fba(action="overexpress", gene_id="gene_C", min_flux=1.0)
mcp__synbio__fba(action="overexpress", gene_id="gene_D", min_flux=0.5)
mcp__synbio__fba(action="maximize", target_reaction="EX_product", min_biomass_fraction=0.1)
```

**Interpret the result:**

- If **infeasible**: one or more interventions are lethal. Remove knockouts one at a time to identify the lethal one. Replace KO with attenuation (reduce flux instead of eliminating).
- If **product flux ≈ baseline** (little improvement): the bottleneck is elsewhere. Examine the top-10 flux reactions in the maximize output to identify where carbon/cofactors are being diverted. Propose new interventions targeting those reactions.
- If **product flux improved significantly**: good. Check if further improvement is possible.

**Round 2+ — Adjust and re-test:**

Based on Round 1 results, adjust the strategy:
- Remove lethal KOs, add alternative targets
- Adjust overexpression levels (too high → infeasible, too low → no effect)
- Try different media conditions
- Target new competing reactions revealed by the top-10 flux list

```
mcp__synbio__fba(action="reset", model_name=<model>)
mcp__synbio__fba(action="add_pathway", reactions=[...])
← apply adjusted interventions →
mcp__synbio__fba(action="maximize", target_reaction="EX_product", min_biomass_fraction=0.1)
```

**Stopping criteria (stop when ANY is met):**

1. Product flux ≥ 95% of baseline theoretical max → already near optimum
2. This round improved product flux by < 5% over the previous round → diminishing returns
3. Model is infeasible and all reasonable adjustments have been tried → report best feasible result
4. Completed 3 rounds → safety cap, stop and report best result

### 4.3 Final Envelope

Once the best strategy is determined, generate the production envelope:

```
← model should still be in the best-strategy state from 4.2 →
mcp__synbio__fba(action="envelope", target_reaction="EX_product", steps=10)
```

If the model was reset in the last iteration, rebuild the best strategy first, then run envelope.

**Print**: Scenario comparison table across all rounds, production envelope for the final strategy.

---

## Final Output

After all phases, output a structured report with 6 sections:

### 1. Biosynthetic Pathway
Complete reaction chain from host precursors to the target product.
Combine Phase 1 (reactions) + Phase 2 (native/missing) + Phase 3 (enzyme source):

| Step | Substrate → Product | EC | Enzyme | Source | UniProt | Native? |
|------|--------------------|----|--------|--------|---------|---------|
| 1 | ... → ... | ... | ... | *S. cerevisiae* | — | Native |
| 2 | ... → ... | ... | ... | *E. coli* | Q57160 | **Heterologous** |
| ... | | | | | | |

Cofactor demands per product molecule: N NADPH, M ATP, K O2, ...

### 2. Host Modifications
Three sub-tables, each with evidence tracing:

**a) Heterologous genes to express** (from Phase 3 + Phase 3.5):

| Gene | Organism | UniProt | MW | Km | kcat | TM? | Mutations | Notes |
|------|----------|---------|-----|-----|------|-----|-----------|-------|
| HpaB | *P. aeruginosa* | Q9HWT7 | 58 kDa | — | — | No | — | PDB: 6QYH |
| IFS2 | *G. max* | Q9SWR5 | 58 kDa | 8µM | 0.5/s | Yes (ER) | multi-copy +3.2x | Needs CPR; ER expansion recommended |

**b) Genes to knock out** (from Phase 2 competing pathways, validated in Phase 4):

| Gene | Yeast ID | Reaction blocked | Competing flux (WT) | Rationale |
|------|----------|-----------------|---------------------|-----------|
| PDC1 | YLR044C | PYRDC (pyr→acald) | 15.95 | Top pyruvate competitor |

**c) Genes to overexpress** (from Phase 2 bottlenecks, validated in Phase 4):

| Gene | Yeast ID | Reaction boosted | Baseline flux | Forced to | Rationale |
|------|----------|-----------------|---------------|-----------|-----------|
| ZWF1 | YNL241C | G6PDH2r | 0.80 | 2.00 | NADPH supply bottleneck |

### 3. Quantitative Predictions
From Phase 4 FBA iterative optimization:

**Optimization rounds:**

| Round | Strategy | Biomass | Product flux | vs Baseline | Notes |
|-------|----------|---------|-------------|-------------|-------|
| Baseline | Pathway only | ... | ... | — | Theoretical ceiling |
| Round 1 | + KO A,B + OE C,D | ... | ... | +X% | KO B was lethal → removed |
| Round 2 | + KO A + OE C,D adjusted | ... | ... | +Y% | Stopped: <5% improvement |

**Production envelope** (growth vs production tradeoff for final strategy).

### 3a. Protein Engineering Notes (from Phase 3.5)
- Enzymes requiring signal peptide removal/modification for expression in yeast
- Known beneficial mutations to apply (feedback-resistant variants, activity-improved variants)
- Interaction partner co-expression requirements (CPR for CYP450, chaperones)
- ER membrane expansion / VHb / heme supply recommendations (if CYP450 present)
- Enzymes recommended for multi-copy integration (low kcat)

### 4. Recommended Conditions
- Carbon source
- Oxygen level (aerobic / microaerobic / anaerobic)
- Special supplements if needed (e.g., lactose feeding for HMO production)

### 5. Known Limitations
What this analysis cannot predict:
- Actual protein expression levels and in vivo solubility (Phase 3.5 provides annotations but not quantitative predictions)
- Product toxicity and feedback inhibition at high titers
- Adaptive laboratory evolution (ALE) outcomes
- Actual titers (FBA gives flux, not g/L)
- Mutation effects beyond what is in the literature (no de novo prediction)
