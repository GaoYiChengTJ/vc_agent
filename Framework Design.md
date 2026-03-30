# Metabolic Engineering Agentic System — 4-Phase Framework Design

## Overview

This document describes a generalized 4-phase agentic framework for designing metabolic engineering strategies in microbial hosts. Given a **(target product, host organism)** pair, the agent systematically discovers biosynthetic pathways, evaluates host capacity, sources heterologous enzymes, and performs simulation-driven optimization through constraint-based modeling.

The framework has been validated on three structurally diverse product classes:

| Case | Product | Class | Pathway Length | Consistency vs. Literature |
|------|---------|-------|---------------|---------------------------|
| 1 | L-Lactic acid | Organic acid | 1 step | ~85% |
| 2 | Lacto-N-neotetraose (LNnT) | Oligosaccharide (HMO) | 3 steps + precursors | ~85% |
| 3 | Biochanin A | Plant isoflavone (CYP450) | 8 steps | ~60% |

---

## Architecture

Each phase answers one question using a distinct tool set. The output of each phase feeds directly into the next.

```
Input: (target_product_name, host_organism)
                │
        ┌───────▼────────┐
        │  Phase 1        │  "What is this molecule and how is it made?"
        │  Target ID      │
        │  (KEGG)         │
        └───────┬────────┘
                │  → compound_id, pathway, reactions[], enzymes[{EC, KO}]
                │
        ┌───────▼────────┐
        │  Phase 2        │  "What does the host already have? What's missing?"
        │  Host Assessment│
        │  (GEM)          │
        └───────┬────────┘
                │  → existing[], missing[], competing_pathways[], baseline_fluxes
                │
        ┌───────▼────────┐
        │  Phase 3        │  "Where do we get the missing enzymes?"
        │  Enzyme Sourcing│
        │  (UniProt/IPR)  │
        └───────┬────────┘
                │  → candidate_enzymes[{UniProt_ID, organism, Km, PDB}]
                │
        ┌───────▼─────────────────────────────────┐
        │  Phase 4                                  │
        │  Simulation & Optimization (FBA)          │
        │                                           │
        │  4.1 Add pathway → verify feasibility     │
        │  4.2 Max theoretical flux                 │
        │  4.3 Production envelope                  │
        │  4.4 Bottleneck identification             │
        │  4.5 Knockout / overexpression scan        │
        │  4.6 Media & condition optimization        │
        │  4.7 Transporter engineering               │
        └───────┬─────────────────────────────────┘
                │
Output: Complete engineering strategy + reasoning chain

─────────────────────────────────────────────────────
Optional post-processing:
        ┌─────────────────┐
        │  Seq. Engineering│  Codon optimization (DNAchisel)
        │  (on demand)     │  Standard operation, not a design decision
        └─────────────────┘
```

### Why 4 Phases, Not More

- **Phases 1–3** each use a distinct external knowledge source (KEGG, GEM, UniProt) to answer a distinct question. They cannot be merged.
- **Phase 4** merges the former "Pathway Validation" and "Optimization Strategies" into a single continuous FBA workflow. Both use the same tool set (COBRApy/FBA), and "validation" is simply the first iteration of "optimization" — separating them creates an artificial boundary.
- **Sequence engineering** (codon optimization) is excluded from the core loop: it is a deterministic post-processing step that involves no design reasoning, and none of the three ground-truth papers mention it as part of their strategy.

---

## Phase 1: Target Identification

**Question**: "What is this molecule? How is it biosynthesized?"

**Goal**: Starting from a product name, discover the full biosynthetic route — every reaction, enzyme, and precursor.

### Logic Flow

```
search_kegg(compound, product_name)
  │
  ├─ HIT → query_kegg_compound(compound_id)
  │          ├─ Extract: PATHWAY[], REACTION[], ENZYME[] lists
  │          ├─ For each pathway → query_kegg_pathway(pathway_id)
  │          ├─ For each module  → query_kegg_pathway(module_id)
  │          └─ For each reaction → query_kegg_reaction(reaction_id)
  │                                  └─ Extract: EC number, ORTHOLOGY (KO ID)
  │
  └─ MISS → Fallback search:
             ├─ search_kegg(pathway, product_keywords)
             ├─ search_kegg(module, product_keywords)
             └─ search_kegg(reaction, product_keywords)
                 └─ Infer pathway from structural analogs or homologous pathways
```

### Decision Points

- **Direct KEGG hit** (L-Lactic acid C00186, Biochanin A C00814): Follow compound → pathway → reaction → enzyme chain directly.
- **No direct hit** (LNnT): Search by pathway/module keywords. LNnT was found via the neolacto-series glycosphingolipid pathway (M00071) which shares the same sugar chain structure.
- **Long pathways**: Trace upstream modules. Biochanin A required chaining M00137 (phenylalanine → naringenin) with M00941 (naringenin → genistein) plus R02931 (genistein → biochanin A).

### Tools Required

| Tool | Function |
|------|----------|
| `search_kegg(database, keyword)` | Search compound/reaction/enzyme/module/pathway |
| `query_kegg_compound(compound_id)` | Get pathways, reactions, enzymes for a compound |
| `query_kegg_pathway(pathway_or_module_id)` | Get pathway/module detail with reaction steps and orthology |
| `query_kegg_reaction(reaction_id)` | Get reaction equation, EC number, KO ID |

### Output Schema

```json
{
  "compound_id": "C00814",
  "pathways": ["map00943", "M00941", "M00137"],
  "reaction_chain": [
    {"id": "R00697", "substrate": "L-Phe", "product": "Cinnamate", "EC": "4.3.1.24", "KO": "K10775"},
    {"id": "R02253", "substrate": "Cinnamate", "product": "p-Coumarate", "EC": "1.14.14.91", "KO": "K00487"},
    "..."
  ],
  "precursors": ["L-Phenylalanine", "Malonyl-CoA", "SAM", "NADPH", "O2"]
}
```

### Validated Performance

| Case | Search Path | Result |
|------|------------|--------|
| L-Lactic acid | `C00186 → R00703 → EC 1.1.1.27` | 1-step, direct |
| LNnT | `No hit → pathway:"lacto" → M00071 → R05971/R05977` | 2-step, indirect |
| Biochanin A | `C00814 → M00137 + M00941 + R02931` | 8-step, chained modules |

### Known Limitations

- KEGG provides **canonical (natural) pathways**, not engineered shortcuts. Example: the TAL shortcut (Tyr → p-coumaric acid directly, bypassing C4H) is widely used in metabolic engineering but is not the route KEGG returns for flavonoid biosynthesis.
- Possible mitigation: maintain a curated "engineering shortcut" knowledge base alongside KEGG.

---

## Phase 2: Host Assessment

**Question**: "What does the host already have? What's missing?"

**Goal**: Map Phase 1 outputs onto the host's genome-scale metabolic model (GEM). Identify which metabolites and reactions exist, which are missing, what competing pathways drain shared precursors, and what baseline flux the host provides.

### Logic Flow

```
For each metabolite/reaction from Phase 1:
  │
  ├─ search_model(model, product_keyword)     → Target product exists?
  ├─ search_model(model, precursor_keyword)   → Precursor metabolites exist?
  ├─ search_model(model, reaction_keyword)    → Biosynthetic reactions exist?
  │
  ├─ Existing reactions → query_gpr(model, gene_id) → Which genes are associated?
  │
  ├─ Missing reactions → Flag as "requires heterologous introduction"
  │
  └─ Baseline FBA:
       ├─ model.optimize() → biomass, exchange fluxes
       ├─ Identify competing pathways (reactions consuming shared precursors)
       │    e.g., PYRDC competes for pyruvate, FAS competes for malonyl-CoA
       └─ Identify product-consuming reactions
            e.g., CYB2 oxidizes L-lactate back to pyruvate
```

### Decision Points

- **Target metabolite exists in model** (e.g., `lac__L_c` for L-lactate): Check whether a *synthesis* reaction also exists. Existence of the metabolite alone is insufficient.
- **Target metabolite absent** (e.g., naringenin, genistein, LNnT): Must add new metabolites AND new reactions.
- **Precursor exists but flux ≈ 0** (e.g., UDP-GlcNAc baseline flux = 0 in LNnT case): Mark as a supply bottleneck even though the metabolite is in the model.

### Tools Required

| Tool | Function |
|------|----------|
| `search_model(model_name, keyword)` | Search metabolites and reactions by keyword |
| `query_gpr(model_name, gene_id)` | Look up gene-protein-reaction associations |
| `simulate_knockout(model, target, knockouts)` | Baseline FBA with optional knockouts |
| `change_media_and_simulate(model, carbon, O2, target)` | Baseline under different conditions |

### Output Schema

```json
{
  "existing_metabolites": ["phe__L_c", "malcoa_c", "amet_c"],
  "missing_metabolites": ["naringenin", "genistein", "biochanin_A"],
  "existing_reactions": [],
  "missing_reactions": ["PAL", "C4H", "4CL", "CHS", "CHI", "IFS", "HID", "IOMT"],
  "competing_pathways": [
    {"reaction": "PYRDC", "gene": "PDC1/5/6", "competes_for": "pyruvate"},
    {"reaction": "FAS*", "gene": "FAS1/FAS2", "competes_for": "malonyl-CoA"}
  ],
  "product_consuming_reactions": [
    {"reaction": "L_LACD2cm", "gene": "CYB2", "consumes": "L-lactate"}
  ],
  "baseline_fluxes": {"ACCOAC": 0.1145, "PPNDH": 0.0385, "METAT": 0.0054}
}
```

### Validated Performance

| Case | Existing | Missing | Key Competing Pathway |
|------|----------|---------|----------------------|
| L-Lactic acid | lac__L_c metabolite, no L-LDH reaction | 1 reaction | PYRDC → ethanol (15.8 flux) |
| LNnT | UDP-GlcNAc, UDP-Gal | 3 reactions + no lactose | Chitin synthase (flux ≈ 0) |
| Biochanin A | Phe, MalCoA, SAM | **8 reactions** | FAS consumes all MalCoA (0.11 flux) |

---

## Phase 3: Enzyme Sourcing

**Question**: "Where do we get the missing enzymes?"

**Goal**: For each missing reaction identified in Phase 2, find the best candidate enzyme from external organisms, with supporting biochemical data (kinetics, structure, protein family).

### Logic Flow

```
For each missing reaction (with known EC number from Phase 1):
  │
  ├─ uniprot_search(ec:{EC} AND reviewed:true)
  │    ├─ Rank candidates by priority:
  │    │    1. Same domain as host (eukaryotic > prokaryotic if host is yeast)
  │    │    2. Previously validated expression in host (literature evidence)
  │    │    3. Crystal structure available (PDB)
  │    │    4. Favorable kinetics (low Km, high kcat)
  │    └─ For top candidates → uniprot_entry(accession) → detailed info
  │
  ├─ Need a transporter? → UniProt keyword search (NOT KEGG)
  │    KEGG does not index transport proteins.
  │    Search: uniprot_search("lactate transporter" OR "formate nitrite transporter")
  │
  └─ Protein family analysis → interpro_query(IPR_id or PF_id)
       └─ Confirm domain architecture, transmembrane topology, species distribution
```

### Decision Points

- **Standard enzyme** (known EC, well-characterized): UniProt reviewed search is sufficient.
- **Transporter**: KEGG cannot find these. Must use UniProt keyword search.
  - L-Lactic acid: PfFNT (O77389) was invisible in KEGG, found via UniProt.
  - LNnT: LAC12 (P07921) found via UniProt gene name search.
- **Non-model organism enzymes**: UniProt `reviewed:true` filter biases toward model species (*G. max*, *A. thaliana*, *E. coli*). Enzymes from *P. lobata* or *P. falciparum* may only appear in unreviewed entries or require broader keyword searches.

### Tools Required

| Tool | Function | Status |
|------|----------|--------|
| `uniprot_search(query, max_results)` | Search UniProt by EC, gene name, organism | **Needs implementation** |
| `uniprot_entry(accession)` | Get detailed protein info (function, kinetics, PDB, GO) | **Needs implementation** |
| `interpro_query(entry_id)` | Get Pfam/InterPro domain info and species distribution | **Needs implementation** |
| KEGG orthology link API | Cross-species gene orthologs for a KO ID | Available via KEGG REST |

### Output Schema

```json
{
  "candidate_enzymes": [
    {
      "reaction": "IFS",
      "EC": "1.14.14.87",
      "uniprot_id": "Q9SWR5",
      "organism": "Glycine max",
      "gene": "IFS2",
      "length": 521,
      "Km_mM": null,
      "PDB": null,
      "notes": "CYP93C family, requires CPR"
    }
  ],
  "transporters": [
    {
      "uniprot_id": "O77389",
      "name": "PfFNT",
      "organism": "Plasmodium falciparum",
      "function": "lactate:proton symporter",
      "length": 309
    }
  ]
}
```

### Validated Performance

| Case | Enzyme Hits | Transporter Hits | Source Accuracy |
|------|-------------|-----------------|----------------|
| L-Lactic acid | L-LDH 30 candidates (EC 1.1.1.27) | PfFNT (O77389) ✓ | Function ✓, Source varies |
| LNnT | LgtA (Q9JXQ6), LgtB (Q51116) ✓ | LAC12 (P07921) ✓ | **100% match** |
| Biochanin A | All 8 enzymes found | — | Function ✓, Source differs (*G. max* vs *P. lobata*) |

### Known Limitations

- UniProt `reviewed:true` biases toward well-studied model organisms. The optimal enzyme may come from an understudied species.
- No automated way to assess "expression success in yeast" — this requires literature mining beyond current tool capabilities.

---

## Phase 4: Simulation & Optimization

**Question**: "Does the pathway work in silico? What are the bottlenecks? How do we optimize?"

**Goal**: Insert the heterologous pathway into the host GEM, validate feasibility, identify bottlenecks, and iteratively test engineering interventions — all within a single continuous FBA workflow.

This phase merges what were previously separate "Pathway Validation" and "Optimization" steps. The rationale: both use FBA on the same augmented model, and validation is simply the first iteration of optimization. Separating them creates an artificial boundary.

### Logic Flow

```
4.1 Construct & validate
  ├─ add_pathway(model, reactions[]) → verify feasibility
  ├─ maximize_product(model, target) → max theoretical flux
  └─ If max flux = 0 → debug: check metabolite names, compartments, connectivity

4.2 Quantify tradeoffs
  ├─ production_envelope(model, target, steps=10)
  │    → Pareto frontier: growth rate vs. product flux
  └─ Per-molecule cost analysis: how much Phe / MalCoA / SAM / NADPH per product?

4.3 Identify bottlenecks
  ├─ Metabolite balance: for each precursor/cofactor
  │    → list producing and consuming reactions, compute net flux
  │    → flag tightest supply
  └─ Identify by-products (ethanol, glycerol, CO2, etc.)

4.4 Eliminate competing pathways
  ├─ From Phase 2 competing_pathways → query_gpr → get gene(s)
  ├─ simulate_knockout(genes) → measure product flux change
  ├─ Compare: full KO vs. attenuation (set upper_bound)
  └─ Check: is KO lethal? (biomass = 0 → need attenuation instead)

4.5 Enhance precursor supply
  ├─ From 4.3 bottleneck metabolites → trace upstream pathway → find gene
  ├─ simulate_overexpression(gene, forced_lb) → test effect
  └─ Caution: forcing one enzyme too high can create new bottlenecks

4.6 Media & condition optimization
  ├─ O₂ level sweep (anaerobic → microaerobic → aerobic)
  ├─ Carbon source co-feeding tests (galactose, glycerol, etc.)
  └─ Quantify: which condition maximizes yield?

4.7 Transporter engineering
  ├─ search_model → native transporter exists?
  │    ├─ Exists → check directionality, consider overexpression
  │    └─ Missing → uniprot_search(transporter + substrate keyword)
  │                  ├─ Organic acids → MCT, FNT family
  │                  ├─ Sugars → sugar efflux, MFS family
  │                  └─ Others → ABC transporters
  └─ Add transporter reaction to model → re-run FBA → measure impact

4.8 Block product degradation
  ├─ Search reactions consuming the target product
  ├─ query_gpr → find gene → simulate_knockout
  └─ Check: is the reaction essential?
```

### Decision Points

- **Max flux = 0** (step 4.1): Substrates unreachable. Check metabolite naming, compartment tags, precursor connectivity.
- **Max flux > 0, biomass = 0** (step 4.2): Strong growth-production tradeoff. Consider two-stage fermentation.
- **Cofactor recycling adequate** (step 4.3): Not a bottleneck (e.g., UDP recycling via NDPK2 in LNnT). If inadequate → cofactor engineering.
- **KO is lethal** (step 4.4): Switch to attenuation (reduce upper bound) instead of full deletion.
- **Overexpression backfires** (step 4.5): Forcing flux too high can reduce product (observed in LNnT GFA1 case). Test incrementally.

### Tools Required

| Tool | Function | Status |
|------|----------|--------|
| `add_heterologous_reaction(model, id, rxn_string)` | Add single reaction + test | Available |
| `add_pathway(model, reactions[])` | Add multiple reactions + verify | **Needs MCP registration** |
| `maximize_product(model, target, min_biomass_frac)` | Max product flux with growth constraint | **Needs MCP registration** |
| `production_envelope(model, target, steps)` | Sweep biomass constraint | **Needs MCP registration** |
| `simulate_knockout(model, target, knockouts)` | Test gene deletions | Available |
| `simulate_overexpression(model, gene, forced_lb)` | Test forced minimum flux | Available |
| `change_media_and_simulate(model, carbon, O2, target)` | Test growth conditions | Available |
| `query_gpr(model, gene_id)` | Gene-reaction associations | Available |
| `search_model(model, keyword)` | Find metabolites/reactions | **Needs MCP registration** |
| `uniprot_search(query)` | Find heterologous transporters | **Needs implementation** |

### Output Schema

```json
{
  "validation": {
    "max_theoretical_flux": 0.789,
    "max_flux_units": "mmol/gDW/h",
    "production_envelope": [
      {"growth_fraction": 0.0, "biomass": 0.0, "product_flux": 0.789},
      {"growth_fraction": 0.5, "biomass": 0.144, "product_flux": 0.413},
      {"growth_fraction": 1.0, "biomass": 0.288, "product_flux": 0.0}
    ],
    "per_molecule_cost": {"Phe": 1, "MalCoA": 3, "SAM": 1, "NADPH": 2}
  },
  "bottlenecks": [
    {"metabolite": "malcoa_c", "demand": 2.14, "supply": 2.15, "margin": "0.5%"}
  ],
  "optimization_results": {
    "knockouts": [
      {"gene": "PDC1/5/6", "effect": "ethanol 15.8 → 0, lactate +15.8"}
    ],
    "overexpression": [
      {"gene": "GFA1", "effect": "UDP-GlcNAc supply +40%"}
    ],
    "media": [
      {"condition": "microaerobic O2=-0.5", "effect": "yield 79% → 85%"}
    ],
    "transporters": [
      {"protein": "PfFNT (O77389)", "effect": "unidirectional lactate export"}
    ]
  }
}
```

### Validated Performance

| Sub-step | L-Lactic acid | LNnT | Biochanin A |
|----------|---------------|------|-------------|
| 4.1–4.3 Validate + bottleneck | Max=20.0, NADH recycling | Max=2.12, UDP-GlcNAc | Max=0.789, MalCoA (3x) |
| 4.4 Competing KO | ✅ ΔPDC1/5/6 | — | ⚠️ FAS lethal, attenuate |
| 4.5 Precursor OE | — | ✅ GFA1/GNA1 | ✅ ACC1/ARO4 |
| 4.6 Media optimization | ✅ Microaerobic | ✅ Gal co-feed (+80%) | — |
| 4.7 Transporter | ✅ PfFNT | ✅ LAC12 | — |
| 4.8 Degradation block | ✅ ΔCYB2 | — | — |

---

## Optional: Sequence Engineering (Post-processing)

Codon optimization is a **standard, deterministic operation** — not a design decision. It is excluded from the core 4-phase reasoning loop for the following reasons:

1. **No design reasoning involved**: Given a CDS, codon optimization is a mechanical transformation (maximize CAI, remove restriction sites, preserve protein). Any molecular biologist does this routinely.
2. **Input dependency**: Requires full CDS nucleotide sequences from external databases (EMBL/GenBank), which the current tool chain does not automatically retrieve.
3. **Not a differentiator**: None of the three ground-truth papers mention codon optimization as part of their engineering strategy. It does not distinguish a good design from a bad one.

When needed, invoke `optimize_sequence(cds)` from `tools/dna_tool.py` (DNAchisel, targeting *S. cerevisiae*: CAI maximization + BsaI site removal + translation preservation).

---

## Tool Inventory

### Currently Available

| Tool | Location | MCP Registered | Used in Phase |
|------|----------|---------------|---------------|
| `search_kegg` | `tools/db_tools.py` | **No** | 1 |
| `query_kegg_compound` | `tools/db_tools.py` | **No** | 1 |
| `query_kegg_pathway` | `tools/db_tools.py` | **No** | 1 |
| `query_kegg_reaction` | `tools/db_tools.py` | **No** | 1 |
| `search_pubmed` | `tools/db_tools.py` | **No** | Auxiliary |
| `search_model` | `tools/fba_tool.py` | **No** | 2, 4 |
| `add_pathway` | `tools/fba_tool.py` | **No** | 4 |
| `maximize_product` | `tools/fba_tool.py` | **No** | 4 |
| `production_envelope` | `tools/fba_tool.py` | **No** | 4 |
| `simulate_knockout` | `tools/fba_tool.py` | **Yes** | 2, 4 |
| `simulate_overexpression` | `tools/fba_tool.py` | **Yes** | 4 |
| `change_media_and_simulate` | `tools/fba_tool.py` | **Yes** | 4 |
| `add_heterologous_reaction` | `tools/fba_tool.py` | **Yes** | 4 |
| `query_gpr` | `tools/fba_tool.py` | **Yes** | 2, 4 |
| `optimize_sequence` | `tools/dna_tool.py` | **Yes** | Optional |

### Needs Implementation

| Tool | API Endpoint | Used in Phase |
|------|-------------|---------------|
| `uniprot_search(query, max_results)` | `GET https://rest.uniprot.org/uniprotkb/search?query={query}&format=tsv` | 3 |
| `uniprot_entry(accession)` | `GET https://rest.uniprot.org/uniprotkb/{accession}?format=json` | 3 |
| `interpro_query(entry_id)` | `GET https://www.ebi.ac.uk/interpro/api/entry/{db}/{id}?format=json` | 3 |

### Needs MCP Registration (code exists, not exposed)

- `search_kegg`, `query_kegg_compound`, `query_kegg_pathway`, `query_kegg_reaction`, `search_pubmed`
- `search_model`, `add_pathway`, `maximize_product`, `production_envelope`

---

## System Prompt Design

The agent's behavior is driven by a structured system prompt that encodes the 4-phase pipeline:

```
You are a metabolic engineering design agent. When the user provides
a (target product, host organism) pair, follow this 4-phase workflow:

Phase 1 [Target ID]: Use KEGG tools to search for the target compound,
  find its compound ID, biosynthetic pathway (pathway/module), key
  reactions, and enzymes (EC/KO). If KEGG has no direct entry, infer
  the pathway from structural analogs or homologous pathways.

Phase 2 [Host Assessment]: Use search_model to check which Phase 1
  metabolites/reactions exist in the host GEM. Run baseline FBA to
  understand current flux distribution. Identify competing pathways
  and product-consuming reactions.

Phase 3 [Enzyme Sourcing]: For each missing reaction, use
  uniprot_search to find candidate enzymes by EC number. Use
  uniprot_entry for kinetics and structure details. IMPORTANT:
  transporters must be searched via UniProt — KEGG does not index them.

Phase 4 [Simulation & Optimization]: This is a continuous FBA workflow:
  4.1 Add heterologous reactions to the GEM and verify feasibility.
  4.2 Compute maximum theoretical flux and production envelope.
  4.3 Perform metabolite balance analysis to identify bottlenecks.
  4.4 Test gene knockouts to eliminate competing pathways.
  4.5 Test overexpression to enhance precursor supply.
  4.6 Sweep media conditions (O₂, carbon source co-feeding).
  4.7 Evaluate transporter engineering for product export.
  4.8 Block product degradation reactions.
  Iterate between sub-steps as needed — each intervention may reveal
  new bottlenecks that require re-analysis.

After completing each phase, summarize findings before proceeding.
Output the complete engineering strategy and reasoning chain at the end.

If the user requests sequence engineering (codon optimization), use
optimize_sequence as a post-processing step.
```

---

## Validated Capabilities and Known Boundaries

### Capabilities (verified across 3 cases)

| Capability | Evidence |
|------------|----------|
| Pathway discovery from product name | 3/3 cases: correct core pathway identified |
| Host metabolic assessment | 3/3 cases: correct existing/missing classification |
| Cross-species enzyme sourcing | 3/3 cases: correct EC/function; source species may vary |
| Quantitative feasibility analysis | 3/3 cases: max flux, production envelope, bottleneck ID |
| Multi-type optimization strategies | KO, OE, media, transporter strategies all validated |
| Reasoning chain traceability | Every prediction traceable to specific tool call |

### Known Boundaries

| Limitation | Exposed By | Root Cause | Potential Solution |
|------------|-----------|------------|-------------------|
| KEGG returns canonical (natural) pathways only | Biochanin A: missed TAL shortcut | KEGG is a biochemistry DB, not an engineering DB | Engineering shortcut knowledge base |
| UniProt `reviewed` bias toward model species | Biochanin A: *P. lobata* enzymes missed | Curation bias | Include unreviewed entries + literature mining |
| Cannot model CYP450 / membrane protein engineering | Biochanin A: missed ER expansion, heme, VHb, chaperone | FBA treats enzymes as infinite-capacity catalysts | CYP450-specific evaluation module |
| FBA cannot distinguish attenuation vs. full KO | L-Lactic acid: PDC attenuation vs. KO | Steady-state stoichiometric model limitation | Enzyme-constrained models (GECKO/sMOMENT) |
| Cannot predict protein engineering strategies | LNnT: missed peptide tag assembly | No protein-protein interaction modeling | Protein structure/interaction tools (PDB, STRING) |
| Process engineering is out of scope | L-Lactic acid: ALE; LNnT: fed-batch | No bioprocess simulation tools | Integrate dFBA or process simulation |
| Cannot predict copy number / expression level effects | Biochanin A: IFS copy number optimization | FBA has no gene expression layer | Integrate expression models or ML prediction |

---

## Appendix: Validation Case Summaries

### Case 1: L-Lactic Acid

- **Pathway**: Pyruvate + NADH → L-Lactate + NAD⁺ (1 step, L-LDH EC 1.1.1.27)
- **Key predictions matched**: L-LDH enzyme, ΔPDC1/5/6, ΔCYB2, microaerobic, PfFNT transporter
- **Key miss**: ALE for acid tolerance, PDC attenuation (vs. KO)
- **Consistency**: ~85% (after supplementary tool round)

### Case 2: Lacto-N-neotetraose (LNnT)

- **Pathway**: Lactose → (LgtA + UDP-GlcNAc) → LNT-II → (LgtB + UDP-Gal) → LNnT (3 steps)
- **Key predictions matched**: LgtA (Q9JXQ6), LgtB (Q51116), LAC12 (P07921), GNA1/AGM1/UAP1, Gal10
- **Key miss**: Peptide tag enzyme assembly, fed-batch optimization
- **Consistency**: ~85%

### Case 3: Biochanin A

- **Pathway**: Phe → Cinnamate → p-Coumarate → p-Coumaroyl-CoA → Naringenin → Genistein → Biochanin A (8 steps)
- **Key predictions matched**: Core 7-step pathway, all enzyme functions (EC), CHS/CHI/HID sources
- **Key miss**: TAL shortcut, *P. lobata* enzyme source, organelle engineering (ER expansion, VHb, heme, chaperones)
- **Consistency**: ~60%
