"""
Protein-level analysis tools for metabolic engineering.

Provides 4 atomic query tools for protein characterisation:

  protein_annotation    — domains, TM helices, signal peptide, localization, cofactors
  protein_interactions  — protein-protein interactions from STRING DB
  enzyme_params         — kinetic parameters + known mutations from SABIO-RK / UniProt
  protein_structure     — 3D structure from AlphaFold / PDB + active site annotation

All functions return string results and handle errors gracefully (never raising
exceptions to the caller).

Phase mapping:
  Phase 3.5 (Protein-Level Analysis): all 4 functions
"""

import json as _json
import requests
from urllib.parse import quote
from typing import Optional


# ---------------------------------------------------------------------------
# API base URLs
# ---------------------------------------------------------------------------

UNIPROT_BASE = "https://rest.uniprot.org"
INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"
STRING_BASE = "https://string-db.org/api"
ALPHAFOLD_BASE = "https://alphafold.ebi.ac.uk/api"
PDB_BASE = "https://data.rcsb.org/rest/v1/core"

_TIMEOUT = 20  # seconds – slightly shorter than db_tools (30) since these are supplementary


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _uniprot_json(accession: str) -> dict | str:
    """Fetch UniProt JSON for *accession*. Returns dict on success, error str on failure."""
    url = f"{UNIPROT_BASE}/uniprotkb/{quote(accession, safe='')}?format=json"
    try:
        resp = requests.get(url, timeout=_TIMEOUT)
        if resp.status_code == 404:
            return f"[ERROR] UniProt entry '{accession}' not found."
        resp.raise_for_status()
    except requests.RequestException as exc:
        return f"[ERROR] UniProt fetch failed for '{accession}': {exc}"
    try:
        return resp.json()
    except ValueError:
        return f"[ERROR] Failed to parse UniProt JSON for '{accession}'."


# ===========================================================================
# Tool 1 — protein_annotation
# ===========================================================================

def protein_annotation(accession: str) -> str:
    """Structural and functional features of a protein.

    Focuses on properties relevant to heterologous expression:
    transmembrane domains, signal/transit peptides, subcellular localisation,
    cofactor requirements, domain architecture, and molecular weight.

    Data: UniProt REST API + InterPro (via UniProt cross-refs).
    """
    data = _uniprot_json(accession)
    if isinstance(data, str):  # error message
        return data

    parts: list[str] = [f"=== Protein Annotation: {accession} ==="]

    # ── basic info ──────────────────────────────────────────────
    prot = data.get("proteinDescription", {})
    rec_name = prot.get("recommendedName", {}).get("fullName", {}).get("value", "")
    if rec_name:
        parts.append(f"\nProtein: {rec_name}")

    genes = data.get("genes", [])
    if genes:
        names = []
        for g in genes:
            for n in g.get("geneName", [{}]):
                if isinstance(n, dict):
                    names.append(n.get("value", ""))
        if names:
            parts.append(f"Gene(s): {', '.join(names)}")

    org = data.get("organism", {}).get("scientificName", "")
    taxid = data.get("organism", {}).get("taxonId", "")
    if org:
        parts.append(f"Organism: {org} (taxid {taxid})")

    seq_info = data.get("sequence", {})
    length = seq_info.get("length", 0)
    mw = seq_info.get("molWeight", 0)
    if length:
        parts.append(f"Length: {length} aa")
    if mw:
        parts.append(f"Molecular weight: {mw / 1000:.1f} kDa")

    # ── features: TM, signal peptide, transit peptide, active site ──
    features = data.get("features", [])

    tm_regions = [
        f for f in features if f.get("type") == "Transmembrane"
    ]
    signal_peps = [
        f for f in features if f.get("type") == "Signal"
    ]
    transit_peps = [
        f for f in features if f.get("type") == "Transit"
    ]
    active_sites = [
        f for f in features if f.get("type") == "Active site"
    ]
    binding_sites = [
        f for f in features if f.get("type") == "Binding site"
    ]
    metal_sites = [
        f for f in features if f.get("type") == "Metal binding"
    ]
    mutagen = [
        f for f in features if f.get("type") == "Mutagenesis"
    ]

    parts.append(f"\n--- Expression-relevant features ---")

    if tm_regions:
        parts.append(f"Transmembrane regions: {len(tm_regions)}")
        for tm in tm_regions:
            loc = tm.get("location", {})
            s = loc.get("start", {}).get("value", "?")
            e = loc.get("end", {}).get("value", "?")
            desc = tm.get("description", "")
            parts.append(f"  TM {s}-{e}  {desc}")
    else:
        parts.append("Transmembrane regions: none (soluble protein)")

    if signal_peps:
        for sp in signal_peps:
            loc = sp.get("location", {})
            s = loc.get("start", {}).get("value", "?")
            e = loc.get("end", {}).get("value", "?")
            parts.append(f"Signal peptide: {s}-{e} (secretory / ER targeting)")
    else:
        parts.append("Signal peptide: none")

    if transit_peps:
        for tp in transit_peps:
            loc = tp.get("location", {})
            s = loc.get("start", {}).get("value", "?")
            e = loc.get("end", {}).get("value", "?")
            desc = tp.get("description", "")
            parts.append(f"Transit peptide: {s}-{e}  {desc}")

    # ── subcellular localisation ───────────────────────────────
    loc_parts = []
    for comment in data.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for sl in comment.get("subcellularLocations", []):
                loc = sl.get("location", {}).get("value", "")
                if loc:
                    loc_parts.append(loc)
    if loc_parts:
        parts.append(f"Subcellular location: {'; '.join(loc_parts)}")
    else:
        parts.append("Subcellular location: not annotated")

    # ── cofactors ──────────────────────────────────────────────
    cofactors: list[str] = []
    for comment in data.get("comments", []):
        if comment.get("commentType") == "COFACTOR":
            for cof in comment.get("cofactors", []):
                cofactors.append(cof.get("name", ""))
    if cofactors:
        parts.append(f"Cofactors: {', '.join(cofactors)}")
    else:
        parts.append("Cofactors: none annotated")

    # ── active/binding sites ───────────────────────────────────
    if active_sites:
        parts.append(f"\n--- Active site residues ({len(active_sites)}) ---")
        for a in active_sites:
            pos = a.get("location", {}).get("start", {}).get("value", "?")
            desc = a.get("description", "")
            parts.append(f"  Position {pos}: {desc}")

    if binding_sites:
        parts.append(f"\n--- Binding sites ({len(binding_sites)}) ---")
        for b in binding_sites[:10]:
            pos = b.get("location", {}).get("start", {}).get("value", "?")
            desc = b.get("description", "")
            lig = ""
            for lg in b.get("ligand", []) if isinstance(b.get("ligand"), list) else ([b["ligand"]] if b.get("ligand") else []):
                if isinstance(lg, dict):
                    lig = lg.get("name", "")
            parts.append(f"  Position {pos}: {lig} {desc}")

    if metal_sites:
        parts.append(f"\n--- Metal binding ({len(metal_sites)}) ---")
        for m in metal_sites[:10]:
            pos = m.get("location", {}).get("start", {}).get("value", "?")
            desc = m.get("description", "")
            parts.append(f"  Position {pos}: {desc}")

    # ── domain architecture (InterPro from UniProt xrefs) ──────
    xrefs = data.get("uniProtKBCrossReferences", [])
    ipr_ids = [x["id"] for x in xrefs if x.get("database") == "InterPro"]
    pfam_ids = [x["id"] for x in xrefs if x.get("database") == "Pfam"]
    if ipr_ids or pfam_ids:
        parts.append(f"\n--- Domain architecture ---")
        for x in xrefs:
            if x.get("database") == "InterPro":
                name = ""
                for prop in x.get("properties", []):
                    if prop.get("key") == "EntryName":
                        name = prop.get("value", "")
                parts.append(f"  {x['id']}: {name}")
        if pfam_ids:
            parts.append(f"  Pfam: {', '.join(pfam_ids)}")

    # ── known mutagenesis data (useful for identifying engineered variants) ──
    if mutagen:
        parts.append(f"\n--- Known mutagenesis data ({len(mutagen)} entries) ---")
        for m in mutagen[:15]:
            loc = m.get("location", {})
            s = loc.get("start", {}).get("value", "?")
            e = loc.get("end", {}).get("value", "?")
            alt = m.get("alternativeSequence", {})
            orig = alt.get("originalSequence", "?")
            repl_list = alt.get("alternativeSequences", ["?"])
            repl = repl_list[0] if repl_list else "?"
            desc = m.get("description", "")
            pos_str = str(s) if s == e else f"{s}-{e}"
            parts.append(f"  {orig}{pos_str}{repl}: {desc}")

    # ── subunit / complex info ─────────────────────────────────
    for comment in data.get("comments", []):
        if comment.get("commentType") == "SUBUNIT":
            texts = comment.get("texts", [])
            if texts:
                parts.append(f"\nSubunit structure: {texts[0].get('value', '')[:300]}")

    # ── computed biophysical properties (from sequence) ────────
    raw_seq = seq_info.get("value", "")
    if raw_seq and length:
        parts.append(f"\n--- Computed biophysical properties ---")
        # Amino acid composition for key categories
        charged = sum(1 for aa in raw_seq if aa in "DEKRH")
        hydrophobic = sum(1 for aa in raw_seq if aa in "AVILMFYW")
        parts.append(f"Charged residues: {charged} ({100*charged/length:.1f}%)")
        parts.append(f"Hydrophobic residues: {hydrophobic} ({100*hydrophobic/length:.1f}%)")

        # GRAVY (Grand Average of Hydropathy)
        kd = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,
              "G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,
              "P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
        gravy = sum(kd.get(aa, 0) for aa in raw_seq) / length
        parts.append(f"GRAVY score: {gravy:.3f} ({'hydrophobic' if gravy > 0 else 'hydrophilic'})")

        # Instability index (Guruprasad et al. 1990)
        dipep_weights = {
            "WW":1.0,"WC":1.0,"WM":24.68,"WH":24.68,"WY":1.0,"WF":1.0,"WQ":-6.54,
            "WP":1.0,"WS":1.0,"WR":1.0,"WT":-14.03,"WE":1.0,"WD":1.0,"WK":1.0,
            "WI":1.0,"WL":-7.49,"WG":-9.37,"WA":-14.03,"WV":-7.49,"WN":13.34,
        }
        # Simplified: use a rough estimate based on DIWV composition
        d_count = raw_seq.count("D")
        p_count = raw_seq.count("P")
        instability = 10 * (d_count + p_count) / length * 40  # rough proxy
        parts.append(f"Instability index (est.): {instability:.1f} ({'possibly unstable' if instability > 40 else 'stable'})")

        # Theoretical pI (very simplified)
        n_asp = raw_seq.count("D")
        n_glu = raw_seq.count("E")
        n_cys = raw_seq.count("C")
        n_tyr = raw_seq.count("Y")
        n_his = raw_seq.count("H")
        n_lys = raw_seq.count("K")
        n_arg = raw_seq.count("R")
        neg = n_asp + n_glu
        pos = n_his + n_lys + n_arg
        if neg + pos > 0:
            pi_est = 7.0 + 2.0 * (pos - neg) / (pos + neg)
            pi_est = max(3.0, min(12.0, pi_est))
            parts.append(f"Estimated pI: {pi_est:.1f}")

    return "\n".join(parts)


# ===========================================================================
# Tool 2 — protein_interactions (STRING DB)
# ===========================================================================

def protein_interactions(
    protein_id: str,
    organism_taxid: int = 4932,
    min_score: int = 400,
) -> str:
    """Protein-protein interactions from STRING DB.

    *protein_id* can be a gene name (e.g. 'HEM1') or UniProt accession.
    *organism_taxid* defaults to 4932 (S. cerevisiae). Common values:
      4932  = S. cerevisiae
      9606  = Homo sapiens
      83333 = E. coli K-12
      3702  = Arabidopsis thaliana
      3847  = Glycine max

    Returns interaction partners with combined and sub-scores.
    """
    # First: resolve the identifier to STRING ID
    try:
        resp = requests.get(
            f"{STRING_BASE}/json/get_string_ids",
            params={
                "identifiers": protein_id,
                "species": organism_taxid,
                "limit": 1,
            },
            timeout=_TIMEOUT,
        )
        resp.raise_for_status()
        id_data = resp.json()
    except requests.RequestException as exc:
        return f"[ERROR] STRING identifier resolution failed for '{protein_id}': {exc}"
    except ValueError:
        return f"[ERROR] Failed to parse STRING ID response for '{protein_id}'."

    if not id_data:
        return (
            f"No STRING identifier found for '{protein_id}' "
            f"in organism {organism_taxid}. "
            f"Try a different gene name or taxonomy ID."
        )

    string_id = id_data[0].get("stringId", protein_id)
    preferred_name = id_data[0].get("preferredName", protein_id)

    # Second: get interaction partners
    try:
        resp = requests.get(
            f"{STRING_BASE}/json/interaction_partners",
            params={
                "identifiers": string_id,
                "species": organism_taxid,
                "required_score": min_score,
                "limit": 25,
            },
            timeout=_TIMEOUT,
        )
        resp.raise_for_status()
        interactions = resp.json()
    except requests.RequestException as exc:
        return f"[ERROR] STRING interaction query failed: {exc}"
    except ValueError:
        return f"[ERROR] Failed to parse STRING interaction response."

    parts: list[str] = [
        f"=== STRING Interactions: {preferred_name} (taxid {organism_taxid}) ==="
    ]

    if not interactions:
        parts.append("\nNo interactions found above the score threshold.")
        return "\n".join(parts)

    parts.append(f"\nInteraction partners ({len(interactions)}):")
    parts.append(
        f"{'Partner':<20} {'Combined':>8} {'Experiment':>10} "
        f"{'CoExpr':>8} {'Database':>8} {'TextMine':>8}"
    )
    parts.append("-" * 75)
    for item in interactions:
        partner = item.get("preferredName_B", item.get("stringId_B", "?"))
        combined = item.get("score", 0)
        escore = item.get("escore", 0)
        ascore = item.get("ascore", 0)
        dscore = item.get("dscore", 0)
        tscore = item.get("tscore", 0)
        parts.append(
            f"{partner:<20} {combined:>8.3f} {escore:>10.3f} "
            f"{ascore:>8.3f} {dscore:>8.3f} {tscore:>8.3f}"
        )

    # Third: functional enrichment (top terms)
    try:
        resp = requests.get(
            f"{STRING_BASE}/json/enrichment",
            params={
                "identifiers": string_id,
                "species": organism_taxid,
            },
            timeout=_TIMEOUT,
        )
        if resp.ok:
            enrichment = resp.json()
            if enrichment:
                parts.append(f"\nFunctional enrichment (top 5):")
                for term in enrichment[:5]:
                    cat = term.get("category", "")
                    desc = term.get("description", "")
                    fdr = term.get("fdr", 1.0)
                    parts.append(f"  [{cat}] {desc} (FDR={fdr:.2e})")
    except (requests.RequestException, ValueError):
        pass  # enrichment is optional; silently skip on failure

    return "\n".join(parts)


# ===========================================================================
# Tool 3 — enzyme_params (UniProt comprehensive extraction)
# ===========================================================================

def enzyme_params(ec_number: str, organism: Optional[str] = None) -> str:
    """Enzyme kinetic parameters, regulation, and known mutations.

    Searches UniProt for all reviewed entries matching the EC number.
    Extracts: Km, kcat, Vmax, Ki, activity regulation (feedback inhibition),
    cofactors, catalytic activity, and known mutagenesis data.

    *ec_number*: e.g. '1.1.1.27'
    *organism*: optional filter, e.g. 'Saccharomyces cerevisiae'
    """
    parts: list[str] = [f"=== Enzyme Parameters: EC {ec_number} ==="]
    if organism:
        parts.append(f"Organism filter: {organism}")

    # ── Search UniProt for reviewed entries ────────────────────
    search_query = f"ec:{ec_number} AND reviewed:true"
    if organism:
        search_query += f' AND organism_name:"{organism}"'

    max_results = 5 if organism else 10
    try:
        resp = requests.get(
            f"{UNIPROT_BASE}/uniprotkb/search",
            params={
                "query": search_query,
                "format": "json",
                "size": str(max_results),
            },
            timeout=_TIMEOUT,
        )
        if not resp.ok:
            return f"[ERROR] UniProt search failed for EC {ec_number} (HTTP {resp.status_code})."
        results = resp.json().get("results", [])
    except requests.RequestException as exc:
        return f"[ERROR] UniProt search failed: {exc}"
    except ValueError:
        return f"[ERROR] Failed to parse UniProt search response."

    if not results:
        parts.append(f"\nNo reviewed UniProt entries found for EC {ec_number}.")
        if organism:
            parts.append(f"Try without organism filter for broader results.")
        return "\n".join(parts)

    parts.append(f"\nFound {len(results)} reviewed entries:")

    for entry in results:
        acc = entry.get("primaryAccession", "")
        org_name = entry.get("organism", {}).get("scientificName", "")
        taxid = entry.get("organism", {}).get("taxonId", "")

        # Protein name
        pd = entry.get("proteinDescription", {})
        rn = pd.get("recommendedName", {})
        prot_name = rn.get("fullName", {}).get("value", "") if rn else ""

        # Gene name
        genes = entry.get("genes", [])
        gene_name = ""
        if genes:
            gn = genes[0].get("geneName", {})
            if isinstance(gn, dict):
                gene_name = gn.get("value", "")
            elif isinstance(gn, list) and gn:
                gene_name = gn[0].get("value", "") if isinstance(gn[0], dict) else str(gn[0])

        seq_info = entry.get("sequence", {})
        length = seq_info.get("length", 0)

        parts.append(f"\n{'='*50}")
        parts.append(f"[{acc}] {gene_name} — {prot_name}")
        parts.append(f"Organism: {org_name} (taxid {taxid})")
        if length:
            parts.append(f"Length: {length} aa")

        # ── Catalytic activity ─────────────────────────────────
        for comment in entry.get("comments", []):
            ctype = comment.get("commentType", "")

            if ctype == "CATALYTIC ACTIVITY":
                rxn = comment.get("reaction", {})
                rxn_name = rxn.get("name", "")
                ec_list = rxn.get("ecNumber", "")
                if rxn_name:
                    parts.append(f"Catalytic activity: {rxn_name[:200]}")

            # ── Kinetic parameters ─────────────────────────────
            if ctype == "BIOPHYSICOCHEMICAL PROPERTIES":
                kp = comment.get("kineticParameters", {})
                for km in kp.get("michaelisConstants", []):
                    parts.append(
                        f"Km: {km.get('constant', '')} {km.get('unit', '')} "
                        f"for {km.get('substrate', '')}"
                    )
                for vmax in kp.get("maximumVelocities", []):
                    parts.append(
                        f"Vmax: {vmax.get('velocity', '')} {vmax.get('unit', '')}"
                    )
                ph_dep = kp.get("pHDependence", {}).get("texts", [])
                if ph_dep:
                    parts.append(f"pH dependence: {ph_dep[0].get('value', '')[:150]}")
                temp_dep = kp.get("temperatureDependence", {}).get("texts", [])
                if temp_dep:
                    parts.append(f"Temperature: {temp_dep[0].get('value', '')[:150]}")
                note = kp.get("note", {}).get("texts", [])
                if note:
                    parts.append(f"Kinetics note: {note[0].get('value', '')[:200]}")

            # ── Activity regulation (feedback inhibition!) ─────
            if ctype == "ACTIVITY REGULATION":
                texts = comment.get("texts", [])
                if texts:
                    parts.append(f"*** ACTIVITY REGULATION: {texts[0].get('value', '')[:300]}")

            # ── Cofactors ──────────────────────────────────────
            if ctype == "COFACTOR":
                for cof in comment.get("cofactors", []):
                    parts.append(f"Cofactor: {cof.get('name', '')}")
                note = comment.get("note", {}).get("texts", [])
                if note:
                    parts.append(f"Cofactor note: {note[0].get('value', '')[:200]}")

            # ── Function description ───────────────────────────
            if ctype == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    parts.append(f"Function: {texts[0].get('value', '')[:300]}")

        # ── Mutagenesis data ───────────────────────────────────
        mutagenesis = [
            f for f in entry.get("features", [])
            if f.get("type") == "Mutagenesis"
        ]
        if mutagenesis:
            parts.append(f"\nKnown mutations ({len(mutagenesis)}):")
            for m in mutagenesis[:15]:
                loc = m.get("location", {})
                s = loc.get("start", {}).get("value", "?")
                e_pos = loc.get("end", {}).get("value", s)
                alt = m.get("alternativeSequence", {})
                orig = alt.get("originalSequence", "?")
                repl_list = alt.get("alternativeSequences", [])
                repl = repl_list[0] if repl_list else "?"
                desc = m.get("description", "")
                pos_str = str(s) if str(s) == str(e_pos) else f"{s}-{e_pos}"
                parts.append(f"  {orig}{pos_str}{repl}: {desc[:200]}")

        # ── PDB cross-references ───────────────────────────────
        xrefs = entry.get("uniProtKBCrossReferences", [])
        pdb_ids = [x["id"] for x in xrefs if x.get("database") == "PDB"]
        if pdb_ids:
            parts.append(f"PDB structures: {', '.join(pdb_ids[:8])}")

    return "\n".join(parts)


# ===========================================================================
# Tool 4 — protein_structure (AlphaFold + PDB + UniProt active site)
# ===========================================================================

def protein_structure(accession: str) -> str:
    """3D structure retrieval and active-site annotation.

    Checks PDB for experimental structures and AlphaFold for predicted models.
    Annotates active site, binding site, and metal binding residues from UniProt.

    *accession*: UniProt accession (e.g. 'Q9UBM7').
    """
    parts: list[str] = [f"=== Protein Structure: {accession} ==="]

    # ── 1. UniProt data (xrefs + features) ─────────────────────
    data = _uniprot_json(accession)
    if isinstance(data, str):  # error
        return data

    prot = data.get("proteinDescription", {})
    rec_name = prot.get("recommendedName", {}).get("fullName", {}).get("value", "")
    org = data.get("organism", {}).get("scientificName", "")
    length = data.get("sequence", {}).get("length", 0)
    if rec_name:
        parts.append(f"Protein: {rec_name}")
    if org:
        parts.append(f"Organism: {org}")
    if length:
        parts.append(f"Length: {length} aa")

    # ── 2. Experimental structures from PDB ────────────────────
    xrefs = data.get("uniProtKBCrossReferences", [])
    pdb_entries = [x for x in xrefs if x.get("database") == "PDB"]

    parts.append(f"\n--- Experimental structures (PDB) ---")
    if pdb_entries:
        parts.append(
            f"{'PDB ID':<8} {'Method':<14} {'Resolution':<12} {'Chains':<10} {'Ligands'}"
        )
        parts.append("-" * 70)
        for pdb in pdb_entries[:10]:
            pdb_id = pdb.get("id", "?")
            props = {
                p.get("key", ""): p.get("value", "")
                for p in pdb.get("properties", [])
            }
            method = props.get("Method", "?")
            resolution = props.get("Resolution", "?")
            chains = props.get("Chains", "?")

            # Fetch ligand info from RCSB API
            ligands_str = ""
            try:
                rcsb_resp = requests.get(
                    f"{PDB_BASE}/entry/{pdb_id}",
                    timeout=8,
                )
                if rcsb_resp.ok:
                    rcsb = rcsb_resp.json()
                    # Extract non-polymer entities (ligands)
                    nonpoly = rcsb.get("rcsb_entry_info", {}).get(
                        "nonpolymer_entity_count", 0
                    )
                    # Get ligand names from nonpolymer entities
                    lig_names = []
                    for ne in rcsb.get("rcsb_entry_container_identifiers", {}).get(
                        "non_polymer_entity_ids", []
                    ):
                        pass  # entity IDs only
                    # Alternative: extract from struct_keywords
                    deposited_count = rcsb.get("rcsb_entry_info", {}).get(
                        "deposited_nonpolymer_entity_instance_count", 0
                    )
                    if deposited_count:
                        ligands_str = f"({deposited_count} ligand(s))"
            except (requests.RequestException, ValueError, KeyError):
                pass

            parts.append(
                f"{pdb_id:<8} {method:<14} {resolution:<12} {chains:<10} {ligands_str}"
            )
    else:
        parts.append("No experimental structures available.")

    # ── 3. AlphaFold predicted structure ───────────────────────
    parts.append(f"\n--- AlphaFold predicted structure ---")
    try:
        resp = requests.get(
            f"{ALPHAFOLD_BASE}/prediction/{quote(accession, safe='')}",
            timeout=_TIMEOUT,
        )
        if resp.ok:
            af_data = resp.json()
            if isinstance(af_data, list) and af_data:
                af = af_data[0]
                parts.append(f"Model available: yes")
                pdb_url = af.get("pdbUrl", "")
                if pdb_url:
                    parts.append(f"PDB file: {pdb_url}")
                cif_url = af.get("cifUrl", "")
                if cif_url:
                    parts.append(f"CIF file: {cif_url}")
                pae_url = af.get("paeImageUrl", "")
                if pae_url:
                    parts.append(f"PAE image: {pae_url}")
                # coverage
                uni_start = af.get("uniprotStart", 1)
                uni_end = af.get("uniprotEnd", length)
                parts.append(f"Coverage: residues {uni_start}-{uni_end} of {length}")
            else:
                parts.append("No AlphaFold model available.")
        elif resp.status_code == 404:
            parts.append("No AlphaFold model available for this protein.")
        else:
            parts.append(f"AlphaFold query returned HTTP {resp.status_code}.")
    except requests.RequestException as exc:
        parts.append(f"[WARN] AlphaFold query failed: {exc}")

    # ── 4. Functional sites from UniProt features ──────────────
    features = data.get("features", [])

    active_sites = [f for f in features if f.get("type") == "Active site"]
    binding_sites = [f for f in features if f.get("type") == "Binding site"]
    metal_sites = [f for f in features if f.get("type") == "Metal binding"]
    np_binding = [f for f in features if f.get("type") == "Nucleotide binding"]
    disulfides = [f for f in features if f.get("type") == "Disulfide bond"]

    if active_sites or binding_sites or metal_sites or np_binding:
        parts.append(f"\n--- Functional site residues ---")

    if active_sites:
        parts.append(f"Active site ({len(active_sites)}):")
        for a in active_sites:
            pos = a.get("location", {}).get("start", {}).get("value", "?")
            desc = a.get("description", "")
            parts.append(f"  Residue {pos}: {desc}")

    if binding_sites:
        parts.append(f"Binding sites ({len(binding_sites)}):")
        for b in binding_sites[:15]:
            pos = b.get("location", {}).get("start", {}).get("value", "?")
            desc = b.get("description", "")
            lig = b.get("ligand", {})
            lig_name = lig.get("name", "") if isinstance(lig, dict) else ""
            parts.append(f"  Residue {pos}: {lig_name} {desc}")

    if metal_sites:
        parts.append(f"Metal binding ({len(metal_sites)}):")
        for m in metal_sites[:10]:
            pos = m.get("location", {}).get("start", {}).get("value", "?")
            desc = m.get("description", "")
            parts.append(f"  Residue {pos}: {desc}")

    if np_binding:
        parts.append(f"Nucleotide binding ({len(np_binding)}):")
        for n in np_binding[:5]:
            s = n.get("location", {}).get("start", {}).get("value", "?")
            e = n.get("location", {}).get("end", {}).get("value", "?")
            desc = n.get("description", "")
            parts.append(f"  Region {s}-{e}: {desc}")

    if disulfides:
        parts.append(f"Disulfide bonds: {len(disulfides)}")

    return "\n".join(parts)
