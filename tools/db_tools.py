"""
External database query tools for pathway and enzyme discovery.

Provides functions to query KEGG, PubMed, UniProt, and InterPro for pathway,
compound, reaction, enzyme, and protein family information. All functions
return string results and handle errors gracefully (never raising exceptions
to the caller).

Phase mapping (see Framework Design.md):
  Phase 1 (Target ID):     search_kegg, query_kegg_compound/pathway/reaction
  Phase 3 (Enzyme Sourcing): uniprot_search, uniprot_entry, interpro_entry
  Auxiliary:                search_pubmed
"""

import json as _json
import requests
from urllib.parse import quote


# ---------------------------------------------------------------------------
# KEGG helpers
# ---------------------------------------------------------------------------

KEGG_BASE = "https://rest.kegg.jp"


def _kegg_get(entry_id: str) -> str:
    """Fetch a KEGG entry by ID. Returns raw text or an error message string."""
    url = f"{KEGG_BASE}/get/{quote(entry_id, safe='')}"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 404:
            return f"[ERROR] KEGG entry '{entry_id}' not found (HTTP 404)."
        resp.raise_for_status()
        return resp.text
    except requests.RequestException as exc:
        return f"[ERROR] Network error querying KEGG for '{entry_id}': {exc}"


def _parse_kegg_sections(raw: str) -> dict[str, list[str]]:
    """Parse flat-file KEGG text into {SECTION_NAME: [lines ...]}."""
    sections: dict[str, list[str]] = {}
    current_key: str | None = None
    for line in raw.splitlines():
        if line.startswith("///"):
            break
        # Section headers start at column 0 with uppercase letters
        if line and line[0] != " ":
            parts = line.split(None, 1)
            current_key = parts[0]
            rest = parts[1] if len(parts) > 1 else ""
            sections.setdefault(current_key, [])
            if rest:
                sections[current_key].append(rest.strip())
        elif current_key is not None:
            sections[current_key].append(line.strip())
    return sections


# ---------------------------------------------------------------------------
# 1. query_kegg_pathway
# ---------------------------------------------------------------------------

def query_kegg_pathway(pathway_id: str) -> str:
    """Query KEGG for a pathway or module entry.

    Supports pathway IDs like 'map00940' and module IDs like 'M00941'.

    Returns a formatted text string with NAME, DESCRIPTION, COMPOUND,
    REACTION, and MODULE information, or an error message on failure.
    """
    raw = _kegg_get(pathway_id)
    if raw.startswith("[ERROR]"):
        return raw

    sections = _parse_kegg_sections(raw)
    parts: list[str] = [f"=== KEGG Entry: {pathway_id} ==="]

    # NAME
    if "NAME" in sections:
        parts.append(f"\nNAME:\n  " + "\n  ".join(sections["NAME"]))

    # DESCRIPTION
    if "DESCRIPTION" in sections:
        parts.append(f"\nDESCRIPTION:\n  " + "\n  ".join(sections["DESCRIPTION"]))

    # DEFINITION (modules use this instead of DESCRIPTION)
    if "DEFINITION" in sections:
        parts.append(f"\nDEFINITION:\n  " + "\n  ".join(sections["DEFINITION"]))

    # CLASS
    if "CLASS" in sections:
        parts.append(f"\nCLASS:\n  " + "\n  ".join(sections["CLASS"]))

    # COMPOUND
    if "COMPOUND" in sections:
        parts.append("\nCOMPOUND entries:")
        for line in sections["COMPOUND"]:
            parts.append(f"  {line}")

    # REACTION
    if "REACTION" in sections:
        parts.append("\nREACTION entries:")
        for line in sections["REACTION"]:
            parts.append(f"  {line}")

    # MODULE
    if "MODULE" in sections:
        parts.append("\nMODULE references:")
        for line in sections["MODULE"]:
            parts.append(f"  {line}")

    # ORTHOLOGY
    if "ORTHOLOGY" in sections:
        parts.append("\nORTHOLOGY:")
        for line in sections["ORTHOLOGY"]:
            parts.append(f"  {line}")

    # PATHWAY (for modules)
    if "PATHWAY" in sections:
        parts.append("\nPATHWAY:")
        for line in sections["PATHWAY"]:
            parts.append(f"  {line}")

    # REFERENCE
    if "REFERENCE" in sections:
        parts.append("\nREFERENCE:")
        for line in sections["REFERENCE"][:10]:  # limit references
            parts.append(f"  {line}")

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# 2. query_kegg_compound
# ---------------------------------------------------------------------------

def query_kegg_compound(compound_id: str) -> str:
    """Query KEGG for a compound entry.

    Extracts NAME, FORMULA, EXACT_MASS, REACTION, PATHWAY, ENZYME, and
    DBLINKS from the result.

    Returns formatted text or an error message string.
    """
    raw = _kegg_get(compound_id)
    if raw.startswith("[ERROR]"):
        return raw

    sections = _parse_kegg_sections(raw)
    parts: list[str] = [f"=== KEGG Compound: {compound_id} ==="]

    for key in ("NAME", "FORMULA", "EXACT_MASS", "MOL_WEIGHT"):
        if key in sections:
            parts.append(f"\n{key}:\n  " + "\n  ".join(sections[key]))

    if "REACTION" in sections:
        parts.append("\nREACTION:")
        for line in sections["REACTION"]:
            parts.append(f"  {line}")

    if "PATHWAY" in sections:
        parts.append("\nPATHWAY:")
        for line in sections["PATHWAY"]:
            parts.append(f"  {line}")

    if "ENZYME" in sections:
        parts.append("\nENZYME:")
        for line in sections["ENZYME"]:
            parts.append(f"  {line}")

    if "DBLINKS" in sections:
        parts.append("\nDBLINKS:")
        for line in sections["DBLINKS"]:
            parts.append(f"  {line}")

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# 3. search_kegg
# ---------------------------------------------------------------------------

def search_kegg(database: str, keyword: str) -> str:
    """Search a KEGG database for a keyword.

    Parameters
    ----------
    database : str
        One of 'compound', 'reaction', 'enzyme', 'module', 'pathway'.
    keyword : str
        Free-text search term.

    Returns a formatted list of up to 30 matches, or an error message string.
    """
    valid_dbs = {"compound", "reaction", "enzyme", "module", "pathway"}
    if database not in valid_dbs:
        return (
            f"[ERROR] Invalid database '{database}'. "
            f"Choose from: {', '.join(sorted(valid_dbs))}"
        )

    url = f"{KEGG_BASE}/find/{database}/{quote(keyword, safe='')}"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
    except requests.RequestException as exc:
        return f"[ERROR] Network error searching KEGG {database} for '{keyword}': {exc}"

    text = resp.text.strip()
    if not text:
        return f"No results found in KEGG {database} for '{keyword}'."

    lines = text.splitlines()
    total = len(lines)
    lines = lines[:30]

    parts: list[str] = [
        f"=== KEGG search: {database} / '{keyword}' "
        f"({total} total result{'s' if total != 1 else ''}, showing {len(lines)}) ==="
    ]
    for line in lines:
        cols = line.split("\t", 1)
        entry_id = cols[0].strip()
        description = cols[1].strip() if len(cols) > 1 else ""
        parts.append(f"  {entry_id:20s} {description}")

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# 4. search_pubmed
# ---------------------------------------------------------------------------

def search_pubmed(query: str, max_results: int = 10) -> str:
    """Search PubMed for articles matching *query* and return abstracts.

    Uses NCBI E-utilities (esearch + efetch).  Returns formatted text with
    title, authors, journal, year, and abstract for each paper, or an error
    message string on failure.
    """
    # Step 1 -- search for PMIDs
    esearch_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    )
    params_search = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "sort": "relevance",
        "retmode": "json",
    }
    try:
        resp_search = requests.get(esearch_url, params=params_search, timeout=30)
        resp_search.raise_for_status()
    except requests.RequestException as exc:
        return f"[ERROR] PubMed search failed: {exc}"

    try:
        data = resp_search.json()
        id_list = data.get("esearchresult", {}).get("idlist", [])
    except (ValueError, KeyError) as exc:
        return f"[ERROR] Failed to parse PubMed search response: {exc}"

    if not id_list:
        return f"No PubMed results found for query: '{query}'"

    # Step 2 -- fetch abstracts
    pmid_str = ",".join(id_list)
    efetch_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    )
    params_fetch = {
        "db": "pubmed",
        "id": pmid_str,
        "rettype": "abstract",
        "retmode": "text",
    }
    try:
        resp_fetch = requests.get(efetch_url, params=params_fetch, timeout=30)
        resp_fetch.raise_for_status()
    except requests.RequestException as exc:
        return (
            f"[ERROR] PubMed efetch failed for PMIDs {pmid_str}: {exc}"
        )

    abstract_text = resp_fetch.text.strip()
    if not abstract_text:
        return f"PubMed returned empty abstracts for PMIDs: {pmid_str}"

    header = (
        f"=== PubMed search: '{query}' "
        f"({len(id_list)} result{'s' if len(id_list) != 1 else ''}) ===\n"
    )
    return header + abstract_text


# ---------------------------------------------------------------------------
# 5. query_kegg_reaction
# ---------------------------------------------------------------------------

def query_kegg_reaction(reaction_id: str) -> str:
    """Query KEGG for a reaction entry.

    Extracts NAME, DEFINITION, EQUATION, ENZYME, and PATHWAY.

    Returns formatted text or an error message string.
    """
    raw = _kegg_get(reaction_id)
    if raw.startswith("[ERROR]"):
        return raw

    sections = _parse_kegg_sections(raw)
    parts: list[str] = [f"=== KEGG Reaction: {reaction_id} ==="]

    for key in ("NAME", "DEFINITION", "EQUATION"):
        if key in sections:
            parts.append(f"\n{key}:\n  " + "\n  ".join(sections[key]))

    if "ENZYME" in sections:
        parts.append("\nENZYME:")
        for line in sections["ENZYME"]:
            parts.append(f"  {line}")

    if "PATHWAY" in sections:
        parts.append("\nPATHWAY:")
        for line in sections["PATHWAY"]:
            parts.append(f"  {line}")

    if "ORTHOLOGY" in sections:
        parts.append("\nORTHOLOGY:")
        for line in sections["ORTHOLOGY"]:
            parts.append(f"  {line}")

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# 6. query_kegg_orthology
# ---------------------------------------------------------------------------

def query_kegg_orthology(ko_id: str) -> str:
    """Query KEGG for an orthology (KO) entry.

    Returns entry name, definition, associated pathways/modules, and a
    summary of organism genes.  Useful for finding cross-species orthologs
    after Phase 1 identifies a KO ID from a reaction.
    """
    raw = _kegg_get(ko_id)
    if raw.startswith("[ERROR]"):
        return raw

    sections = _parse_kegg_sections(raw)
    parts: list[str] = [f"=== KEGG Orthology: {ko_id} ==="]

    for key in ("NAME", "DEFINITION"):
        if key in sections:
            parts.append(f"\n{key}:\n  " + "\n  ".join(sections[key]))

    for key in ("PATHWAY", "MODULE", "DBLINKS"):
        if key in sections:
            parts.append(f"\n{key}:")
            for line in sections[key][:15]:
                parts.append(f"  {line}")

    if "GENES" in sections:
        gene_count = len(sections["GENES"])
        parts.append(f"\nGENES: {gene_count} organism entries (first 20):")
        for line in sections["GENES"][:20]:
            parts.append(f"  {line}")

    return "\n".join(parts)


# ===========================================================================
# UniProt REST API  (Phase 3 — Enzyme Sourcing)
# ===========================================================================

UNIPROT_BASE = "https://rest.uniprot.org"


def uniprot_search(query: str, max_results: int = 15) -> str:
    """Search UniProt for proteins matching *query*.

    The query uses UniProt query syntax, for example:
      - ``ec:1.1.1.27 AND reviewed:true``
      - ``gene:lgtA AND organism_name:neisseria AND reviewed:true``
      - ``name:"formate nitrite transporter" AND organism_id:5833``

    Returns a TSV-formatted table with accession, protein name, organism,
    gene names, length, and EC number.
    """
    fields = "accession,protein_name,organism_name,gene_names,length,ec"
    url = (
        f"{UNIPROT_BASE}/uniprotkb/search"
        f"?query={quote(query, safe='')}"
        f"&format=tsv&fields={fields}&size={max_results}"
    )
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 400:
            return f"[ERROR] UniProt query syntax error for: '{query}'"
        resp.raise_for_status()
    except requests.RequestException as exc:
        return f"[ERROR] UniProt search failed: {exc}"

    text = resp.text.strip()
    if not text or text.count("\n") == 0:
        return f"No UniProt results for query: '{query}'"

    lines = text.splitlines()
    header = lines[0]
    rows = lines[1:]
    if not rows:
        return f"No UniProt results for query: '{query}'"

    parts = [
        f"=== UniProt search: '{query}' ({len(rows)} results) ===",
        header,
    ]
    parts.extend(rows[:max_results])
    return "\n".join(parts)


def uniprot_entry(accession: str) -> str:
    """Fetch detailed information for a single UniProt entry.

    Returns function, catalytic activity, kinetic parameters, cofactors,
    subcellular location, cross-references (PDB, Pfam, InterPro), and key
    publications.
    """
    url = f"{UNIPROT_BASE}/uniprotkb/{quote(accession, safe='')}?format=json"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 404:
            return f"[ERROR] UniProt entry '{accession}' not found."
        resp.raise_for_status()
    except requests.RequestException as exc:
        return f"[ERROR] UniProt fetch failed for '{accession}': {exc}"

    try:
        data = resp.json()
    except ValueError:
        return f"[ERROR] Failed to parse UniProt JSON for '{accession}'."

    parts: list[str] = [f"=== UniProt Entry: {accession} ==="]

    # Protein name
    prot = data.get("proteinDescription", {})
    rec_name = prot.get("recommendedName", {}).get("fullName", {}).get("value", "")
    if rec_name:
        parts.append(f"\nProtein: {rec_name}")

    # Gene names
    genes = data.get("genes", [])
    if genes:
        names = []
        for g in genes:
            for n in g.get("geneName", [{}]):
                if isinstance(n, dict):
                    names.append(n.get("value", ""))
        if names:
            parts.append(f"Gene(s): {', '.join(names)}")

    # Organism
    org = data.get("organism", {}).get("scientificName", "")
    if org:
        parts.append(f"Organism: {org}")

    # Sequence length
    seq_info = data.get("sequence", {})
    if seq_info.get("length"):
        parts.append(f"Length: {seq_info['length']} aa")

    # Function comments
    for comment in data.get("comments", []):
        ctype = comment.get("commentType", "")
        texts = comment.get("texts", [])
        if ctype == "FUNCTION" and texts:
            parts.append(f"\nFunction: {texts[0].get('value', '')}")
        if ctype == "CATALYTIC ACTIVITY":
            rxn = comment.get("reaction", {})
            if rxn.get("name"):
                parts.append(f"\nCatalytic activity: {rxn['name']}")
            ec_list = [e.get("value", "") for e in rxn.get("ecNumber", []) if isinstance(e, dict)]
            if not ec_list:
                ec_val = rxn.get("ecNumber", "")
                if ec_val and isinstance(ec_val, str):
                    ec_list = [ec_val]
            if ec_list:
                parts.append(f"EC: {', '.join(ec_list)}")
        if ctype == "BIOPHYSICOCHEMICAL PROPERTIES":
            kp = comment.get("kineticParameters", {})
            for km in kp.get("michaelisConstants", []):
                parts.append(
                    f"Km: {km.get('constant', '')} {km.get('unit', '')} "
                    f"for {km.get('substrate', '')}"
                )
            for vmax in kp.get("maximumVelocities", []):
                parts.append(f"Vmax: {vmax.get('velocity', '')} {vmax.get('unit', '')}")
        if ctype == "COFACTOR":
            for cof in comment.get("cofactors", []):
                parts.append(f"Cofactor: {cof.get('name', '')}")
        if ctype == "SUBCELLULAR LOCATION":
            locs = []
            for sl in comment.get("subcellularLocations", []):
                loc = sl.get("location", {}).get("value", "")
                if loc:
                    locs.append(loc)
            if locs:
                parts.append(f"\nSubcellular location: {', '.join(locs)}")
        if ctype == "SUBUNIT" and texts:
            parts.append(f"Subunit: {texts[0].get('value', '')[:200]}")

    # Cross-references (PDB, Pfam, InterPro)
    xrefs = data.get("uniProtKBCrossReferences", [])
    pdb_ids = [x["id"] for x in xrefs if x.get("database") == "PDB"]
    pfam_ids = [x["id"] for x in xrefs if x.get("database") == "Pfam"]
    ipr_ids = [x["id"] for x in xrefs if x.get("database") == "InterPro"]
    if pdb_ids:
        parts.append(f"\nPDB: {', '.join(pdb_ids[:10])}")
    if pfam_ids:
        parts.append(f"Pfam: {', '.join(pfam_ids)}")
    if ipr_ids:
        parts.append(f"InterPro: {', '.join(ipr_ids[:10])}")

    # Key references (first 3)
    refs = data.get("references", [])
    if refs:
        parts.append("\nKey references:")
        for ref in refs[:3]:
            cit = ref.get("citation", {})
            title = cit.get("title", "")
            journal = cit.get("journal", "")
            year = cit.get("publicationDate", "")
            if title:
                parts.append(f"  - {title} ({journal}, {year})")

    return "\n".join(parts)


# ===========================================================================
# InterPro / Pfam  (Phase 3 — Enzyme Sourcing)
# ===========================================================================

INTERPRO_BASE = "https://www.ebi.ac.uk/interpro/api"


def interpro_entry(entry_id: str) -> str:
    """Query InterPro or Pfam for a protein family / domain entry.

    Accepts InterPro IDs (e.g. 'IPR000292') or Pfam IDs (e.g. 'PF01226').
    Returns the entry name, description, type, member database counts, and
    species distribution summary.
    """
    if entry_id.startswith("PF"):
        db = "pfam"
    elif entry_id.startswith("IPR"):
        db = "interpro"
    else:
        db = "interpro"

    url = f"{INTERPRO_BASE}/entry/{db}/{quote(entry_id, safe='')}?format=json"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 404:
            return f"[ERROR] InterPro entry '{entry_id}' not found."
        resp.raise_for_status()
    except requests.RequestException as exc:
        return f"[ERROR] InterPro query failed for '{entry_id}': {exc}"

    try:
        data = resp.json()
    except ValueError:
        return f"[ERROR] Failed to parse InterPro JSON for '{entry_id}'."

    meta = data.get("metadata", data)
    parts: list[str] = [f"=== InterPro Entry: {entry_id} ==="]

    name = meta.get("name", "")
    if isinstance(name, dict):
        name = name.get("name", name.get("short", "N/A"))
    parts.append(f"Name: {name}")
    parts.append(f"Type: {meta.get('type', 'N/A')}")
    parts.append(f"Source database: {meta.get('source_database', 'N/A')}")

    desc = meta.get("description", [])
    if isinstance(desc, list):
        for block in desc[:2]:
            if isinstance(block, dict):
                for item in block.get("p", block.get("text", [])):
                    text = item if isinstance(item, str) else (item.get("text", "") if isinstance(item, dict) else "")
                    if text:
                        parts.append(f"\nDescription: {text[:500]}")

    counters = meta.get("counters", {})
    if counters:
        parts.append(f"\nProteins: {counters.get('proteins', 'N/A')}")
        parts.append(f"Proteomes: {counters.get('proteomes', 'N/A')}")
        parts.append(f"Structures: {counters.get('structures', 'N/A')}")

    return "\n".join(parts)
