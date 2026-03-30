"""
External database query tools for pathway discovery.

Provides functions to query KEGG and PubMed for pathway, compound, reaction,
and literature information. All functions return string results and handle
errors gracefully (never raising exceptions to the caller).
"""

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
