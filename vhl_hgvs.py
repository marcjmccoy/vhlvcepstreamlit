import re

def parse_vhl_hgvs(hgvs_full: str):
    """
    Shared VHL‑specific HGVS parser.

    Accepts:
      - 'NM_000551.4(VHL):c.160A>T (p.Met54Leu)'
      - 'NM_000551.4:c.160A>T (p.Met54Leu)'

    Returns:
        transcript (str or None)
        cdna (str or None)
        protein (str or None, 'p.'‑prefixed)
    """
    if not isinstance(hgvs_full, str):
        return None, None, None

    # Transcript, allowing optional gene name in parentheses.
    m_tx = re.match(r"^([A-Z0-9_.]+)\(", hgvs_full)
    if not m_tx:
        m_tx = re.match(r"^([A-Z0-9_.]+):", hgvs_full)
    transcript = m_tx.group(1) if m_tx else "NM_000551.4"

    # cDNA HGVS.
    m_c = re.search(r"(c\.[^ \)]+)", hgvs_full)
    cdna = m_c.group(1) if m_c else None

    # Protein HGVS.
    m_p = re.search(r"\(p\.([^)\s]+)\)", hgvs_full)
    protein = f"p.{m_p.group(1)}" if m_p else None

    return transcript, cdna, protein
