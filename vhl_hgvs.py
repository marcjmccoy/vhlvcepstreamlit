# vhl_hgvs.py

import re

def parse_vhl_hgvs(hgvs_full: str):
    """
    Shared VHL‑specific HGVS parser.

    Accepts user input like:
      - 'NM_000551.4(VHL):c.160A>T (p.Met54Leu)'
      - 'NM_000551.4:c.160A>T (p.Met54Leu)'

    Returns:
        transcript: e.g. 'NM_000551.4'
        cdna:       e.g. 'c.160A>T'
        protein:    e.g. 'p.Met54Leu' or None
    """
    if not isinstance(hgvs_full, str):
        return None, None, None

    # Transcript (strip optional '(VHL)').
    m_tx = re.match(r"^([A-Z0-9_.]+)\(VHL\)", hgvs_full)
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