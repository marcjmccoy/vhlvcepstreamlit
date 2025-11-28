import re

def parse_vhl_hgvs(hgvs_full: str):
    if not isinstance(hgvs_full, str):
        return None, None, None

    m_tx = re.match(r"^([A-Z0-9_.]+)\(", hgvs_full)
    if not m_tx:
        m_tx = re.match(r"^([A-Z0-9_.]+):", hgvs_full)
    transcript = m_tx.group(1) if m_tx else "NM_000551.4"

    m_c = re.search(r"(c\.[^ \)]+)", hgvs_full)
    cdna = m_c.group(1) if m_c else None

    m_p = re.search(r"\(p\.([^)\s]+)\)", hgvs_full)
    protein = f"p.{m_p.group(1)}" if m_p else None

    return transcript, cdna, protein