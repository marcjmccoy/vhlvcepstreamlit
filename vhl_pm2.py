import requests


# VHL VCEP PM2_Supporting GroupMax FAF threshold (gnomAD v4).[web:7][web:69]
GNOMAD_PM2_MAX_FAF: float = 0.00000156  # 0.000156%


# ClinGen LDH base URL (front end to gnomAD v4).[web:70]
CLINGEN_LDH_BASE: str = "https://ldh.clinicalgenome.org/ldh/"


def _safe_float(value):
    """
    Convert value to float, returning None if conversion fails.
    """
    try:
        if value is None:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _safe_int(value, default=0):
    """
    Convert value to int, returning `default` if conversion fails.
    """
    try:
        if value is None:
            return default
        return int(value)
    except (TypeError, ValueError):
        return default


def _query_clingen_gnomad_v4(cdna_hgvs: str, transcript: str = "NM_000551.4"):
    """
    Query ClinGen LDH for gnomAD v4 presence and GroupMax FAF for a VHL variant.

    Args:
        cdna_hgvs:
            cDNA HGVS (e.g. 'c.451A>G').
        transcript:
            RefSeq transcript (e.g. 'NM_000551.4').

    Returns:
        present_in_gnomad (bool):
            True if allele count (AC) > 0 in gnomAD v4.
        groupmax_faf (float or None):
            GroupMax Filtering Allele Frequency (FAF) if available, otherwise None.
    """
    if not cdna_hgvs:
        return False, None

    full_hgvs = f"{transcript}:{cdna_hgvs}"

    params = {
        "q": full_hgvs,
        "dataset": "gnomAD_v4",  # dataset key used in LDH for gnomAD v4.[web:70]
        "rows": 1,
    }

    try:
        resp = requests.get(
            CLINGEN_LDH_BASE + "search",
            params=params,
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json() or {}
    except Exception:
        # On any API failure, treat as no reliable population data.
        return False, None

    docs = data.get("docs") or []
    if not docs:
        return False, None

    doc = docs[0]

    # gnomAD/LDH can expose GroupMax FAF with slightly different field names.[web:69][web:74]
    faf_raw = (
        doc.get("jointGrpMaxFAF95")
        or doc.get("grpmaxFAF")
        or doc.get("groupmax_faf")
        or doc.get("groupMaxFaf")
    )
    faf = _safe_float(faf_raw)

    # Allele count can also be named differently.
    ac_raw = doc.get("AC") or doc.get("allele_count") or doc.get("ac")
    ac = _safe_int(ac_raw, default=0)

    present = ac > 0
    return present, faf


def classify_vhl_pm2(hgvs_full: str):
    """
    VHL VCEP PM2_Supporting classifier using gnomAD v4 via ClinGen LDH.

    Criteria (Supporting strength) per VHL VCEP GN078:[web:7]
      - Variant is absent from gnomAD; OR
      - Variant has GroupMax Filtering Allele Frequency (gnomAD v4)
        ≤ 0.00000156 (0.000156%); OR
      - Variant is present but no GroupMax FAF is calculated
        (e.g., single observation).

    Args:
        hgvs_full:
            Full VHL HGVS including transcript, e.g. 'NM_000551.4:c.451A>G'.

    Returns:
        dict with keys:
            - strength: "PM2_Supporting" or None
            - context: human‑readable explanation
            - present_in_gnomad: bool
            - groupmax_faf: float or None
    """
    # _parse_vhl_hgvs should be provided elsewhere in your code base.
    # It must return (transcript, cdna_hgvs, protein_hgvs or None).
    transcript, cdna, _ = _parse_vhl_hgvs(hgvs_full)

    if cdna is None or transcript is None:
        return {
            "strength": None,
            "context": (
                "Could not parse transcript and cDNA HGVS from input; "
                "PM2_Supporting not applied."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    present, faf = _query_clingen_gnomad_v4(cdna, transcript=transcript)

    # 1) Absent from gnomAD v4.
    if not present:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{transcript}:{cdna} is absent from gnomAD v4 (no alternate "
                "alleles observed), meeting the VHL VCEP PM2_Supporting "
                "population criterion."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    # 2) Present with GroupMax FAF at or below threshold.
    if faf is not None and faf <= GNOMAD_PM2_MAX_FAF:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{transcript}:{cdna} is present in gnomAD v4 with "
                f"GroupMax FAF={faf:.3e}, which is ≤1.56×10⁻⁶, so "
                "PM2_Supporting is applicable per VHL VCEP."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": faf,
        }

    # 3) Present but no GroupMax FAF computed (e.g., single observation).
    if present and faf is None:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{transcript}:{cdna} is observed in gnomAD v4 but lacks a "
                "GroupMax FAF estimate (likely a single‑allele site); under "
                "VHL VCEP this still qualifies for PM2_Supporting."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": None,
        }

    # 4) FAF above the allowed threshold – PM2_Supporting not met.
    return {
        "strength": None,
        "context": (
            f"{transcript}:{cdna} has GroupMax FAF={faf:.3e} in gnomAD v4, "
            "which exceeds the VHL VCEP PM2_Supporting cutoff of 1.56×10⁻⁶; "
            "PM2_Supporting not applied."
        ),
        "present_in_gnomad": True,
        "groupmax_faf": faf,
    }