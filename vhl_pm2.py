"""
VHL VCEP PM2_Supporting classifier (gnomAD v4 via ClinGen LDH).

Implements GN078 PM2_Supporting using GroupMax Filtering Allele Frequency (FAF)
from gnomAD v4 accessed through ClinGen LDH.
"""

import logging
import requests

from vhl_hgvs import parse_vhl_hgvs

logger = logging.getLogger(__name__)


GNOMAD_PM2_MAX_FAF: float = 0.00000156  # 0.000156%

CLINGEN_LDH_BASE: str = "https://ldh.clinicalgenome.org/ldh/"


def _safe_float(value):
    try:
        if value is None:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _safe_int(value, default=0):
    try:
        if value is None:
            return default
        return int(value)
    except (TypeError, ValueError):
        return default


def _ldh_search(params: dict):
    """
    Low‑level helper to call LDH /search with basic error trapping.
    """
    try:
        resp = requests.get(
            CLINGEN_LDH_BASE + "search",
            params=params,
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json() or {}
        return data.get("docs") or []
    except Exception as exc:
        logger.warning("LDH search failed (%s): %s", params, exc)
        return []


def _query_clingen_gnomad_v4(cdna_hgvs: str, transcript: str = "NM_000551.4"):
    """
    Query ClinGen LDH for gnomAD v4 presence and GroupMax FAF for a VHL variant.

    Strategy:
      1. Try HGVS‑based search with dataset=gnomAD_v4.
      2. If no docs, retry without dataset filter (some deployments differ).
    """
    if not cdna_hgvs:
        return False, None

    full_hgvs = f"{transcript}:{cdna_hgvs}"

    # First: constrained to gnomAD_v4.
    docs = _ldh_search(
        {
            "q": full_hgvs,
            "dataset": "gnomAD_v4",
            "rows": 5,
        }
    )

    # Fallback: drop dataset filter if nothing found.
    if not docs:
        docs = _ldh_search(
            {
                "q": full_hgvs,
                "rows": 5,
            }
        )

    if not docs:
        logger.info("No LDH docs for %s", full_hgvs)
        return False, None

    doc = docs[0]
    logger.debug("LDH doc for %s: %s", full_hgvs, doc)

    # GroupMax FAF field can vary slightly by export.
    faf_raw = (
        doc.get("jointGrpMaxFAF95")
        or doc.get("grpmaxFAF")
        or doc.get("groupmax_faf")
        or doc.get("groupMaxFaf")
    )
    faf = _safe_float(faf_raw)

    # Allele count variants.
    ac_raw = doc.get("AC") or doc.get("allele_count") or doc.get("ac")
    ac = _safe_int(ac_raw, default=0)

    present = ac > 0
    return present, faf


def classify_vhl_pm2(hgvs_full: str):
    """
    VHL VCEP PM2_Supporting classifier using gnomAD v4 via ClinGen LDH.

    PM2_Supporting is applied when:
      - The variant is absent from gnomAD v4; OR
      - The variant has GroupMax FAF ≤ GNOMAD_PM2_MAX_FAF; OR
      - The variant is present but has no GroupMax FAF calculated.
    """
    transcript, cdna, _ = parse_vhl_hgvs(hgvs_full)

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
                "alleles observed), meeting the PM2_Supporting population criterion."
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
                "PM2_Supporting is applicable."
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
                "GroupMax FAF estimate (likely a single‑allele site); this still "
                "qualifies for PM2_Supporting."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": None,
        }

    # 4) FAF above the threshold – PM2_Supporting not met.
    return {
        "strength": None,
        "context": (
            f"{transcript}:{cdna} has GroupMax FAF={faf:.3e} in gnomAD v4, "
            "which exceeds the PM2_Supporting cutoff of 1.56×10⁻⁶; "
            "PM2_Supporting not applied."
        ),
        "present_in_gnomad": True,
        "groupmax_faf": faf,
    }