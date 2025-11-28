# vhl_pm2.py

"""
VHL VCEP PM2_Supporting classifier using gnomAD v4 via MyVariant.info.

Implements the VHL VCEP specification:

- PM2_Supporting can be applied for variants either absent from gnomAD or with
  <= 0.00000156 (0.000156%) GroupMax Filtering Allele Frequency in gnomAD
  (based on gnomAD v4 release).
- If no GroupMax Filtering Allele Frequency is calculated (e.g. single
  observation), PM2_Supporting may also be applied.

This module uses MyVariant.info to access gnomAD v4 fields and applies PM2 only
under those conditions.
"""

import logging
import requests

from vhl_hgvs import parse_vhl_hgvs, format_vhl_hgvs


logger = logging.getLogger(__name__)


# VHL VCEP PM2_Supporting GroupMax FAF threshold (gnomAD v4)
GNOMAD_PM2_MAX_FAF: float = 0.00000156  # 0.000156%

# MyVariant.info base URL
MYVARIANT_BASE = "https://myvariant.info/v1"


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


def _query_myvariant_gnomad_v4(hgvs_tx: str):
    """
    Query MyVariant.info for gnomAD v4 data for a transcript HGVS.

    Args:
        hgvs_tx:
            Transcript HGVS string, e.g. 'NM_000551.4:c.1A>T'.

    Returns:
        present_in_gnomad (bool):
            True if any non-zero AC or AF is reported in gnomAD v4.
        groupmax_faf (float or None):
            GroupMax Filtering Allele Frequency (FAF) if available, else None.
    """
    try:
        resp = requests.get(
            f"{MYVARIANT_BASE}/variant/{hgvs_tx}",
            params={"fields": "gnomad_v4"},
            timeout=10,
        )
        if resp.status_code == 404:
            return False, None
        resp.raise_for_status()
        data = resp.json() or {}
    except Exception as exc:
        logger.warning("MyVariant lookup failed for %s: %s", hgvs_tx, exc)
        return False, None

    g4 = data.get("gnomad_v4") or {}

    # Prefer GroupMax FAF; if not present, try alternate names or None.
    faf = (
        _safe_float(g4.get("groupmax_faf"))
        or _safe_float(g4.get("group_max_faf"))
        or _safe_float(g4.get("faf"))
    )

    # Presence: any non-zero AC or AF implies present in gnomAD v4.
    ac = _safe_int(g4.get("ac"), default=0)
    af = _safe_float(g4.get("af"))

    present = (ac > 0) or (af is not None and af > 0.0)
    return present, faf


def classify_vhl_pm2(hgvs_full: str):
    """
    VHL VCEP PM2_Supporting classifier using gnomAD v4 via MyVariant.info.

    PM2_Supporting is applied when:
      - The variant is absent from gnomAD v4; OR
      - The variant has GroupMax FAF <= GNOMAD_PM2_MAX_FAF; OR
      - The variant is present but no GroupMax FAF is calculated.
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

    display_hgvs = format_vhl_hgvs(transcript, cdna)   # e.g. NM_000551.4(VHL):c.1A>T
    myvariant_hgvs = f"{transcript}:{cdna}"            # e.g. NM_000551.4:c.1A>T

    present, faf = _query_myvariant_gnomad_v4(myvariant_hgvs)

    # 1) Absent from gnomAD v4.
    if not present:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{display_hgvs} is absent from gnomAD v4 (no alternate alleles "
                "observed in MyVariant/gnomAD v4), meeting the PM2_Supporting "
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
                f"{display_hgvs} is present in gnomAD v4 with GroupMax FAF={faf:.3e}, "
                "which is ≤1.56×10⁻⁶, so PM2_Supporting is applicable."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": faf,
        }

    # 3) Present but no GroupMax FAF computed.
    if present and faf is None:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{display_hgvs} is observed in gnomAD v4 via MyVariant.info but lacks "
                "a GroupMax Filtering Allele Frequency estimate (e.g. single "
                "observation); this still qualifies for PM2_Supporting."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": None,
        }

    # 4) Present with FAF above threshold – PM2_Supporting not met.
    return {
        "strength": None,
        "context": (
            f"{display_hgvs} has GroupMax FAF={faf:.3e} in gnomAD v4 via "
            "MyVariant.info, which exceeds the PM2_Supporting cutoff of "
            "1.56×10⁻⁶; PM2_Supporting not applied."
        ),
        "present_in_gnomad": True,
        "groupmax_faf": faf,
    }