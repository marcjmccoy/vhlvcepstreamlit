"""
VHL VCEP PM2_Supporting classifier (gnomAD v4 via ClinGen LDH).

Implements GN078 PM2_Supporting using GroupMax Filtering Allele Frequency (FAF)
from gnomAD v4 accessed through ClinGen LDH, resolving variants via
the ClinGen Allele Registry CA ID.
"""

import logging
import requests

from vhl_hgvs import parse_vhl_hgvs

logger = logging.getLogger(__name__)

GNOMAD_PM2_MAX_FAF: float = 0.00000156  # 0.000156%

CLINGEN_LDH_BASE: str = "https://ldh.clinicalgenome.org/ldh/"
ALLELE_REGISTRY_BASE: str = (
    "https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry"
)


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


def _get_ca_id_from_hgvs(hgvs: str):
    """
    Resolve an HGVS string to a ClinGen CA ID via the Allele Registry.

    Returns:
        e.g. 'CA123456789' or None if not found / error.
    """
    try:
        resp = requests.get(
            f"{ALLELE_REGISTRY_BASE}/alleles",
            params={"hgvs": hgvs},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json() or {}
    except Exception as exc:
        logger.warning("Allele Registry lookup failed for %s: %s", hgvs, exc)
        return None

    # Registry can return dict keyed by CA or list of allele objects.
    if isinstance(data, dict):
        ca_keys = [k for k in data.keys() if isinstance(k, str) and k.startswith("CA")]
        return ca_keys[0] if ca_keys else None

    if isinstance(data, list) and data:
        cand = data[0].get("id") or data[0].get("@id")
        if isinstance(cand, str) and cand.startswith("CA"):
            return cand

    return None


def _query_clingen_gnomad_v4_with_ca(ca_id: str):
    """
    Query ClinGen LDH for gnomAD v4 presence and GroupMax FAF using a CA ID.

    Returns:
        present_in_gnomad (bool)
        groupmax_faf (float or None)
    """
    if not ca_id:
        return False, None

    try:
        resp = requests.get(
            f"{CLINGEN_LDH_BASE}allele/{ca_id}",
            params={"dataset": "gnomAD_v4"},
            timeout=10,
        )
        resp.raise_for_status()
        doc = resp.json() or {}
    except Exception as exc:
        logger.warning("LDH allele lookup failed for %s: %s", ca_id, exc)
        return False, None

    # GroupMax FAF field variants.
    faf_raw = (
        doc.get("jointGrpMaxFAF95")
        or doc.get("grpmaxFAF")
        or doc.get("groupmax_faf")
        or doc.get("groupMaxFaf")
    )
    faf = _safe_float(faf_raw)

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

    full_hgvs = f"{transcript}:{cdna}"

    # Step 1: resolve to CA ID.
    ca_id = _get_ca_id_from_hgvs(full_hgvs)
    if ca_id is None:
        logger.info("No CA ID for %s; treating as no population data.", full_hgvs)
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{full_hgvs} could not be resolved to a ClinGen CA ID; "
                "treating as absent from population datasets and applying "
                "PM2_Supporting conservatively."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    # Step 2: query LDH for gnomAD v4 by CA ID.
    present, faf = _query_clingen_gnomad_v4_with_ca(ca_id)

    # 1) Absent from gnomAD v4.
    if not present:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{full_hgvs} (CA ID {ca_id}) is absent from gnomAD v4 "
                "(no alternate alleles observed), meeting the "
                "PM2_Supporting population criterion."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    # 2) Present with GroupMax FAF at or below threshold.
    if faf is not None and faf <= GNOMAD_PM2_MAX_FAF:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{full_hgvs} (CA ID {ca_id}) is present in gnomAD v4 with "
                f"GroupMax FAF={faf:.3e}, which is ≤1.56×10⁻⁶, so "
                "PM2_Supporting is applicable."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": faf,
        }

    # 3) Present but no GroupMax FAF computed.
    if present and faf is None:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{full_hgvs} (CA ID {ca_id}) is observed in gnomAD v4 but lacks "
                "a GroupMax FAF estimate (likely a single‑allele site); this "
                "still qualifies for PM2_Supporting."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": None,
        }

    # 4) FAF above the threshold – PM2_Supporting not met.
    return {
        "strength": None,
        "context": (
            f"{full_hgvs} (CA ID {ca_id}) has GroupMax FAF={faf:.3e} in gnomAD v4, "
            "which exceeds the PM2_Supporting cutoff of 1.56×10⁻⁶; "
            "PM2_Supporting not applied."
        ),
        "present_in_gnomad": True,
        "groupmax_faf": faf,
    }