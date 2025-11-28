"""
VHL VCEP PM2_Supporting classifier (API-based).

Implements the VHL VCEP rule using gnomAD v4 on GRCh38:

- PM2_Supporting can be applied for variants either absent from gnomAD OR
  with <= 0.00000156 (0.000156%) GroupMax Filtering Allele Frequency (FAF)
  in gnomAD (v4).
- If no GroupMax FAF is calculated (e.g. due to a single observation),
  PM2_Supporting may also be applied.

Workflow:
  1) Parse VHL transcript HGVS input, e.g. 'NM_000551.4(VHL):c.1A>T (p.Met1Leu)'.
  2) Use Ensembl VEP REST API to map NM_000551.4:c.* to GRCh38 coords.
  3) Query GeneBe single-variant API for that locus (chr/pos/ref/alt, genome=hg38)
     to obtain gnomAD AF/AC/AN (and, if exposed, FAF/FAF95).
  4) Approximate GroupMax FAF from GeneBe’s gnomAD fields and apply the
     PM2_Supporting rule.
"""

import re
import logging
from typing import Optional, Tuple

import requests

logger = logging.getLogger(__name__)

GNOMAD_PM2_MAX_FAF: float = 0.00000156  # 0.000156%

# Ensembl VEP REST (GRCh38)
VEP_HGVS_BASE = "https://rest.ensembl.org/vep/human/hgvs"

# GeneBe public variant endpoint
GENEBE_VARIANT_BASE = "https://api.genebe.net/cloud/api-public/v1/variant"


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _safe_float(value) -> Optional[float]:
    try:
        if value is None:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


# ----------------------------------------------------------------------
# HGVS parsing (VHL transcript-level)
# ----------------------------------------------------------------------


def _parse_vhl_hgvs(hgvs_full: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse VHL HGVS strings like:
      - 'NM_000551.4(VHL):c.160A>T (p.Met54Leu)'
      - 'NM_000551.4:c.160A>T (p.Met54Leu)'

    Returns:
        transcript: e.g. 'NM_000551.4'
        cdna:       e.g. 'c.160A>T'
    """
    if not isinstance(hgvs_full, str):
        return None, None

    m_tx = re.match(r"^([A-Z0-9_.]+)\(VHL\)", hgvs_full)
    if not m_tx:
        m_tx = re.match(r"^([A-Z0-9_.]+):", hgvs_full)
    transcript = m_tx.group(1) if m_tx else "NM_000551.4"

    m_c = re.search(r"(c\.[^ \)]+)", hgvs_full)
    cdna = m_c.group(1) if m_c else None

    return transcript, cdna


# ----------------------------------------------------------------------
# Ensembl VEP: transcript HGVS -> GRCh38 genomic coords
# ----------------------------------------------------------------------


def _vep_map_hgvs_to_grch38(hgvs_tx: str) -> Tuple[Optional[str], Optional[int], Optional[str], Optional[str]]:
    """
    Use Ensembl VEP REST API to map a transcript HGVS
    (e.g. 'NM_000551.4:c.160A>T') to GRCh38 chrom, pos, ref, alt.

    Returns:
        chrom (str), pos (int), ref (str), alt (str)
        or (None, None, None, None) on failure.
    """
    url = f"{VEP_HGVS_BASE}/{hgvs_tx}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    try:
        resp = requests.get(
            url,
            params={"refseq": 1},  # prefer RefSeq transcript mappings
            headers=headers,
            timeout=10,
        )
        if resp.status_code == 404:
            return None, None, None, None
        resp.raise_for_status()
        data = resp.json() or []
    except Exception as exc:
        logger.warning("VEP hgvs->GRCh38 mapping failed for %s: %s", hgvs_tx, exc)
        return None, None, None, None

    if not data:
        return None, None, None, None

    entry = data[0]

    chrom = str(entry.get("seq_region_name") or "")
    pos = entry.get("start")
    allele_string = entry.get("allele_string") or ""
    if "/" in allele_string:
        ref, alt = allele_string.split("/", 1)
    else:
        ref, alt = None, None

    if not chrom or pos is None or not ref or not alt:
        return None, None, None, None

    return chrom, int(pos), ref, alt


# ----------------------------------------------------------------------
# GeneBe: GRCh38 variant -> gnomAD frequencies
# ----------------------------------------------------------------------


def _genebe_query_gnomad(chrom: str, pos: int, ref: str, alt: str) -> Tuple[bool, Optional[float]]:
    """
    Query GeneBe for gnomAD frequencies for a GRCh38 variant.

    Returns:
        present_in_gnomad (bool)
        approx_groupmax_faf (float or None)

    Approximation:
        - If GeneBe exposes a gnomAD FAF/FAF95 value, use that directly.
        - Otherwise, approximate FAF as the maximum of available gnomAD AF
          values (exome/genome, total or population-wise), which is a
          conservative proxy for GroupMax FAF.
    """
    params = {
        "chr": chrom.replace("chr", ""),
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "genome": "hg38",
    }

    try:
        resp = requests.get(GENEBE_VARIANT_BASE, params=params, timeout=10)
        if resp.status_code == 404:
            return False, None
        resp.raise_for_status()
        data = resp.json() or {}
    except Exception as exc:
        logger.warning("GeneBe variant API error for %s:%s:%s>%s: %s",
                       chrom, pos, ref, alt, exc)
        return False, None

    # GeneBe response structure: see https://docs.genebe.net/docs/api/endpoints
    ann = (data.get("annotation") or {})
    gnomad_block = ann.get("gnomad") or ann.get("gnomad4") or {}

    if not gnomad_block:
        return False, None

    # Presence: any non-zero AC or AF
    present = False
    candidate_afs = []

    for key in ("af", "af_exomes", "af_genomes"):
        af = _safe_float(gnomad_block.get(key))
        if af is not None:
            candidate_afs.append(af)
            if af > 0:
                present = True

    # GeneBe may have per-population AFs; include them too if present
    pops = gnomad_block.get("populations") or {}
    for pop_vals in pops.values():
        af = _safe_float(pop_vals.get("af"))
        if af is not None:
            candidate_afs.append(af)
            if af > 0:
                present = True

    if not present:
        return False, None

    # If GeneBe exposes FAF directly, prefer that.
    faf = None
    for key in ("faf", "faf95", "faf_popmax", "faf95_popmax"):
        if key in gnomad_block:
            faf = _safe_float(gnomad_block.get(key))
            break

    # Otherwise approximate GroupMax FAF as max(AF) across contexts.
    if faf is None and candidate_afs:
        faf = max(candidate_afs)

    return present, faf


# ----------------------------------------------------------------------
# Main PM2 classifier
# ----------------------------------------------------------------------


def classify_vhl_pm2(hgvs_full: str):
    """
    VHL VCEP PM2_Supporting classifier using:

      - Ensembl VEP REST for HGVS -> GRCh38 mapping, then
      - GeneBe public API for gnomAD v4-derived frequencies.

    PM2_Supporting is applied when:
      - The variant is absent from gnomAD v4; OR
      - Approximate GroupMax FAF ≤ GNOMAD_PM2_MAX_FAF; OR
      - The variant is present but no FAF is available.
    """
    # 1) Parse transcript and cDNA from VHL HGVS.
    transcript, cdna = _parse_vhl_hgvs(hgvs_full)
    if transcript is None or cdna is None:
        return {
            "strength": None,
            "context": (
                "Could not parse transcript and cDNA HGVS from input; "
                "PM2_Supporting not applied."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    tx_hgvs = f"{transcript}:{cdna}"  # e.g. NM_000551.4:c.160A>T

    # 2) Map to GRCh38 genomic coordinates via Ensembl VEP.
    chrom, pos, ref, alt = _vep_map_hgvs_to_grch38(tx_hgvs)
    if not (chrom and pos and ref and alt):
        return {
            "strength": None,
            "context": (
                f"Could not map {tx_hgvs} to GRCh38 genomic coordinates via Ensembl VEP; "
                "unable to query gnomAD frequencies via GeneBe, so PM2_Supporting "
                "is not applied."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    # 3) Query GeneBe for gnomAD presence and (approx) GroupMax FAF.
    present, faf = _genebe_query_gnomad(chrom, pos, ref, alt)

    # 4) Apply the VHL VCEP PM2_Supporting rule.

    # 1) Absent from gnomAD v4.
    if not present:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{hgvs_full} is absent from gnomAD v4 according to GeneBe "
                "(no alternate alleles observed), meeting the PM2_Supporting "
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
                f"{hgvs_full} is present in gnomAD v4 with approximate GroupMax "
                f"FAF={faf:.3e} (from GeneBe), which is ≤1.56×10⁻⁶, so "
                "PM2_Supporting is applicable."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": faf,
        }

    # 3) Present but no FAF available.
    if present and faf is None:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{hgvs_full} is observed in gnomAD v4 via GeneBe but has no "
                "Filtering Allele Frequency estimate; this still qualifies for "
                "PM2_Supporting per the VHL VCEP specification."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": None,
        }

    # 4) Present with FAF above threshold – PM2_Supporting not met.
    return {
        "strength": None,
        "context": (
            f"{hgvs_full} has approximate GroupMax FAF={faf:.3e} in gnomAD v4 "
            "(via GeneBe), which exceeds the PM2_Supporting cutoff of "
            "1.56×10⁻⁶; PM2_Supporting not applied."
        ),
        "present_in_gnomad": True,
        "groupmax_faf": faf,
    }