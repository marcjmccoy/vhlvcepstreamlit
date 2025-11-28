"""
VHL VCEP PM2_Supporting classifier.

Implements the VHL VCEP rule using gnomAD v4 on GRCh38:

- PM2_Supporting can be applied for variants either absent from gnomAD or
  with <= 0.00000156 (0.000156%) GroupMax Filtering Allele Frequency (FAF)
  in gnomAD v4.
- If no GroupMax FAF is calculated (e.g. single observation), PM2_Supporting
  may also be applied.

Workflow:
  1) Parse VHL transcript HGVS input, e.g. 'NM_000551.4(VHL):c.1A>T (p.Met1Leu)'.
  2) Use Ensembl VEP REST API to map NM_000551.4:c.* to GRCh38 genomic coords.
  3) Look up that locus in a local VHL-only gnomAD v4 VCF via LocalGnomadVHL.
  4) Apply the PM2_Supporting rule.
"""

import re
import logging
import requests

from vhl_gnomad_local import LocalGnomadVHL

logger = logging.getLogger(__name__)

# VHL VCEP PM2_Supporting GroupMax FAF threshold (gnomAD v4)
GNOMAD_PM2_MAX_FAF: float = 0.00000156  # 0.000156%

# Ensembl VEP REST (GRCh38)
VEP_HGVS_BASE = "https://rest.ensembl.org/vep/human/hgvs"

# Path to your local VHL-only gnomAD v4 VCF (exomes or merged)
LOCAL_GNOMAD_VHL_VCF = "gnomad.exomes.v4.1.VHL.vcf"

# Single shared instance
_local_gnomad = LocalGnomadVHL(LOCAL_GNOMAD_VHL_VCF)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _safe_float(value):
    try:
        if value is None:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


# ----------------------------------------------------------------------
# HGVS parsing (VHL transcript-level)
# ----------------------------------------------------------------------


def _parse_vhl_hgvs(hgvs_full: str):
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


def _vep_map_hgvs_to_grch38(hgvs_tx: str):
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
# Main PM2 classifier
# ----------------------------------------------------------------------


def classify_vhl_pm2(hgvs_full: str):
    """
    VHL VCEP PM2_Supporting classifier using:

      - Ensembl VEP REST for HGVS -> GRCh38 mapping, then
      - Local gnomAD v4 VHL VCF (via LocalGnomadVHL) for AC and GroupMax FAF.

    PM2_Supporting is applied when:
      - The variant is absent from gnomAD v4; OR
      - GroupMax FAF ≤ GNOMAD_PM2_MAX_FAF; OR
      - The variant is present but no GroupMax FAF is calculated.
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
                "unable to query local gnomAD v4 VHL VCF, so PM2_Supporting is not applied."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    # 3) Query local gnomAD v4 VHL VCF for presence and GroupMax FAF.
    present, faf = _local_gnomad.lookup(chrom, pos, ref, alt)

    # 4) Apply the VHL VCEP PM2_Supporting rule.

    # 1) Absent from gnomAD v4.
    if not present:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{hgvs_full} is absent from local gnomAD v4 VHL VCF "
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
                f"{hgvs_full} is present in local gnomAD v4 VHL VCF "
                f"with GroupMax FAF={faf:.3e}, which is ≤1.56×10⁻⁶, "
                "so PM2_Supporting is applicable."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": faf,
        }

    # 3) Present but no GroupMax FAF computed.
    if present and faf is None:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{hgvs_full} is observed in local gnomAD v4 VHL VCF but lacks a "
                "GroupMax Filtering Allele Frequency estimate (e.g. single "
                "observation); this still qualifies for PM2_Supporting."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": None,
        }

    # 4) Present with FAF above threshold – PM2_Supporting not met.
    return {
        "strength": None,
        "context": (
            f"{hgvs_full} has GroupMax FAF={faf:.3e} in local gnomAD v4 VHL VCF, "
            "which exceeds the PM2_Supporting cutoff of 1.56×10⁻⁶; "
            "PM2_Supporting not applied."
        ),
        "present_in_gnomad": True,
        "groupmax_faf": faf,
    }