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
  3) Query gnomAD v4 REST API for that locus (chrom-pos-ref-alt).
  4) Apply the PM2_Supporting rule.
"""

import re
import logging
import requests

logger = logging.getLogger(__name__)

GNOMAD_PM2_MAX_FAF: float = 0.00000156  # 0.000156%

# Ensembl VEP REST (GRCh38)
VEP_HGVS_BASE = "https://rest.ensembl.org/vep/human/hgvs"

# gnomAD REST variant endpoint (v4 / GRCh38)
GNOMAD_VARIANT_BASE = "https://gnomad.broadinstitute.org/api/variant"


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


def _safe_int(value, default=0):
    try:
        if value is None:
            return default
        return int(value)
    except (TypeError, ValueError):
        return default


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
# gnomAD v4 REST query by GRCh38 coordinates
# ----------------------------------------------------------------------


def _gnomad_query_variant_grch38(chrom: str, pos: int, ref: str, alt: str):
    """
    Query gnomAD REST API for a GRCh38 variant and return
    presence and GroupMax FAF (if available) from the v4 dataset.

    Uses the “variantId” pattern 'chrom-pos-ref-alt', e.g. '3-10142007-A-T'.

    Returns:
        present_in_gnomad (bool)
        groupmax_faf (float or None)
    """
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"

    try:
        resp = requests.get(
            GNOMAD_VARIANT_BASE,
            params={"variant": variant_id, "dataset": "gnomad_r4"},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json() or {}
    except Exception as exc:
        logger.warning("gnomAD REST error for %s: %s", variant_id, exc)
        return False, None

    variant = data.get("variant") or {}
    if not variant:
        return False, None

    genome = variant.get("genome") or {}
    ac = _safe_int(genome.get("ac"), default=0)
    present = ac > 0

    faf_block = genome.get("faf95") or {}
    faf_popmax = _safe_float(faf_block.get("popmax"))

    return present, faf_popmax


# ----------------------------------------------------------------------
# Main PM2 classifier
# ----------------------------------------------------------------------


def classify_vhl_pm2(hgvs_full: str):
    """
    VHL VCEP PM2_Supporting classifier using:

      - Ensembl VEP REST for HGVS -> GRCh38 mapping, then
      - gnomAD v4 REST API for AC and GroupMax FAF.

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
                "unable to query gnomAD v4, so PM2_Supporting is not applied."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    # 3) Query gnomAD v4 for presence and GroupMax FAF.
    present, faf = _gnomad_query_variant_grch38(chrom, pos, ref, alt)

    # 4) Apply the VHL VCEP PM2_Supporting rule.

    # 1) Absent from gnomAD v4.
    if not present:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{hgvs_full} is absent from gnomAD v4 "
                f"(no alternate alleles observed at {chrom}-{pos}-{ref}-{alt}), "
                "meeting the PM2_Supporting population criterion."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    # 2) Present with GroupMax FAF at or below threshold.
    if faf is not None and faf <= GNOMAD_PM2_MAX_FAF:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{hgvs_full} is present in gnomAD v4 at {chrom}-{pos}-{ref}-{alt} "
                f"with GroupMax FAF={faf:.3e}, which is ≤1.56×10⁻⁶, so "
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
                f"{hgvs_full} is observed in gnomAD v4 at {chrom}-{pos}-{ref}-{alt} "
                "but lacks a GroupMax Filtering Allele Frequency estimate "
                "(e.g. single observation); this still qualifies for PM2_Supporting."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": None,
        }

    # 4) Present with FAF above threshold – PM2_Supporting not met.
    return {
        "strength": None,
        "context": (
            f"{hgvs_full} has GroupMax FAF={faf:.3e} in gnomAD v4 at "
            f"{chrom}-{pos}-{ref}-{alt}, which exceeds the "
            "PM2_Supporting cutoff of 1.56×10⁻⁶; PM2_Supporting not applied."
        ),
        "present_in_gnomad": True,
        "groupmax_faf": faf,
    }