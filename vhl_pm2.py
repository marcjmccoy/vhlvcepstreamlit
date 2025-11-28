"""
VHL VCEP PM2_Supporting classifier.

Implements the VHL VCEP rule using gnomAD v4 on GRCh38:

- PM2_Supporting can be applied for variants either absent from gnomAD or
  with <= 0.00000156 (0.000156%) GroupMax Filtering Allele Frequency
  (GroupMax FAF) in gnomAD v4.
- If no GroupMax FAF is calculated (e.g. single observation), PM2_Supporting
  may also be applied.

This module:
  1) Parses VHL transcript HGVS input, e.g. 'NM_000551.4(VHL):c.1A>T (p.Met1Leu)'.
  2) Uses MyVariant.info to map NM_000551.4:c.* to a GRCh38 genomic HGVS.
  3) Queries the gnomAD v4 GraphQL API at that locus for AC and FAF95.popmax.
  4) Applies the PM2_Supporting rule.
"""

import re
import logging
import requests

logger = logging.getLogger(__name__)

# VHL VCEP PM2_Supporting GroupMax FAF threshold (gnomAD v4)
GNOMAD_PM2_MAX_FAF: float = 0.00000156  # 0.000156%

# gnomAD GraphQL endpoint (v4 on GRCh38)
GNOMAD_API_BASE = "https://gnomad.broadinstitute.org/api"

# MyVariant.info base URL
MYVARIANT_BASE = "https://myvariant.info/v1"


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
# MyVariant: c.HGVS -> GRCh38 genomic HGVS
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


def _parse_genomic_hgvs_grch38(hgvs_g: str):
    """
    Parse GRCh38 genomic HGVS like 'chr3:g.10191340A>T' into
    chrom='3', pos=10191340, ref='A', alt='T'.
    """
    m = re.match(r"^chr(\w+):g\.(\d+)([ACGT])>([ACGT])$", hgvs_g)
    if not m:
        return None, None, None, None
    chrom = m.group(1)
    pos = int(m.group(2))
    ref = m.group(3)
    alt = m.group(4)
    return chrom, pos, ref, alt


def _get_grch38_coords_from_cdna(hgvs_tx: str):
    """
    Use MyVariant.info to map a transcript HGVS (NM_000551.4:c.*)
    to a GRCh38 genomic HGVS and return chrom, pos, ref, alt.

    Returns:
        chrom, pos, ref, alt or (None, None, None, None) if mapping fails.
    """
    try:
        resp = requests.get(
            f"{MYVARIANT_BASE}/query",
            params={
                "q": hgvs_tx,
                "fields": "hg38",
                "assembly": "hg38",
                "size": 1,
            },
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json() or {}
    except Exception as exc:
        logger.warning("MyVariant hg38 mapping failed for %s: %s", hgvs_tx, exc)
        return None, None, None, None

    hits = data.get("hits") or []
    if not hits:
        return None, None, None, None

    hg38_block = hits[0].get("hg38") or {}
    hgvs_g = hg38_block.get("hgvs")
    if not hgvs_g:
        return None, None, None, None

    return _parse_genomic_hgvs_grch38(hgvs_g)


# ----------------------------------------------------------------------
# gnomAD v4 (GraphQL) query by GRCh38 coordinates
# ----------------------------------------------------------------------


def _gnomad_query_variant_grch38(chrom: str, pos: int, ref: str, alt: str):
    """
    Query gnomAD GraphQL API for a GRCh38 variant and return
    presence and GroupMax FAF (if available) from the v4 dataset.

    Returns:
        present_in_gnomad (bool)
        groupmax_faf (float or None)
    """
    query = """
    query Variant($chrom: String!, $pos: Int!, $ref: String!, $alt: String!) {
      variant(dataset: gnomad_r4, variantId: $chrom:$pos:$ref:$alt) {
        genome {
          ac
          an
          faf95 {
            popmax
          }
        }
      }
    }
    """

    variables = {
        "chrom": str(chrom),
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
    }

    try:
        resp = requests.post(
            GNOMAD_API_BASE,
            json={"query": query, "variables": variables},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json() or {}
    except Exception as exc:
        logger.warning(
            "gnomAD API error for %s:%s:%s>%s: %s", chrom, pos, ref, alt, exc
        )
        return False, None

    variant = (data.get("data") or {}).get("variant")
    if not variant:
        # No record in v4 – treat as absent.
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
    VHL VCEP PM2_Supporting classifier using gnomAD v4 via:
      - MyVariant.info (for GRCh38 mapping), then
      - gnomAD GraphQL API (for AC and GroupMax FAF).

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

    tx_hgvs = f"{transcript}:{cdna}"  # e.g. NM_000551.4:c.1A>T

    # 2) Map to GRCh38 genomic coordinates via MyVariant.info.
    chrom, pos, ref, alt = _get_grch38_coords_from_cdna(tx_hgvs)
    if not (chrom and pos and ref and alt):
        return {
            "strength": None,
            "context": (
                f"Could not map {tx_hgvs} to GRCh38 genomic coordinates via MyVariant.info; "
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
                f"(no alternate alleles observed at chr{chrom}:g.{pos}{ref}>{alt}), "
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
                f"{hgvs_full} is present in gnomAD v4 at chr{chrom}:g.{pos}{ref}>{alt} "
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
                f"{hgvs_full} is observed in gnomAD v4 at chr{chrom}:g.{pos}{ref}>{alt} "
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
            f"chr{chrom}:g.{pos}{ref}>{alt}, which exceeds the "
            "PM2_Supporting cutoff of 1.56×10⁻⁶; PM2_Supporting not applied."
        ),
        "present_in_gnomad": True,
        "groupmax_faf": faf,
    }