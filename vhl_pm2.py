import logging
import requests

logger = logging.getLogger(__name__)

# VHL VCEP PM2_Supporting GroupMax FAF threshold (gnomAD v4)
GNOMAD_PM2_MAX_FAF: float = 0.00000156  # 0.000156%

# gnomAD REST base for GRCh38 / v4.x
GNOMAD_API_BASE = "https://gnomad.broadinstitute.org/api"


def _safe_float(value):
    try:
        if value is None:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _gnomad_query_variant_grch38(chrom: str, pos: int, ref: str, alt: str):
    """
    Query gnomAD GraphQL API for a GRCh38 variant and return
    presence and GroupMax FAF (if available) from the v4 dataset.

    Args:
        chrom: chromosome string, e.g. '3'
        pos:   1-based genomic position (int)
        ref:   reference allele
        alt:   alternate allele

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
        "chrom": chrom,
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
        logger.warning("gnomAD API error for %s:%s:%s>%s: %s",
                       chrom, pos, ref, alt, exc)
        return False, None

    variant = (data.get("data") or {}).get("variant")
    if not variant:
        # No record in v4 – treat as absent.
        return False, None

    genome = variant.get("genome") or {}
    ac = genome.get("ac") or 0
    present = ac > 0

    faf_block = genome.get("faf95") or {}
    faf_popmax = _safe_float(faf_block.get("popmax"))

    return present, faf_popmax


def classify_vhl_pm2(
    hgvs_full: str,
    chrom: str | None = None,
    pos: int | None = None,
    ref: str | None = None,
    alt: str | None = None,
):
    """
    VHL VCEP PM2_Supporting classifier using gnomAD v4 via the gnomAD GraphQL API.

    The VCEP rule:
      - Apply PM2_Supporting if the variant is ABSENT from gnomAD v4; OR
      - If GroupMax Filtering Allele Frequency (FAF) ≤ 1.56×10⁻⁶; OR
      - If the variant is present but no GroupMax FAF is calculated.

    To honour this exactly, this function expects a GRCh38 genomic
    representation (chrom, pos, ref, alt) for querying gnomAD v4.
    `hgvs_full` is used only for explanatory context.

    Args:
        hgvs_full:
            VHL HGVS for display (e.g. 'NM_000551.4(VHL):c.1A>T').
        chrom, pos, ref, alt:
            GRCh38 genomic representation matching gnomAD v4, e.g.
            chrom='3', pos=10191340, ref='A', alt='T'.

    Returns:
        dict with keys:
            - strength: "PM2_Supporting" or None
            - context: explanation string
            - present_in_gnomad: bool
            - groupmax_faf: float or None
    """
    if not (chrom and pos and ref and alt):
        return {
            "strength": None,
            "context": (
                "No GRCh38 genomic coordinates were provided; cannot query gnomAD v4 "
                "for GroupMax FAF, so PM2_Supporting is not applied."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    present, faf = _gnomad_query_variant_grch38(chrom, pos, ref, alt)

    # 1) Absent from gnomAD v4.
    if not present:
        return {
            "strength": "PM2_Supporting",
            "context": (
                f"{hgvs_full} is absent from gnomAD v4 (no alternate alleles observed), "
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
                f"{hgvs_full} is present in gnomAD v4 with GroupMax FAF={faf:.3e}, "
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
                f"{hgvs_full} is observed in gnomAD v4 but lacks a GroupMax Filtering "
                "Allele Frequency estimate (e.g. single observation); this still "
                "qualifies for PM2_Supporting."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": None,
        }

    # 4) Present with FAF above threshold – PM2_Supporting not met.
    return {
        "strength": None,
        "context": (
            f"{hgvs_full} has GroupMax FAF={faf:.3e} in gnomAD v4, which exceeds the "
            "PM2_Supporting cutoff of 1.56×10⁻⁶; PM2_Supporting not applied."
        ),
        "present_in_gnomad": True,
        "groupmax_faf": faf,
    }