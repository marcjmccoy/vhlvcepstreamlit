import re
import logging
from typing import Optional, Tuple

import requests

logger = logging.getLogger(__name__)

# BS1 threshold: 0.00156% GroupMax FAF in gnomAD v4
GNOMAD_BS1_MIN_FAF: float = 0.0000156  # 0.00156%

# Reuse these from your BA1 module or keep them here if BS1 is separate
ENSEMBL_REST_SERVER = "https://rest.ensembl.org"
ENSEMBL_VEP_HGVS_EXT = "/vep/human/hgvs"
GNOMAD_GRAPHQL_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET_ID = "gnomad_r4"


def _safe_float(x) -> Optional[float]:
    try:
        if x is None:
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def _normalize_vhl_hgvs(hgvs_full: str) -> str:
    """
    Normalize full VHL HGVS like:
      'NM_000551.4(VHL):c.1A>T (p.Met1Leu)'
    to a clean transcript+cDNA string:
      'NM_000551.4:c.1A>T'
    """
    if not isinstance(hgvs_full, str):
        return hgvs_full

    core = re.sub(r"\s*\(p\.[^)]+\)\s*$", "", hgvs_full)

    m_tx = re.match(r"^([A-Z0-9_.]+)\(VHL\)", core)
    if not m_tx:
        m_tx = re.match(r"^([A-Z0-9_.]+):", core)
    transcript = m_tx.group(1) if m_tx else "NM_000551.4"

    m_c = re.search(r"(c\.[^ \)]+)", core)
    cdna = m_c.group(1) if m_c else None

    if cdna is None:
        return core
    return f"{transcript}:{cdna}"


def _resolve_hgvs_to_grch38(hgvs_tx: str) -> Optional[Tuple[str, int, str, str]]:
    """
    Resolve transcript HGVS (e.g. 'NM_000551.4:c.1A>T') to GRCh38 genomic
    coordinates: (chrom, pos, ref, alt) using Ensembl VEP /hgvs.
    """
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    data = {
        "hgvs_notations": [hgvs_tx],
        "assembly_name": "GRCh38",
    }

    try:
        resp = requests.post(
            ENSEMBL_REST_SERVER + ENSEMBL_VEP_HGVS_EXT,
            headers=headers,
            json=data,
            timeout=20,
        )
        resp.raise_for_status()
        decoded = resp.json()
    except Exception as exc:
        logger.warning("Ensembl VEP HGVS resolve failed for %s: %s", hgvs_tx, exc)
        return None

    try:
        rec_list = decoded
        if not rec_list:
            return None
        rec = rec_list[0]

        if rec.get("assembly_name") != "GRCh38":
            logger.warning(
                "HGVS %s resolved to %s, not GRCh38",
                hgvs_tx,
                rec.get("assembly_name"),
            )

        chrom = str(rec["seq_region_name"])
        pos = int(rec["start"])
        allele_string = rec["allele_string"]
        ref, alt = allele_string.split("/")

        return chrom, pos, ref, alt
    except Exception as exc:
        logger.warning("Unexpected VEP response for %s: %s", hgvs_tx, exc)
        return None


def _lookup_gnomad_v4_grch38(
    chrom: str, pos: int, ref: str, alt: str
) -> Tuple[bool, Optional[float]]:
    """
    Look up a GRCh38 SNV/indel in gnomAD v4 via the gnomAD browser GraphQL
    endpoint and return:
        (present_in_gnomad, groupmax_faf)

    Uses exome/genome faf95.popmax (FAF95 popmax). Intended for light use.
    """
    # gnomAD v4 variantId format is "3-10142007-A-T" (no "chr" prefix).
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"

    query = """
    query VariantQuery($variantId: String!, $datasetId: DatasetId!) {
      variant(variantId: $variantId, dataset: $datasetId) {
        variantId
        exome {
          ac
          an
          faf95 {
            popmax
          }
        }
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
        "variantId": variant_id,
        "datasetId": GNOMAD_DATASET_ID,
    }

    try:
        resp = requests.post(
            GNOMAD_GRAPHQL_URL,
            json={"query": query, "variables": variables},
            headers={"Content-Type": "application/json"},
            timeout=20,
        )
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:
        logger.warning("gnomAD GraphQL error for %s: %s", variant_id, exc)
        return False, None

    if "errors" in data:
        logger.info("gnomAD GraphQL errors for %s: %s", variant_id, data["errors"])
        return False, None

    v = (data.get("data") or {}).get("variant")
    if v is None:
        # variant absent
        return False, None

    exome = v.get("exome") or {}
    genome = v.get("genome") or {}

    # Determine presence from AC
    ac_exome = exome.get("ac")
    ac_genome = genome.get("ac")
    present = any(isinstance(x, int) and x > 0 for x in (ac_exome, ac_genome))

    # Collect FAF95 popmax values
    faf_candidates = []

    exome_faf = (exome.get("faf95") or {}).get("popmax")
    genome_faf = (genome.get("faf95") or {}).get("popmax")

    faf_candidates.append(_safe_float(exome_faf))
    faf_candidates.append(_safe_float(genome_faf))

    faf_values = [f for f in faf_candidates if f is not None]

    if not present:
        # no alternate alleles observed even if record exists
        return False, None

    if not faf_values:
        # present but no FAF95 popmax reported
        return True, None

    groupmax_faf = max(faf_values)
    return True, groupmax_faf


def classify_vhl_bs1(hgvs_full: str):
    """
    BS1 classifier using GRCh38-anchored, gnomAD v4 lookup.

    VHL VCEP rule:
      - Apply BS1 (Strong) if variant is present in gnomAD v4 AND
        GroupMax FAF (faf95.popmax) >= 0.0000156 (0.00156%).
      - Otherwise, BS1 not applied.
    """
    # Normalize HGVS to transcript:cDNA
    hgvs_tx = _normalize_vhl_hgvs(hgvs_full)

    # Resolve to GRCh38 genomic coordinates via Ensembl VEP /hgvs
    resolved = _resolve_hgvs_to_grch38(hgvs_tx)
    if resolved is None:
        return {
            "strength": None,
            "context": (
                f"Could not resolve {hgvs_tx} to GRCh38 genomic coordinates; "
                "BS1 not applied."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    chrom, pos, ref, alt = resolved

    # Look up the GRCh38 variant in gnomAD v4 via GraphQL
    present, faf = _lookup_gnomad_v4_grch38(chrom, pos, ref, alt)

    # Variant completely absent from gnomAD v4: cannot satisfy BS1
    if present is False:
        return {
            "strength": None,
            "context": (
                f"{hgvs_full} ({chrom}-{pos}-{ref}-{alt}, GRCh38) is absent from "
                "gnomAD v4; BS1 not applied."
            ),
            "present_in_gnomad": False,
            "groupmax_faf": None,
        }

    # Present but no FAF reported: cannot test BS1 threshold
    if present is True and faf is None:
        return {
            "strength": None,
            "context": (
                f"{hgvs_full} ({chrom}-{pos}-{ref}-{alt}, GRCh38) is present in "
                "gnomAD v4 but has no GroupMax Filtering Allele Frequency; "
                "BS1 not applied."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": None,
        }

    # Present with FAF: check BS1 cutoff
    if faf is not None and faf >= GNOMAD_BS1_MIN_FAF:
        return {
            "strength": "BS1",  # or "BS1_Strong" if you encode strength level
            "context": (
                f"{hgvs_full} ({chrom}-{pos}-{ref}-{alt}, GRCh38) has GroupMax FAF "
                f"≈{faf:.3e} in gnomAD v4, ≥1.56×10⁻⁵ (0.00156%); BS1 applied."
            ),
            "present_in_gnomad": True,
            "groupmax_faf": faf,
        }

    # Present but below BS1 FAF threshold
    return {
        "strength": None,
        "context": (
            f"{hgvs_full} ({chrom}-{pos}-{ref}-{alt}, GRCh38) has GroupMax FAF "
            f"≈{faf:.3e} in gnomAD v4, below 1.56×10⁻⁵ (0.00156%); BS1 not applied."
        ),
        "present_in_gnomad": True,
        "groupmax_faf": faf,
    }