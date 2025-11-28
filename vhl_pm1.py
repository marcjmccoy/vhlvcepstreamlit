"""
VHL VCEP PM1 / PM1_Supporting classifier.

PM1 (Moderate):
  - Putative missense variants that are known germline hotspots; and/or
  - Located in a key pVHL functional domain (63–192 AA); and/or
  - Somatic hotspots with ≥10 instances for the same AA in cancerhotspots.org,
    provided the codon is not a germline hotspot.

PM1_Supporting:
  - Putative missense variants seen in somatic databases, with 1–9 instances
    for the same AA in cancerhotspots.org (and not germline hotspots).

Only simple missense substitutions are eligible.
"""

import re

from vhl_hgvs import parse_vhl_hgvs


# Germline missense hotspot codons (Stebbins + Chiorean; keep synced to VCEP table).
GERMLINE_HOTSPOT_POSITIONS = {
    167, 162, 178, 98, 78, 86,   # Stebbins
    65, 76, 80, 88, 96, 112, 117,
    161, 170, 176,
}

# Somatic hotspot codons from cancerhotspots.org / Walsh (not treated as
# somatic if also in GERMLINE_HOTSPOT_POSITIONS). Values are tumor counts.
SOMATIC_HOTSPOT_COUNTS = {
    65: 5,
    68: 12,
    74: 11,
    80: 6,
    88: 15,
    89: 10,
    96: 4,
    111: 14,
    114: 13,
    115: 10,
    121: 11,
    135: 17,
    151: 10,
    158: 11,
    169: 10,
}

# Critical pVHL functional domain for PM1: 63–192 AA.
PM1_FUNCTIONAL_DOMAINS = [(63, 192)]


def _parse_protein_codon(hgvs_protein: str):
    """Extract the amino‑acid position from a protein HGVS string."""
    if not hgvs_protein:
        return None
    if not hgvs_protein.startswith("p."):
        hgvs_protein = "p." + hgvs_protein

    m = re.search(r"p\.[A-Za-z]{3}(\d+)[A-Za-z*=?]+", hgvs_protein)
    return int(m.group(1)) if m else None


def _is_simple_missense(hgvs_protein: str) -> bool:
    """
    Determine whether the protein HGVS represents a simple missense substitution.
    Excludes nonsense, synonymous, frameshift, and in‑frame indels.
    """
    if not hgvs_protein:
        return False
    if not hgvs_protein.startswith("p."):
        hgvs_protein = "p." + hgvs_protein

    if any(x in hgvs_protein for x in ["*", "Ter", "=", "fs", "del", "ins", "dup"]):
        return False

    return bool(re.fullmatch(r"p\.[A-Za-z]{3}\d+[A-Za-z]{3}", hgvs_protein))


def _residue_in_pm1_domain(codon: int) -> bool:
    """Return True if the codon is within any PM1‑relevant pVHL functional domain."""
    if codon is None:
        return False
    return any(start <= codon <= end for start, end in PM1_FUNCTIONAL_DOMAINS)


def classify_vhl_pm1(hgvs_full: str):
    """
    VHL VCEP PM1 / PM1_Supporting classifier.

    Returns:
        dict with keys:
            - strength: "PM1", "PM1_Supporting", or None
            - context: human‑readable explanation string
    """
    _, _, protein = parse_vhl_hgvs(hgvs_full)
    codon = _parse_protein_codon(protein)

    # Restrict PM1 to simple missense substitutions.
    if not _is_simple_missense(protein):
        return {
            "strength": None,
            "context": (
                "PM1 is specified only for simple missense variants in VHL; "
                "this protein HGVS does not represent a missense substitution."
            ),
        }

    if codon is None:
        return {
            "strength": None,
            "context": (
                "Could not parse an amino‑acid position from the protein HGVS; "
                "PM1 not applied."
            ),
        }

    # 1) Germline missense hotspots (Moderate).
    if codon in GERMLINE_HOTSPOT_POSITIONS:
        return {
            "strength": "PM1",
            "context": (
                f"Amino acid {codon} is a curated germline missense hotspot in VHL; "
                "assign PM1 (Moderate)."
            ),
        }

    # 2) Somatic hotspots from cancerhotspots.org / Walsh (not germline).
    som_count = SOMATIC_HOTSPOT_COUNTS.get(codon, 0)
    if som_count >= 10:
        return {
            "strength": "PM1",
            "context": (
                f"Amino acid {codon} is a somatic hotspot with ≥10 tumors reported "
                "in cancerhotspots.org/Walsh and is not a germline hotspot; "
                "assign PM1 (Moderate)."
            ),
        }
    if 0 < som_count < 10:
        return {
            "strength": "PM1_Supporting",
            "context": (
                f"Amino acid {codon} is a somatic hotspot with <10 tumors reported "
                "in cancerhotspots.org/Walsh and is not a germline hotspot; "
                "assign PM1_Supporting."
            ),
        }

    # 3) Critical functional domains (63–192 AA).
    if _residue_in_pm1_domain(codon):
        return {
            "strength": "PM1",
            "context": (
                f"Amino acid {codon} lies within a key functional domain of pVHL "
                "(AA 63–192) that is enriched for pathogenic missense and lacks "
                "benign variation; assign PM1 (Moderate)."
            ),
        }

    # 4) No PM1 evidence.
    return {
        "strength": None,
        "context": (
            f"Amino acid {codon} is not in the curated germline or somatic hotspot "
            "lists and lies outside defined key functional domains; PM1 not applied."
        ),
    }