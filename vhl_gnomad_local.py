from typing import Tuple, Optional, Sequence
from cyvcf2 import VCF

class LocalGnomadVHL:
    """
    Lightweight wrapper around a VHL-only gnomAD v4.x VCF.

    Assumes VCF has INFO fields including at least:
      - AC  (array, per ALT allele)
      - AN
      - some per-allele FAF95 popmax-like field (see header; common names
        include FAF95_POPMAX, FAF95_popmax, faf95_popmax, FAF95).
    """

    def __init__(self, vcf_path: str):
        self.vcf_path = vcf_path
        self.vcf = VCF(vcf_path)

    @staticmethod
    def _get_per_allele_value(info_val, alt_index: int) -> Optional[float]:
        """
        Extract a per-allele float from an INFO value that may be:
          - a single number
          - a comma-separated string
          - a list/tuple
        """
        if info_val is None:
            return None

        # cyvcf2 often gives a list/tuple for per-allele values
        if isinstance(info_val, (list, tuple)):
            if 0 <= alt_index < len(info_val):
                try:
                    return float(info_val[alt_index])
                except (TypeError, ValueError):
                    return None
            return None

        # Comma-separated string
        if isinstance(info_val, str) and "," in info_val:
            parts = info_val.split(",")
            if 0 <= alt_index < len(parts):
                try:
                    return float(parts[alt_index])
                except (TypeError, ValueError):
                    return None

        # Single scalar value
        try:
            return float(info_val)
        except (TypeError, ValueError):
            return None

    def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> Tuple[bool, Optional[float]]:
        """
        Look up a variant by chrom (string), pos (1-based int), ref, alt.

        Returns:
            present_in_gnomad (bool)
            groupmax_faf (float or None)
        """
        # cyvcf2 expects bare chromosome name (no 'chr')
        chrom = chrom.replace("chr", "")

        region = f"{chrom}:{pos}-{pos}"
        for rec in self.vcf(region):
            # Exact REF match
            if rec.REF != ref:
                continue

            # ALT may be multi-allelic; find index of our alt
            try:
                alt_index = list(rec.ALT).index(alt)
            except ValueError:
                continue  # this record is different alt

            # AC can be scalar or per-allele
            ac_val = rec.INFO.get("AC")
            ac = self._get_per_allele_value(ac_val, alt_index)
            present = bool(ac is not None and ac > 0)

            # Try several possible FAF95 popmax field names
            faf = None
            faf_keys: Sequence[str] = (
                "FAF95_POPMAX",
                "FAF95_popmax",
                "faf95_popmax",
                "FAF95",
                "faf95",
            )
            for key in faf_keys:
                if key in rec.INFO:
                    faf = self._get_per_allele_value(rec.INFO.get(key), alt_index)
                    break

            return present, faf

        # No matching record at this position for this ref/alt
        return False, None
