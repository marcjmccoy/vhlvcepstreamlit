import re

############################# PS2 ###################################################
def classify_vhl_ps2(
    hgvs_str,
    is_de_novo=False,  # True if both maternity and paternity are confirmed
    phenotype=None,    # 'danish', 'consistent', 'nonspecific', or None
    panel_neg=None,    # dict of {gene: 'neg'} or None if not done
    family_history=False
):
    """
    Classifies PS2 (VHL: de novo) for variant per ACMG VCEP rules.
    Args:
        hgvs_str: string, HGVS format ("NM_000551.4(VHL):c.191G>C (p.Arg64Pro)")
        is_de_novo: bool, True only if both maternity and paternity confirmed
        phenotype: string, one of 'danish', 'consistent', 'nonspecific', or None
        panel_neg: dict, e.g. {'SDHB': 'neg', 'SDHC': 'neg', ...}
        family_history: bool, VHL family history (False if none)
    Returns: dict with 'strength', 'context'
    """

    # Variant type logic
    is_missense = re.search(r'\(p\.[A-Za-z]{3}\d+[A-Za-z]{3}\)', hgvs_str) is not None
    is_nonsense = "Ter" in hgvs_str or "*" in hgvs_str or "stop" in hgvs_str
    is_frameshift = "fs" in hgvs_str
    is_splice = "+" in hgvs_str or "-" in hgvs_str or re.search(r'splice', hgvs_str)
    is_inframe_indel = "del" in hgvs_str and not "fs" in hgvs_str or "dup" in hgvs_str and not "fs" in hgvs_str
    is_synonymous = "synonymous" in hgvs_str or "silent" in hgvs_str

    # PS2 applies only to de novo (with confirmed maternity and paternity), no family history
    if not is_de_novo or family_history:
        return {
            "strength": None,
            "context": "PS2 only applies to confirmed de novo (no family history)."
        }

    # Phenotype consistency scoring
    score = 0
    context_detail = []
    if phenotype == "danish":
        score = 2
        context_detail.append("Phenotype highly specific by Danish criteria.")
    elif phenotype == "consistent":
        # Check panel negative for relevant genes if RCC/Pheo
        if panel_neg:
            if ("SDHB" in panel_neg and "SDHC" in panel_neg and
                "SDHD" in panel_neg and all(panel_neg.get(gene, 'neg') == 'neg'
                                            for gene in ['SDHB', 'SDHC', 'SDHD'])):
                score = 1
                context_detail.append("Phenotype consistent and panel required genes negative (RCC/Pheo).")
            else:
                score = 0.5
                context_detail.append("Phenotype consistent but incomplete negative panel.")
        else:
            score = 1
            context_detail.append("Phenotype consistent (no panel info available).")
    elif phenotype == "nonspecific":
        score = 0.5
        context_detail.append("Phenotype not highly specific.")
    else:
        score = 0
        context_detail.append("Phenotype does not fit specified categories.")

    # Final assignment per SVI guidance: max 1 pt per proband for confirmed de novo
    if score >= 4:
        strength = "PS2_VeryStrong"
    elif score >= 2:
        strength = "PS2"
    elif score >= 1:
        strength = "PS2_Moderate"
    elif score >= 0.5:
        strength = "PS2_Supporting"
    else:
        strength = None

    # If the variant is not relevant for PS2
    if not (is_missense or is_nonsense or is_frameshift or is_splice):
        return {
            "strength": None,
            "context": "Variant type is not eligible for PS2 scoring."
        }

    return {
        "strength": strength,
        "context": "; ".join(context_detail) + f" (Score: {score})"
    }