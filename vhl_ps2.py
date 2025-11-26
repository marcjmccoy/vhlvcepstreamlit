import re

def classify_vhl_ps2(
    hgvs_str,
    is_de_novo=False,
    phenotype=None,  # 'danish', 'consistent', 'nonspecific', or None
    panel_neg=None,  # dict of {gene: 'neg'} or None
    family_history=False
):
    panel_genes = ["SDHB", "SDHC", "SDHD", "RET", "MAX", "FH", "TMEM127", "NF1", "SDHA", "SDHAF2", "VHL"]

    # Variant type recognition
    is_missense = re.search(r'\(p\.[A-Za-z]{3}\d+[A-Za-z]{3}\)', hgvs_str) is not None
    is_nonsense = "Ter" in hgvs_str or "*" in hgvs_str or "stop" in hgvs_str
    is_frameshift = "fs" in hgvs_str
    is_splice = "+" in hgvs_str or "-" in hgvs_str or re.search(r'splice', hgvs_str)
    is_inframe_indel = "del" in hgvs_str and not "fs" in hgvs_str or "dup" in hgvs_str and not "fs" in hgvs_str

    # Eligibility checks
    if not is_de_novo:
        return {
            "strength": None,
            "context": "PS2 evidence applies only when both maternity and paternity of the proband are confirmed, indicating a true de novo event."
        }
    if family_history:
        return {
            "strength": None,
            "context": "PS2 cannot be scored because there is family history of VHL disease."
        }
    if not (is_missense or is_nonsense or is_frameshift or is_splice or is_inframe_indel):
        return {
            "strength": None,
            "context": "This variant type (not missense, nonsense, frameshift, splice, or in-frame indel) is not eligible for PS2 scoring."
        }

    score = 0
    context_detail = []

    # Panel negativity handling
    if panel_neg:
        neg_list = [gene for gene, result in panel_neg.items() if result == "neg"]
        missing_list = [gene for gene in panel_genes if gene not in panel_neg or panel_neg[gene] != "neg"]
        neg_panel_str = f"Panel genes negative: {', '.join(neg_list)}." if neg_list else "No panel genes were marked as negative."

    else:
        neg_panel_str = "No panel testing information provided."
        missing_list = panel_genes

    # Phenotype scoring and context
    if phenotype == "danish":
        score = 2
        context_detail.append(
            "Phenotype is highly specific for VHL by Danish criteria, which gives strong evidence for a pathogenic de novo VHL variant."
        )
        context_detail.append(
            "No family history, maternity and paternity confirmed, and variant type is eligible. "
            f"{neg_panel_str} This meets ACMG VCEP SVI criteria for 'Strong' strength (2 points)."
        )
    elif phenotype == "consistent":
        # RCC or Pheo—require panel
        if panel_neg:
            # For Pheo, SDHB/C/D and RET required; For RCC, SDHB/C/D required
            # Give credit for all negatives, but only full if all required are negative
            pheo_genes = {'SDHB', 'SDHC', 'SDHD', 'RET'}
            rcc_genes = {'SDHB', 'SDHC', 'SDHD'}
            missing_pheo = [g for g in pheo_genes if g not in neg_list]
            missing_rcc = [g for g in rcc_genes if g not in neg_list]

            if not missing_pheo or not missing_rcc:
                score = 1
                context_detail.append(
                    "Phenotype is consistent with VHL (such as isolated RCC or pheochromocytoma), and required panel genes for differential diagnosis are negative "
                    f"[Missing for RCC: {', '.join(missing_rcc) if missing_rcc else 'none'}; Missing for Pheo: {', '.join(missing_pheo) if missing_pheo else 'none'}]."
                )
                context_detail.append(
                    f"{neg_panel_str} This satisfies the VCEP criteria for 'Moderate' strength (1 point) with consistent phenotype and sufficient panel exclusion."
                )
            else:
                score = 0.5
                context_detail.append(
                    f"Phenotype is consistent, but panel testing is incomplete (the following required genes were not shown negative: {', '.join((missing_rcc if missing_rcc else []) + (missing_pheo if missing_pheo else [])}). "
                    f"{neg_panel_str} Evidence downgraded to 'Supporting' (0.5 points) due to incomplete panel testing."
                )
        else:
            score = 1
            context_detail.append(
                "Phenotype is consistent with VHL, but no panel testing information is available. Strength assigned 'Moderate' (1 point) per ACMG guidance."
            )
    elif phenotype == "nonspecific":
        score = 0.5
        context_detail.append(
            "Phenotype is not highly specific for VHL, lacking strong indication for VHL-associated features; therefore, evidence is supporting only ('Supporting', 0.5 points)."
        )
        context_detail.append(neg_panel_str)
    else:
        score = 0
        context_detail.append(
            "Phenotype does not fit established VHL criteria; PS2 cannot be scored."
        )

    # Strength conversion
    if score >= 4:
        strength = "PS2_VeryStrong"
        context_detail.append(
            "Multiple probands required to reach 'Very Strong' level (≥4 points); single proband scores cannot achieve this."
        )
    elif score >= 2:
        strength = "PS2"
        context_detail.append(
            "'Strong' evidence granted: Danish phenotype and negative testing (2 points)."
        )
    elif score >= 1:
        strength = "PS2_Moderate"
        context_detail.append(
            "'Moderate' evidence granted: Consistent phenotype and negative panel testing (1 point)."
        )
    elif score >= 0.5:
        strength = "PS2_Supporting"
        context_detail.append(
            "'Supporting' evidence: Phenotype less specific or panel testing incomplete (0.5 points)."
        )
    else:
        strength = None
        context_detail.append(
            "Score did not reach any PS2 evidence threshold for VHL."
        )

    return {
        "strength": strength,
        "context": " ".join(context_detail) + f" (Score: {score})"
    }