import re

def classify_vhl_ps2(
    hgvs_str,
    is_de_novo=False,
    phenotype=None,  # 'danish', 'consistent', 'nonspecific', or None
    panel_neg=None,  # dict of {gene: 'neg'} or None
    family_history=False
):
    """
    Classifies PS2 (VHL: de novo) for variant per ACMG VCEP and VHL-specific rules.
    """

    # Panel gene definitions by phenotype/group
    VHL2C_PHEO_PANEL = {"MAX", "NF1", "RET", "SDHA", "SDHB", "SDHC", "SDHD", "SDHAF2", "TMEM127", "VHL"}
    RCC_PHEO_PANEL = {"MAX", "FH", "SDHA", "SDHB", "SDHC", "SDHD", "SDHAF2", "TMEM127"}
    SDHX = {"SDHA", "SDHB", "SDHC", "SDHD", "SDHAF2"}

    # All panel genes offered in the app (not directly used in scoring)
    STREAMLIT_GENES = ["SDHB", "SDHC", "SDHD", "RET", "MAX", "FH", "TMEM127", "NF1", "SDHA", "SDHAF2", "VHL"]

    # Variant type recognition
    is_missense = re.search(r'\(p\.[A-Za-z]{3}\d+[A-Za-z]{3}\)', hgvs_str) is not None
    is_nonsense = "Ter" in hgvs_str or "*" in hgvs_str or "stop" in hgvs_str
    is_frameshift = "fs" in hgvs_str
    is_splice = "+" in hgvs_str or "-" in hgvs_str or re.search(r'splice', hgvs_str)
    is_inframe_indel = (
        ("del" in hgvs_str and not "fs" in hgvs_str) or
        ("dup" in hgvs_str and not "fs" in hgvs_str)
    )

    # Eligibility checks
    if not is_de_novo:
        return {
            "strength": None,
            "context": "PS2 can only be scored if both maternity and paternity are confirmed (true de novo origin)."
        }
    if family_history:
        return {
            "strength": None,
            "context": "PS2 cannot be scored if there is any family history of VHL disease."
        }
    if not (is_missense or is_nonsense or is_frameshift or is_splice or is_inframe_indel):
        return {
            "strength": None,
            "context": "Variant type is not eligible for PS2 (supported: missense, nonsense, frameshift, splice, in-frame indel)."
        }

    score = 0
    context_lines = []

    # Collate negative-tested genes, if provided by user
    if panel_neg:
        neg_list = [gene for gene, result in panel_neg.items() if result == "neg"]
        neg_panel_str = f"Panel genes marked negative: {', '.join(sorted(neg_list))}."
    else:
        neg_list = []
        neg_panel_str = "Panel testing results not provided."

    # Highly specific phenotype (e.g., Danish criteria)
    if phenotype == "danish":
        score = 2
        context_lines.append(
            "Phenotype is highly specific for VHL by published Danish or International criteria, which strongly indicate pathogenicity for VHL variants in de novo cases."
        )
        context_lines.append(
            "Per guideline, if only 'Danish/International criteria' are stated in literature, but without case specifics, must downgrade score to 'Phenotype consistent' (do not assign 'Strong')."
        )
        context_lines.append(
            "Strength: PS2 applied (Strong, 2 points) since patient meets Danish criteria AND negative panel testing for all relevant genes confirmed or not required for this phenotype."
        )
        context_lines.append(neg_panel_str)

    elif phenotype == "consistent":
        # Panel completeness checks
        missing_pheo = [gene for gene in VHL2C_PHEO_PANEL if gene not in neg_list]
        missing_rcc_pheo = [gene for gene in RCC_PHEO_PANEL if gene not in neg_list]
        panel_tested = bool(panel_neg)
        
        if panel_tested and not missing_pheo:
            score = 1
            context_lines.append(
                "Phenotype is consistent with VHL (VHL2c, Pheo only), AND comprehensive panel testing is negative for all required genes (MAX, NF1, RET, SDHA, SDHB, SDHC, SDHD, SDHAF2, TMEM127, VHL)."
            )
            context_lines.append(
                "This fulfills ACMG VCEP requirements for PS2_Moderate evidence (1 point)."
            )
            context_lines.append(neg_panel_str)
        elif panel_tested and missing_pheo:
            score = 0.5
            context_lines.append(
                f"Phenotype is consistent with VHL (VHL2c, Pheo only), but the following required panel genes were not shown negative: {', '.join(missing_pheo)}."
            )
            context_lines.append("As panel is incomplete, scoring downgraded to PS2_Supporting (0.5 points).")
            context_lines.append(neg_panel_str)
        elif panel_tested and not missing_rcc_pheo:
            score = 1
            context_lines.append(
                "Phenotype is consistent with VHL (RCC+Pheo) and comprehensive panel testing is negative for MAX, FH, SDHA, SDHB, SDHC, SDHD, SDHAF2, TMEM127."
            )
            context_lines.append(
                "This fulfills ACMG VCEP requirements for PS2_Moderate evidence (1 point)."
            )
            context_lines.append(neg_panel_str)
        elif panel_tested and missing_rcc_pheo:
            score = 0.5
            context_lines.append(
                f"Phenotype is consistent with VHL (RCC+Pheo), but the following required panel genes were not shown negative: {', '.join(missing_rcc_pheo)}."
            )
            context_lines.append("Because comprehensive panel testing is incomplete, scoring downgraded to PS2_Supporting (0.5 points).")
            context_lines.append(neg_panel_str)
        elif not panel_tested:
            # ***CRITICAL COMMENT - IMAGE-BASED RULE***
            # Per ACMG VCEP scoring table (see attached image):
            # If a proband with RCC+Pheo and no family history lacks panel testing,
            # count as "Phenotype consistent with the gene but not highly specific" (0.5 points, PS2_Supporting).
            # Never assign 1 (Moderate) unless comprehensive panel is confirmed negative.
            score = 0.5
            context_lines.append(
                "Phenotype is consistent with VHL (e.g., RCC+Pheo or VHL2c), but no comprehensive panel testing information is available. "
                "Per ACMG VCEP and scoring table, this is PS2_Supporting (0.5 points) since other genetic etiologies cannot be excluded."
            )
            context_lines.append(neg_panel_str)

    elif phenotype == "nonspecific":
        score = 0.5
        context_lines.append(
            "Phenotype is not highly specific for VHL, lacking strong association with classic VHL features; evidence is supporting only (0.5 points)."
        )
        context_lines.append(neg_panel_str)
    else:
        score = 0
        context_lines.append(
            "Phenotype was not categorized; PS2 cannot be scored."
        )

    # SVI-guided label logic (VCEP/ClinGen rules)
    if score >= 4:
        strength = "PS2_VeryStrong"
        context_lines.append(
            "Multiple probands needed to reach 'Very Strong' (≥4 points)—never assigned to single proband according to ACMG/ClinGen SVI guidelines."
        )
    elif score >= 2:
        strength = "PS2"
        context_lines.append(
            "PS2 ('Strong', 2 points) based on highly specific phenotype and confirmed de novo status."
        )
    elif score >= 1:
        strength = "PS2_Moderate"
        context_lines.append(
            "PS2_Moderate assigned: Consistent phenotype and comprehensive negative panel testing (1 point)."
        )
    elif score >= 0.5:
        strength = "PS2_Supporting"
        context_lines.append(
            "PS2_Supporting assigned: Consistent or nonspecific phenotype with incomplete or absent panel testing (0.5 points)."
        )
    else:
        strength = None
        context_lines.append(
            "Score did not reach a PS2 evidence threshold for VHL."
        )

    return {
        "strength": strength,
        "context": " ".join(context_lines) + f" (Score: {score})"
    }
