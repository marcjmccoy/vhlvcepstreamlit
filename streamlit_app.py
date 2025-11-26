import streamlit as st
import re
from vhl_pvs1 import classify_vhl_pvs1  # Import PVS1 classifier
from vhl_ps1 import classify_vhl_ps1    # Import PS1 classifier
from vhl_ps2 import classify_vhl_ps2    # Import PS2 classifier

# ---- Sidebar: visually enhanced rules for combining ----
st.sidebar.markdown("""
<style>
    .sidebar-title {font-size: 22px; font-weight: bold; color: #155a90;}
    .category {font-size: 17px; font-weight: bold; margin-top: 12px;}
    .combo-rule {margin-bottom: 6px;}
    .rule-label {font-weight: bold; color: #2c2c2c;}
    .pathogenic {background: #ffe4e4; border-radius: 7px; padding: 7px;}
    .likely_pathogenic {background: #fff3d6; border-radius: 7px; padding: 7px;}
    .benign {background: #e8fae8; border-radius: 7px; padding: 7px;}
    .likely_benign {background: #f8fbee; border-radius: 7px; padding: 7px;}
</style>
<div class='sidebar-title'>🧩 Rules for Combining Criteria</div>

<div class='category pathogenic'>🔴 <b>Pathogenic</b>
<ul style="margin-bottom:6px;">
    <li class='combo-rule'><span class='rule-label'>1 Very Strong</span> <small>(PVS1, PS2_VeryStrong, PS4_VeryStrong)</small> <br/>
        <b>AND</b> <span class='rule-label'>1 Strong</span> <small>(PS1, PS2, PS4, PP1_Strong)</small></li>
    <li class='combo-rule'><span class='rule-label'>1 Very Strong</span> <b>AND</b> <span class='rule-label'>2 Moderate</span> <small>(PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)</small></li>
    <li class='combo-rule'><span class='rule-label'>1 Very Strong</span> <b>AND</b> <span class='rule-label'>1 Moderate</span> <b>AND</b> <span class='rule-label'>1 Supporting</span> <small>(PS2_Supporting, PS3_Supporting, PS4_Supporting, PM1_Supporting, PM2_Supporting, PP1, PP3)</small></li>
    <li class='combo-rule'><span class='rule-label'>1 Very Strong</span> <b>AND</b> <span class='rule-label'>2 Supporting</span></li>
    <li class='combo-rule'><span class='rule-label'>&ge; 2 Strong</span></li>
    <li class='combo-rule'><span class='rule-label'>1 Strong</span> <b>AND</b> <span class='rule-label'>&ge; 3 Moderate</span></li>
    <li class='combo-rule'><span class='rule-label'>1 Strong</span> <b>AND</b> <span class='rule-label'>2 Moderate</span> <b>AND</b> <span class='rule-label'>2 Supporting</span></li>
    <li class='combo-rule'><span class='rule-label'>1 Strong</span> <b>AND</b> <span class='rule-label'>1 Moderate</span> <b>AND</b> <span class='rule-label'>&ge; 4 Supporting</span></li>
</ul>
</div>

<div class='category likely_pathogenic'>🟠 <b>Likely Pathogenic</b>
<ul style="margin-bottom:6px;">
    <li class='combo-rule'><span class='rule-label'>1 Very Strong</span> <b>AND</b> <span class='rule-label'>1 Moderate</span></li>
    <li class='combo-rule'><span class='rule-label'>1 Strong</span> <b>AND</b> <span class='rule-label'>&ge; 2 Supporting</span></li>
    <li class='combo-rule'><span class='rule-label'>&ge; 3 Moderate</span></li>
    <li class='combo-rule'><span class='rule-label'>2 Moderate</span> <b>AND</b> <span class='rule-label'>&ge; 2 Supporting</span></li>
    <li class='combo-rule'><span class='rule-label'>1 Moderate</span> <b>AND</b> <span class='rule-label'>&ge; 4 Supporting</span></li>
</ul>
</div>

<div class='category benign'>🟢 <b>Benign</b>
<ul style="margin-bottom:6px;">
    <li class='combo-rule'><span class='rule-label'>&ge; 2 Strong</span> <small>(BS1, BS2, BS4, BP2_Strong)</small></li>
    <li class='combo-rule'><span class='rule-label'>1 Stand Alone</span> <small>(BA1)</small></li>
</ul>
</div>

<div class='category likely_benign'>🟡 <b>Likely Benign</b>
<ul style="margin-bottom:6px;">
    <li class='combo-rule'><span class='rule-label'>1 Strong</span> <b>AND</b> <span class='rule-label'>1 Supporting</span></li>
    <li class='combo-rule'><span class='rule-label'>&ge; 2 Supporting</span></li>
</ul>
</div>
""", unsafe_allow_html=True)

# ---- Main panel ----
st.title("VHL/VCEP Classifier (v1.2 TEST) 🧬")

st.markdown(
    """
    Enter a variant in HGVS (cDNA) notation to predict PVS1, PS1, and PS2 strength according to VHL VCEP guidelines. Classifiers for additional criteria are under development. Context explains why each strength label was applied. Refer to the original criteria [here](https://cspec.genome.network/cspec/ui/svi/doc/GN078).
    <br/><b style='color:red;'>Note:</b> PM3, PP2, PP4, PP5, BP1, BP6 are <b>not recommended</b> for use by the ClinGen Sequence Variant Interpretation VCEP Review Committee for VHL.
    """, unsafe_allow_html=True
)

st.markdown(
    """
    **Example variants:**  
    NM_000551.4(VHL):c.263+1G>A  
    NM_000551.4(VHL):c.119_340del  
    NM_000551.4(VHL):c.123+5G>A  
    NM_000551.4(VHL):c.1A>T (p.Met1Leu)  
    NM_000551.4(VHL):c.160A>T (p.Met54Leu)  
    NM_000551.4(VHL):c.150G>A (p.Pro50=)  
    NM_000551.4(VHL):c.408del (p.Phe136fs)  
    NM_000551.4(VHL):c.120_122del (p.Val41del)  
    NM_000551.4(VHL):c.190_195dup (p.Lys64_Leu65dup)  
    NM_000551.4(VHL):c.263G>A (p.Trp88Ter)
    """
)

hgvs_input = st.text_input(
    "HGVS c. Notation (see examples above)", ""
)

with st.expander("PVS1/PS1 Options"):
    st.markdown("### PVS1 Options")
    exon_skipping = st.checkbox("Exon Skipping", value=False)
    cryptic_disrupts_rf = st.radio(
        "Cryptic Splice Disrupts Reading Frame?", (None, True, False))
    cryptic_preserves_rf = st.radio(
        "Cryptic Splice Preserves Reading Frame?", (None, True, False))
    duplication_type = st.selectbox(
        "Duplication Type", (None, "tandem", "not_in_tandem"))

with st.expander("PS2 (De Novo) Options"):
    st.markdown(
        """
        **Scoring Clarity for PS2:**  
        - If <b>ANY family history</b> of VHL disease is present, PS2 cannot be assigned and results will indicate this restriction.
        """, unsafe_allow_html=True
    )
    family_history = st.checkbox(
        "Any family history of VHL?",
        value=False,
        help="Check if any features of VHL or known pathogenic variant(s) present in family. PS2 cannot be scored in such cases."
    )
    if family_history:
        st.warning("PS2 cannot be assigned if any family history of VHL disease is present. This overrides all other options.")
    else:
        st.markdown(
            """
            To achieve <b>Moderate (1 point)</b> PS2 evidence:  
            1. Confirmed de novo status (both maternity and paternity testing done)  
            2. For <i>Consistent</i> or <i>RCC+Pheo</i> phenotypes, all required genes must be negative.  
            <br/>
            <b>If panel is incomplete or untested, only Supporting (0.5 points) can be assigned, even with confirmed de novo status.</b>
            """, unsafe_allow_html=True
        )
        is_de_novo = st.checkbox(
            "Confirmed de novo (both maternity and paternity tested)",
            value=False,
            help="Check only if BOTH maternity and paternity testing confirm de novo status."
        )

        st.markdown("**Select phenotype category:**")
        phenotype = st.selectbox(
            "Phenotype category",
            (None, "danish", "consistent", "nonspecific"),
            format_func=lambda x: {
                None: "Not specified",
                "danish": "Highly specific (Danish/International criteria)",
                "consistent": "Consistent with gene (e.g. VHL2c, RCC+Pheo as per specs)",
                "nonspecific": "Nonspecific phenotype"
            }.get(x, x)
        )

        st.markdown(
            """
            <b>Panel Results (should be negative to maximize evidence):</b>  
            • For VHL2c (pheo only): All [MAX, NF1, RET, SDHA, SDHB, SDHC, SDHD, SDHAF2, TMEM127, VHL] should be marked negative.  
            • For RCC+Pheo: All [MAX, FH, SDHA, SDHB, SDHC, SDHD, SDHAF2, TMEM127] should be negative.  
            <br/>If incomplete, only PS2_Supporting can be assigned.
            """, unsafe_allow_html=True
        )

        st.markdown("**Mark negative testing results below only if confirmed:**")
        panel_genes = [
            "SDHB", "SDHC", "SDHD", "RET", "MAX", "FH", "TMEM127",
            "NF1", "SDHA", "SDHAF2", "VHL"
        ]
        panel_neg = {}
        cols = st.columns(4)
        for i, gene in enumerate(panel_genes):
            with cols[i % 4]:
                panel_neg[gene] = "neg" if st.checkbox(
                    f"{gene} negative", value=False, key=gene
                ) else None

# ----- Run classifiers and show results -----
if hgvs_input:
    pvs1_result = classify_vhl_pvs1(
        hgvs_input,
        exon_skipping=exon_skipping,
        cryptic_disrupts_rf=True if cryptic_disrupts_rf is True else None,
        cryptic_preserves_rf=True if cryptic_preserves_rf is True else None,
        duplication_type=duplication_type if duplication_type not in (None, "None") else None,
    )
    ps1_result = classify_vhl_ps1(hgvs_input)
    ps2_result = classify_vhl_ps2(
        hgvs_input,
        is_de_novo=is_de_novo if not family_history else False,
        phenotype=phenotype if not family_history and phenotype is not None else None,
        panel_neg={gene: result for gene, result in panel_neg.items() if result == "neg"} if not family_history else {},
        family_history=family_history
    )
    st.header("Classification Results")
    st.subheader("PVS1 Result")
    st.write(pvs1_result)
    st.subheader("PS1 Result")
    st.write(ps1_result)
    st.subheader("PS2 Result")
    st.write(ps2_result)
