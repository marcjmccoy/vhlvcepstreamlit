import streamlit as st
import re
from vhl_pvs1 import classify_vhl_pvs1  # Import PVS1 classifier
from vhl_ps1 import classify_vhl_ps1    # Import PS1 classifier
from vhl_ps2 import classify_vhl_ps2    # Import PS2 classifier

st.set_page_config(
    page_title="VHL/VCEP Classifier",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("VHL/VCEP Classifier (v1.2 TEST) 🧬")

st.markdown(
    """
    Enter a VHL variant in HGVS (cDNA) notation to predict PVS1, PS1, and PS2 strengths per VHL VCEP guidelines. Classifiers for additional criteria are under development. Context explains why each strength label was applied.  
    *Criteria not recommended for VHL: PM3, PP2, PP4, PP5, BP1, BP6.*  
    Refer to original criteria [here](https://cspec.genome.network/cspec/ui/svi/doc/GN078).
    """
)

# Sidebar: Input & options
st.sidebar.header("Enter Variant")
examples = [
    "NM_000551.4(VHL):c.263+1G>A",
    "NM_000551.4(VHL):c.119_340del",
    "NM_000551.4(VHL):c.123+5G>A",
    "NM_000551.4(VHL):c.1A>T (p.Met1Leu)",
    "NM_000551.4(VHL):c.160A>T (p.Met54Leu)",
    "NM_000551.4(VHL):c.150G>A (p.Pro50=)",
    "NM_000551.4(VHL):c.408del (p.Phe136fs)",
    "NM_000551.4(VHL):c.120_122del (p.Val41del)",
    "NM_000551.4(VHL):c.263G>A (p.Trp88Ter)"
]

selected_example = st.sidebar.selectbox(
    "Try an example variant:", [""] + examples, help="Quick select for demonstration"
)
hgvs_input = st.sidebar.text_input(
    "HGVS c. Notation", value=selected_example, help="Paste or enter HGVS c. Notation"
)

st.sidebar.markdown("### PVS1 Options")
exon_skipping = st.sidebar.checkbox("Exon Skipping", value=False, help="Check if variant causes exon skipping.")
cryptic_disrupts_rf = st.sidebar.radio(
    "Cryptic Splice Disrupts Reading Frame?", (None, True, False), help="Set if cryptic splice disrupts frame."
)
cryptic_preserves_rf = st.sidebar.radio(
    "Cryptic Splice Preserves Reading Frame?", (None, True, False), help="Set if cryptic splice preserves frame."
)
duplication_type = st.sidebar.selectbox(
    "Duplication Type", (None, "tandem", "not_in_tandem"), help="Specify type if variant is a duplication."
)

st.sidebar.markdown("### PS2 (De Novo) Options")
family_history = st.sidebar.checkbox(
    "Any family history of VHL?", value=False,
    help="If true, PS2 evidence cannot be assigned."
)
is_de_novo = False
phenotype = None
panel_neg = {}

if not family_history:
    is_de_novo = st.sidebar.checkbox("Confirmed de novo (both maternity and paternity tested)", value=False)
    st.sidebar.markdown("Select phenotype category:")
    phenotype = st.sidebar.selectbox(
        "Phenotype category", (None, "danish", "consistent", "nonspecific"),
        format_func=lambda x: {
            None: "Not specified",
            "danish": "Highly specific (Danish/International criteria)",
            "consistent": "Consistent with gene (e.g. VHL2c, RCC+Pheo as per specs)",
            "nonspecific": "Nonspecific phenotype"
        }.get(x, x)
    )
    st.sidebar.markdown("Mark negative panels below only if confirmed:")
    panel_genes = [
        "SDHB", "SDHC", "SDHD", "RET", "MAX", "FH", "TMEM127",
        "NF1", "SDHA", "SDHAF2", "VHL"
    ]
    panel_neg = {
        gene: "neg" if st.sidebar.checkbox(
            f"{gene} negative", value=False, key=gene
        ) else None for gene in panel_genes
    }

# Main Results Area
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

    tab1, tab2, tab3 = st.tabs(["PVS1 Result", "PS1 Result", "PS2 Result"])
    with tab1:
        if 'Strong' in str(pvs1_result):
            st.success(f"PVS1: {pvs1_result}")
        else:
            st.info(f"PVS1: {pvs1_result}")
    with tab2:
        if 'Strong' in str(ps1_result):
            st.success(f"PS1: {ps1_result}")
        else:
            st.info(f"PS1: {ps1_result}")
    with tab3:
        if 'Strong' in str(ps2_result):
            st.success(f"PS2: {ps2_result}")
        elif family_history:
            st.error("PS2 cannot be assigned due to family history of VHL.")
        else:
            st.info(f"PS2: {ps2_result}")

    st.divider()  # horizontal rule for visual separation

# --- Rules Panel at the very bottom ---
st.markdown("## ⚖️ Rules for Combining ACMG/AMP Criteria")
st.markdown(
    """
    <style>
    .criteria-pathogenic {background: #ffe5e8; padding: 10px; border-radius: 8px;}
    .criteria-likely-pathogenic {background: #fff7e1; padding: 10px; border-radius: 8px;}
    .criteria-benign {background: #e1ffe1; padding: 10px; border-radius: 8px;}
    .criteria-likely-benign {background: #e6f0ff; padding: 10px; border-radius: 8px;}
    </style>
    """,
    unsafe_allow_html=True
)

with st.expander(":red_circle: Pathogenic"):
    st.markdown(
        """
        <div class="criteria-pathogenic">
        • 1 Very Strong (PVS1, PS2_Very Strong, PS4_Very Strong) **AND** ≥ 1 Strong (PS1, PS2, PS4, PP1_Strong)<br>
        • 1 Very Strong **AND** ≥ 2 Moderate (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)<br>
        • 1 Very Strong **AND** 1 Moderate **AND** 1 Supporting<br>
        • 1 Very Strong **AND** ≥ 2 Supporting<br>
        • ≥ 2 Strong<br>
        • 1 Strong **AND** ≥ 3 Moderate<br>
        • 1 Strong **AND** 2 Moderate **AND** ≥ 2 Supporting<br>
        • 1 Strong **AND** 1 Moderate **AND** ≥ 4 Supporting
        </div>
        """,
        unsafe_allow_html=True
    )

with st.expander(":orange_circle: Likely Pathogenic"):
    st.markdown(
        """
        <div class="criteria-likely-pathogenic">
        • 1 Very Strong **AND** 1 Moderate<br>
        • 1 Strong **AND** 1 Moderate<br>
        • 1 Strong **AND** ≥ 2 Supporting<br>
        • ≥ 3 Moderate<br>
        • 2 Moderate **AND** ≥ 2 Supporting<br>
        • 1 Moderate **AND** ≥ 4 Supporting<br>
        • 1 Strong **AND** 2 Moderate
        </div>
        """,
        unsafe_allow_html=True
    )

with st.expander(":green_circle: Benign"):
    st.markdown(
        """
        <div class="criteria-benign">
        • ≥ 2 Strong (BS1, BS2, BS4, BP2_Strong)<br>
        • 1 Stand Alone (BA1)
        </div>
        """,
        unsafe_allow_html=True
    )

with st.expander(":large_blue_circle: Likely Benign"):
    st.markdown(
        """
        <div class="criteria-likely-benign">
        • 1 Strong **AND** 1 Supporting<br>
        • ≥ 2 Supporting
        </div>
        """,
        unsafe_allow_html=True
    )