import streamlit as st
import re
import pandas as pd
from vhl_pvs1 import classify_vhl_pvs1
from vhl_ps1 import classify_vhl_ps1
from vhl_ps2 import classify_vhl_ps2

# Configure app ONCE, in main entry file
st.set_page_config(page_title="VHL/VCEP Classifier", page_icon="🧬", layout="wide")

# ---------- Sidebar (About only) ----------
st.sidebar.markdown(
    """
### About VHL VCEP

The **VHL Variant Curation Expert Panel (VCEP)** is part of ClinGen’s effort to provide expert-level clinical validity for variants in the *VHL* gene. The committee is chaired by **Dr. Raymond Kim**, whose work at Princess Margaret Cancer Centre, SickKids Hospital, and the Early Cancer Detection Program focuses on hereditary cancer genetics and variant interpretation. The panel has developed a disease-specific annotation protocol using Hypothes.is to enable community curation of VHL variants and accelerate resolution of variants of uncertain significance.

📄 Read more [here:](Developing a disease-specific annotation protocol for VHL gene curation using Hypothes.is](https://pmc.ncbi.nlm.nih.gov/articles/PMC9825735/)
"""
)

# ---------- Main: Classifier ----------
st.title("VHL/VCEP Classifier (v1.2 TEST) 🧬")

st.markdown(
    """
Enter a variant in HGVS (cDNA) notation to predict **PVS1**, **PS1**, and **PS2** strength according to VHL VCEP guidelines.  
Classifiers for additional criteria are under development. Context explains why each strength label was applied.  
Refer to the original criteria [here](https://cspec.genome.network/cspec/ui/svi/doc/GN078).

**Note:** PM3, PP2, PP4, PP5, BP1, BP6 are *not recommended* for use by the ClinGen Sequence Variant Interpretation VCEP Review Committee in the case of VHL.
"""
)

hgvs_input = st.text_input("HGVS c. Notation (see examples below)", "")

with st.expander("PVS1/PS1 Options"):
    st.markdown("### PVS1 Options")
    exon_skipping = st.checkbox("Exon Skipping", value=False)
    cryptic_disrupts_rf = st.radio(
        "Cryptic Splice Disrupts Reading Frame?", (None, True, False)
    )
    cryptic_preserves_rf = st.radio(
        "Cryptic Splice Preserves Reading Frame?", (None, True, False)
    )
    duplication_type = st.selectbox(
        "Duplication Type", (None, "tandem", "not_in_tandem")
    )

with st.expander("PS2 (De Novo) Options"):
    st.markdown(
        """
**Scoring Clarity for PS2**

- If **any family history** of VHL disease is present, PS2 evidence cannot be assigned and results will indicate this restriction.
"""
    )
    family_history = st.checkbox(
        "Any family history of VHL?",
        value=False,
        help="Check if any features of VHL or known pathogenic variant(s) present in family. PS2 cannot be scored in such cases.",
    )

    if family_history:
        st.warning(
            "PS2 cannot be assigned if any family history of VHL disease is present. This overrides all other options."
        )
        is_de_novo = False
        phenotype = None
        panel_neg = {}
    else:
        st.markdown(
            """
To achieve **Moderate (1 point)** PS2 evidence:

1. Confirmed de novo status (both maternity and paternity testing done)  
2. For *Consistent* or *RCC+Pheo* phenotypes, all required genes must be negative.

If panel is incomplete or untested, only **Supporting (0.5 points)** evidence can be assigned, even with confirmed de novo status.
"""
        )
        is_de_novo = st.checkbox(
            "Confirmed de novo (both maternity and paternity tested)",
            value=False,
            help="Check only if BOTH maternity and paternity testing confirm de novo status.",
        )

        st.markdown("**Select phenotype category:**")
        phenotype = st.selectbox(
            "Phenotype category",
            (None, "danish", "consistent", "nonspecific"),
            format_func=lambda x: {
                None: "Not specified",
                "danish": "Highly specific (Danish/International criteria)",
                "consistent": "Consistent with gene (e.g. VHL2c, RCC+Pheo as per specs)",
                "nonspecific": "Nonspecific phenotype",
            }.get(x, x),
        )

        st.markdown(
            """
**Panel Results (should be negative to maximize evidence):**

- For VHL2c (pheo only): all [MAX, NF1, RET, SDHA, SDHB, SDHC, SDHD, SDHAF2, TMEM127, VHL] should be marked negative.  
- For RCC+Pheo: all [MAX, FH, SDHA, SDHB, SDHC, SDHD, SDHAF2, TMEM127] should be negative.

If incomplete, only PS2_Supporting can be assigned.
"""
        )

        st.markdown("**Mark negative testing results below only if confirmed:**")
        panel_genes = [
            "SDHB",
            "SDHC",
            "SDHD",
            "RET",
            "MAX",
            "FH",
            "TMEM127",
            "NF1",
            "SDHA",
            "SDHAF2",
            "VHL",
        ]
        panel_neg = {}
        cols = st.columns(4)
        for i, gene in enumerate(panel_genes):
            with cols[i % 4]:
                panel_neg[gene] = (
                    "neg"
                    if st.checkbox(f"{gene} negative", value=False, key=gene)
                    else None
                )

# ---------- Example variants table ----------
st.markdown("### Example variants")

example_variants = pd.DataFrame(
    {
        "Example": [
            "NM_000551.4(VHL):c.263+1G>A",
            "NM_000551.4(VHL):c.119_340del",
            "NM_000551.4(VHL):c.123+5G>A",
            "NM_000551.4(VHL):c.1A>T (p.Met1Leu)",
            "NM_000551.4(VHL):c.160A>T (p.Met54Leu)",
            "NM_000551.4(VHL):c.150G>A (p.Pro50=)",
            "NM_000551.4(VHL):c.408del (p.Phe136fs)",
            "NM_000551.4(VHL):c.120_122del (p.Val41del)",
            "NM_000551.4(VHL):c.190_195dup (p.Lys64_Leu65dup)",
            "NM_000551.4(VHL):c.263G>A (p.Trp88Ter)",
        ],
        "Context": [
            "Canonical donor splice variant to exercise PVS1 splice rules.",
            "Multi-exon deletion to demonstrate very strong truncating PVS1 evidence.",
            "Non-canonical splice variant to test nuanced PVS1 splice strength.",
            "Start-loss variant at the initiation codon for special PVS1 handling.",
            "Missense in a functionally important region to explore PS1/PM1 logic.",
            "Synonymous variant with no amino acid change for benign/BP-style scenarios.",
            "Single-base frameshift to test NMD boundary and truncating logic.",
            "In-frame deletion for PM4 and domain-specific protein effects.",
            "In-frame duplication to check tandem vs non-tandem duplication rules.",
            "Classic nonsense variant to illustrate straightforward very strong PVS1.",
        ],
    }
)

def blue_table_style(df):
    return df.style.set_properties(
        **{
            "background-color": "#ffffff",
            "border-color": "#a3b8e6",
        }
    )

st.dataframe(
    blue_table_style(example_variants),
    use_container_width=True,
)

# ---------- Run classifiers only if variant entered ----------
if hgvs_input:
    pvs1_result = classify_vhl_pvs1(
        hgvs_input,
        exon_skipping=exon_skipping,
        cryptic_disrupts_rf=True if cryptic_disrupts_rf is True else None,
        cryptic_preserves_rf=True if cryptic_preserves_rf is True else None,
        duplication_type=duplication_type
        if duplication_type not in (None, "None")
        else None,
    )
    ps1_result = classify_vhl_ps1(hgvs_input)
    ps2_result = classify_vhl_ps2(
        hgvs_input,
        is_de_novo=is_de_novo if not family_history else False,
        phenotype=phenotype if not family_history and phenotype is not None else None,
        panel_neg={
            gene: result
            for gene, result in panel_neg.items()
            if result == "neg"
        }
        if not family_history
        else {},
        family_history=family_history,
    )

    st.header("Classification Results")
    st.subheader("PVS1 Result")
    st.write(pvs1_result)
    st.subheader("PS1 Result")
    st.write(ps1_result)
    st.subheader("PS2 Result")
    st.write(ps2_result)