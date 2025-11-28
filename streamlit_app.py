import streamlit as st
import pandas as pd

from vhl_pvs1 import classify_vhl_pvs1
from vhl_ps1 import classify_vhl_ps1
from vhl_ps2 import classify_vhl_ps2
from vhl_pm4 import classify_vhl_pm4


# ---------------- Page config ----------------
st.set_page_config(
    page_title="VHL/VCEP Classifier",
    page_icon="🧬",
    layout="wide",
)

# ---------------- CSS for wrapped tables ----------------
st.markdown(
    """
    <style>
    .vhl-table table {
        width: 100%;
        border-collapse: collapse;
    }
    .vhl-table th, .vhl-table td {
        white-space: normal !important;
        word-wrap: break-word;
        max-width: 350px;
        text-align: left;
        vertical-align: top;
        padding: 0.25rem 0.5rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


def show_wrapped_table(df: pd.DataFrame) -> None:
    """Render a DataFrame with wrapped text in cells."""
    html = df.to_html(index=False, escape=False)
    st.markdown(f'<div class="vhl-table">{html}</div>', unsafe_allow_html=True)


# ---------------- Sidebar: About ----------------
st.sidebar.markdown(
    """
### About VHL VCEP

The **VHL Variant Curation Expert Panel (VCEP)** is part of ClinGen’s effort to provide
expert-level clinical validity for variants in the *VHL* gene. The committee is chaired
by **Dr. Raymond Kim**, whose work at Princess Margaret Cancer Centre and the Early Cancer
Detection Program focuses on hereditary cancer genetics and variant interpretation. His team
has developed a disease-specific annotation protocol using Hypothes.is to enable community
curation of VHL variants and accelerate resolution of variants of uncertain significance.

- 📄 Read more in the [VHL VCEP protocol paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9825735/)
- 📄 Reference the [ClinGen VHL Expert Panel Guideline](https://cspec.genome.network/cspec/ui/svi/doc/GN078)
"""
)


# ---------------- Helper: Example variants table ----------------
def build_example_variants_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Variant": [
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
                "NM_000551.4(VHL):c.189_191del",
            ],
            "Rationale": [
                "Canonical +1 donor splice change expected to abolish splicing in a loss-of-function VHL gene.",
                "Deletion spanning key coding regions across multiple exons in a 3-exon gene.",
                "+5 intronic variant near a donor site with uncertain splice impact.",
                "Loss of the canonical initiation codon that may or may not be rescued by alternative starts.",
                "Missense change at a residue within or near a functionally important VHL region.",
                "Synonymous variant with no amino-acid change.",
                "Single-base deletion causing a frameshift in the distal coding region.",
                "Three-base in-frame deletion removing one conserved residue.",
                "In-frame tandem duplication adding two residues without a frameshift.",
                "Nonsense variant introducing a premature stop in the central part of VHL.",
                "In-frame deletion that removes codon 63, which lies within the critical beta domain (AA 63–155) of pVHL"
            ],
        }
    )


# ---------------- Main title and intro ----------------
st.title("VHL/VCEP Classifier (v2.0 TEST) 🧬")

st.markdown(
    """
Enter a variant in HGVS (cDNA) notation to predict **PVS1**, **PS1**, and **PS2** strength
according to VHL VCEP guidelines. Classifiers for additional criteria are under development.
The context fields in the output explain why each strength label was applied. Refer to the
original [ACMG/AMP Variant Interpretation Guidelines for VHL Version 1.1.0](https://cspec.genome.network/cspec/ui/svi/doc/GN078)

**Note:** PM3, PP2, PP4, PP5, BP1, BP6 are *not recommended* for use by the ClinGen SVI VCEP Review Committee for VHL.
"""
)

# ---------------- Example variants (above input) ----------------
st.markdown("### Example variants")
example_variants = build_example_variants_df()
show_wrapped_table(example_variants)

# ---------------- HGVS input ----------------
hgvs_input = st.text_input("HGVS c. Notation (see examples above)", "")

# ---------------- PVS1 / PS1 options ----------------
with st.expander("PVS1 / PS1 Options"):
    st.markdown("#### PVS1 Options")

    exon_skipping = st.checkbox("Exon Skipping", value=False)

    cryptic_disrupts_rf = st.radio(
        "Cryptic splice site disrupts reading frame?",
        (None, True, False),
        index=0,
    )

    cryptic_preserves_rf = st.radio(
        "Cryptic splice site preserves reading frame?",
        (None, True, False),
        index=0,
    )

    duplication_type = st.selectbox(
        "Duplication type",
        (None, "tandem", "not_in_tandem"),
        index=0,
    )

# ---------------- PS2 (De novo) options ----------------
with st.expander("PS2 (De Novo) Options"):
    st.markdown(
        """
**Scoring clarity for PS2 (de novo):**

- If **any family history** of VHL disease is present, PS2 evidence cannot be assigned.
- When there is no family history, PS2 strength depends on:
  - Confirmed de novo status (maternity and paternity tested).
  - Phenotype category.
  - Completeness and negativity of the relevant gene panel.
"""
    )

    family_history = st.checkbox(
        "Any family history of VHL?",
        value=False,
        help=(
            "Check if any features of VHL or known pathogenic variant(s) are present "
            "in the family. PS2 cannot be scored in such cases."
        ),
    )

    if family_history:
        st.warning(
            "PS2 cannot be assigned if any family history of VHL disease is present. "
            "This overrides all other options."
        )
        is_de_novo = False
        phenotype = None
        panel_neg = {}
    else:
        st.markdown(
            """
To achieve **Moderate (1 point)** PS2 evidence:

1. Confirmed de novo status (both maternity and paternity testing performed).
2. For *Highly specific* or *Consistent* phenotypes, all required genes must be negative
   on panel testing.

If the panel is incomplete or not performed, only **Supporting (0.5 points)** PS2
evidence can be assigned, even with confirmed de novo status.
"""
        )

        is_de_novo = st.checkbox(
            "Confirmed de novo (both maternity and paternity tested)",
            value=False,
            help=(
                "Check only if BOTH maternity and paternity testing confirm that the "
                "variant is de novo."
            ),
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
**Panel results (should be negative to maximize evidence):**

- For VHL2c (pheochromocytoma-only): all of
  [MAX, NF1, RET, SDHA, SDHB, SDHC, SDHD, SDHAF2, TMEM127, VHL] should be negative.
- For RCC+Pheo: all of
  [MAX, FH, SDHA, SDHB, SDHC, SDHD, SDHAF2, TMEM127] should be negative.

If these panels are incomplete or unknown, only PS2_Supporting can be assigned.
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
                checked = st.checkbox(f"{gene} negative", value=False, key=gene)
                panel_neg[gene] = "neg" if checked else None


# ---------------- Run classifiers ----------------
def run_classifiers(
    hgvs: str,
    exon_skipping: bool,
    cryptic_disrupts_rf,
    cryptic_preserves_rf,
    duplication_type,
    is_de_novo: bool,
    phenotype,
    panel_neg: dict,
    family_history: bool,
):
    # PVS1
    pvs1_result = classify_vhl_pvs1(
        hgvs,
        exon_skipping=exon_skipping,
        cryptic_disrupts_rf=True if cryptic_disrupts_rf is True else None,
        cryptic_preserves_rf=True if cryptic_preserves_rf is True else None,
        duplication_type=duplication_type
        if duplication_type not in (None, "None")
        else None,
    )

    # PS1
    ps1_result = classify_vhl_ps1(hgvs)

    # PS2
    effective_is_de_novo = is_de_novo if not family_history else False
    effective_phenotype = (
        phenotype if (not family_history and phenotype is not None) else None
    )
    effective_panel = (
        {gene: result for gene, result in panel_neg.items() if result == "neg"}
        if not family_history
        else {}
    )

    ps2_result = classify_vhl_ps2(
        hgvs,
        is_de_novo=effective_is_de_novo,
        phenotype=effective_phenotype,
        panel_neg=effective_panel,
        family_history=family_history,
    )

    # PM4
    pm4_result = classify_vhl_pm4(
        hgvs
    )

    return pvs1_result, ps1_result, ps2_result, pm4_result

# ---------------- Trigger classification ----------------
if hgvs_input:
    pvs1_result, ps1_result, ps2_result, pm4_result = run_classifiers(
        hgvs_input,
        exon_skipping=exon_skipping,
        cryptic_disrupts_rf=cryptic_disrupts_rf,
        cryptic_preserves_rf=cryptic_preserves_rf,
        duplication_type=duplication_type,
        is_de_novo=is_de_novo,
        phenotype=phenotype,
        panel_neg=panel_neg,
        family_history=family_history,
    )

    st.header("Classification results")

    st.subheader("PVS1")
    st.write(pvs1_result)

    st.subheader("PS1")
    st.write(ps1_result)

    st.subheader("PS2")
st.write(ps2_result)

    st.subheader("PM4")
    st.write(pm4_result)