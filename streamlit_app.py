import streamlit as st
import re
from vhl_pvs1 import classify_vhl_pvs1  # Import PVS1 classifier

st.title("VHL Variant Classifiers")

st.markdown(
    """
    Enter a variant in HGVS (cDNA) notation to predict PVS1 and PS1 strength and context for VHL.
    """
)

hgvs_input = st.text_input(
    "HGVS c. Notation (e.g., c.200C>T, c.100_200del, c.123+1G>A):",
    ""
)

with st.expander("Advanced options"):
    exon_skipping = st.checkbox("Exon Skipping", value=False)
    cryptic_disrupts_rf = st.radio("Cryptic Splice Disrupts Reading Frame?", (None, True, False))
    cryptic_preserves_rf = st.radio("Cryptic Splice Preserves Reading Frame?", (None, True, False))
    duplication_type = st.selectbox("Duplication Type", (None, "tandem", "not_in_tandem"))
    initiation_codon = st.selectbox("Initiation Codon", (None, "Met1", "Met54"))

if hgvs_input:
    pvs1_result = classify_vhl_pvs1(
        hgvs_input,
        exon_skipping=exon_skipping,
        cryptic_disrupts_rf=True if cryptic_disrupts_rf is True else None,
        cryptic_preserves_rf=True if cryptic_preserves_rf is True else None,
        duplication_type=duplication_type if duplication_type not in (None, "None") else None,
        initiation_codon=initiation_codon if initiation_codon not in (None, "None") else None
    )

    st.header("Classification Results")
    st.subheader("PVS1 Result")
    st.write(pvs1_result)