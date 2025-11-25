import streamlit as st
import re
from vhl_pvs1 import classify_vhl_pvs1  # Import PVS1 classifier
from vhl_ps1 import classify_vhl_ps1    # Import PS1 classifier

st.title("VHL/VCEP Classifer (v1.0 TEST) 🧬")

st.markdown(
    """
    Enter a variant in HGVS (cDNA) notation to predict PVS1 and PS1 strength in accordance with VHL VCEP guidelines (Classifers for remaining criteria will be uploading in the coming days! This is a test!).
    Please note strength returns "null" if the criteria is not relevant. Context is meant to orient you as to why the strength label was applied. For your reference, the original
    criteria are [here](https://cspec.genome.network/cspec/ui/svi/doc/GN078#pmid_35709961).

    Some suggested variants to try entering:

    - NM_000551.4(VHL):c.263+1G>A
    - NM_000551.4(VHL):c.119_340del
    - NM_000551.4(VHL):c.123+5G>A
    - NM_000551.4(VHL):c.1A>T (p.Met1Leu)
    - NM_000551.4(VHL):c.160A>T (p.Met54Leu)
    - NM_000551.4(VHL):c.150G>A (p.Pro50=)
    - NM_000551.4(VHL):c.408del (p.Phe136fs)
    - NM_000551.4(VHL):c.120_122del (p.Val41del)
    - NM_000551.4(VHL):c.190_195dup (p.Lys64_Leu65dup)
    - NM_000551.4(VHL):c.263G>A (p.Trp88Ter)
    """
)

hgvs_input = st.text_input(
    "use HGVS c. Notation in alignment with the above examples",
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
    ps1_result = classify_vhl_ps1(hgvs_input)  # Adjust arguments if needed for your PS1 classifier

    st.header("Classification Results")
    st.subheader("PVS1 Result")
    st.write(pvs1_result)
    st.subheader("PS1 Result")
    st.write(ps1_result)