import streamlit as st

# DO NOT call st.set_page_config() in page files when using pages/ directory

st.title("Rules for Combining Criteria")

st.markdown(
    """
These rules summarize how ACMG/AMP evidence strengths are combined to derive final variant classifications.
They reflect the standard framework applied by the VHL VCEP when integrating PVS, PS, PM, PP, BS, BP, and BA codes. [web:1][web:13]
"""
)

st.markdown("## 🔴 Pathogenic (P)")
st.markdown(
    """
- **1 Very Strong** (PVS1, PS2_VeryStrong, PS4_VeryStrong) **AND ≥ 1 Strong** (PS1, PS2, PS4, PP1_Strong)  
- **1 Very Strong** (PVS1, PS2_VeryStrong, PS4_VeryStrong) **AND ≥ 2 Moderate**  
  - Moderate: (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)  
- **1 Very Strong** (PVS1, PS2_VeryStrong, PS4_VeryStrong) **AND 1 Moderate AND 1 Supporting**  
  - Moderate: (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)  
  - Supporting: (PS2_Supporting, PS3_Supporting, PS4_Supporting, PM1_Supporting, PM2_Supporting, PP1, PP3)  
- **1 Very Strong** (PVS1, PS2_VeryStrong, PS4_VeryStrong) **AND ≥ 2 Supporting**  
  - Supporting: (PS2_Supporting, PS3_Supporting, PS4_Supporting, PM1_Supporting, PM2_Supporting, PP1, PP3)  
- **≥ 2 Strong** (PS1, PS2, PS4, PP1_Strong)  
- **1 Strong** (PS1, PS2, PS4, PP1_Strong) **AND ≥ 3 Moderate**  
  - Moderate: (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)  
- **1 Strong** (PS1, PS2, PS4, PP1_Strong) **AND 2 Moderate AND ≥ 2 Supporting**  
  - Moderate: (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)  
  - Supporting: (PS2_Supporting, PS3_Supporting, PS4_Supporting, PM1_Supporting, PM2_Supporting, PP1, PP3)  
- **1 Strong** (PS1, PS2, PS4, PP1_Strong) **AND 1 Moderate AND ≥ 4 Supporting**  
  - Moderate: (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)  
  - Supporting: (PS2_Supporting, PS3_Supporting, PS4_Supporting, PM1_Supporting, PM2_Supporting, PP1, PP3)
"""
)

st.markdown("## 🟠 Likely Pathogenic (LP)")
st.markdown(
    """
- **1 Very Strong** (PVS1, PS2_VeryStrong, PS4_VeryStrong) **AND 1 Moderate**  
  - Moderate: (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)  
- **1 Strong** (PS1, PS2, PS4, PP1_Strong) **AND 1 Moderate**  
  - Moderate: (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)  
- **1 Strong** (PS1, PS2, PS4, PP1_Strong) **AND ≥ 2 Supporting**  
  - Supporting: (PS2_Supporting, PS3_Supporting, PS4_Supporting, PM1_Supporting, PM2_Supporting, PP1, PP3)  
- **≥ 3 Moderate**  
  - Moderate: (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)  
- **2 Moderate AND ≥ 2 Supporting**  
  - Moderate: (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)  
  - Supporting: (PS2_Supporting, PS3_Supporting, PS4_Supporting, PM1_Supporting, PM2_Supporting, PP1, PP3)  
- **1 Moderate AND ≥ 4 Supporting**  
  - Moderate: (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)  
  - Supporting: (PS2_Supporting, PS3_Supporting, PS4_Supporting, PM1_Supporting, PM2_Supporting, PP1, PP3)  
- **1 Strong** (PS1, PS2, PS4, PP1_Strong) **AND 2 Moderate**  
  - Moderate: (PS2_Moderate, PS4_Moderate, PM1, PM4, PM5, PM6, PP1_Moderate)
"""
)

st.markdown("## 🟢 Benign (B)")
st.markdown(
    """
- **≥ 2 Strong** (BS1, BS2, BS4, BP2_Strong)  
- **1 Stand Alone** (BA1)
"""
)

st.markdown("## 🟡 Likely Benign (LB)")
st.markdown(
    """
- **1 Strong** (BS1, BS2, BS4, BP2_Strong) **AND 1 Supporting**  
  - Supporting: (BS2_Supporting, BS3_Supporting, BS4_Supporting, BP2, BP3, BP4, BP5, BP7)  
- **≥ 2 Supporting**  
  - Supporting: (BS2_Supporting, BS3_Supporting, BS4_Supporting, BP2, BP3, BP4, BP5, BP7)
"""
)
