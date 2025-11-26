import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path


st.title("MAVE / SGE summary")

st.markdown(
    """
This page summarizes VHL saturation genome editing (SGE/MAVE) data from Findlay et al.,
showing functional scores across the gene, colored by ClinVar-level annotation.
"""
)


# --- Load data ---
@st.cache_data
def load_mave_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["clinvar_simple"] = df["clinvar_simple"].astype(str).str.strip()
    return df


# CSV is in the project root, one level above `pages/`
BASE_DIR = Path(__file__).parent.parent
data_path = BASE_DIR / "Findlay_MAVE_processed.csv"
df = load_mave_data(data_path)


# --- Download data ---
st.download_button(
    label="Download processed MAVE CSV",
    data=df.to_csv(index=False),
    file_name="Findlay_MAVE_processed.csv",
    mime="text/csv",
)


# --- Basic preprocessing for plotting ---
if "alt_pos" not in df.columns:
    st.error("Column 'alt_pos' is missing from the CSV. Please check the input file.")
    st.stop()

df["genomic_position"] = df["alt_pos"]


# --- ClinVar category handling ---
clinvar_order = [
    "Absent",
    "Benign",
    "Likely benign",
    "Uncertain significance",
    "Conflicting interpretations of pathogenicity",
    "Likely pathogenic",
    "Pathogenic",
]
df["clinvar_simple"] = pd.Categorical(
    df["clinvar_simple"], categories=clinvar_order, ordered=True
)


# --- Sidebar controls ---
with st.sidebar:
    st.header("MAVE / SGE filters")

    # Consequence
    consequence_options = (
        sorted(df["consequence"].dropna().unique())
        if "consequence" in df.columns
        else []
    )
    selected_consequences = st.multiselect(
        "Consequence",
        options=consequence_options,
        default=consequence_options,
    )

    # SGE region
    region_options = (
        sorted(df["sge_region"].dropna().unique())
        if "sge_region" in df.columns
        else []
    )
    selected_regions = st.multiselect(
        "SGE region",
        options=region_options,
        default=region_options,
    )

    # ClinVar
    clinvar_options = [c for c in clinvar_order if c in df["clinvar_simple"].unique()]
    selected_clinvar = st.multiselect(
        "ClinVar annotation",
        options=clinvar_options,
        default=clinvar_options,
    )

    # Scaling toggle
    scale_genomic = st.checkbox("Genomic position at scale", value=False)


# --- Apply filters ---
mask = pd.Series(True, index=df.index)

if selected_consequences:
    mask &= df["consequence"].isin(selected_consequences)
if selected_regions:
    mask &= df["sge_region"].isin(selected_regions)
if selected_clinvar:
    mask &= df["clinvar_simple"].isin(selected_clinvar)

df_plot = df[mask].copy()

if df_plot.empty:
    st.warning("No variants match the current filters.")
    st.stop()


# --- X-axis scaling ---
if scale_genomic:
    x_min, x_max = df_plot["genomic_position"].min(), df_plot["genomic_position"].max()
    if x_max == x_min:
        df_plot["x_position"] = 0.5
    else:
        df_plot["x_position"] = (df_plot["genomic_position"] - x_min) / (x_max - x_min)
    xaxis_title = "Scaled genomic position"
else:
    df_plot["x_position"] = df_plot["genomic_position"]
    xaxis_title = "Unscaled genomic position"


# --- Plot ---
fig = px.scatter(
    df_plot,
    x="x_position",
    y="function_score_final",
    color="clinvar_simple",
    hover_name=(
        "Preferred Variant Title"
        if "Preferred Variant Title" in df_plot.columns
        else None
    ),
    hover_data={
        "sge_region": "sge_region" in df_plot.columns,
        "consequence": "consequence" in df_plot.columns,
        "tier_class": "tier_class" in df_plot.columns,
        "function_score_final": True,
        "rna_score": "rna_score" in df_plot.columns,
        "x_position": False,
        "genomic_position": True,
    },
    color_discrete_sequence=px.colors.qualitative.Set2,
)

fig.update_layout(
    template="plotly_dark",
    xaxis_title=xaxis_title,
    yaxis_title="Function score",
    legend_title="ClinVar annotation",
    height=500,
    margin=dict(l=40, r=40, t=60, b=40),
)

st.plotly_chart(fig, use_container_width=True)
