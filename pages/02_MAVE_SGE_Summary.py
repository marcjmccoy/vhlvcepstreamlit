import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path

st.title("MAVE / SGE summary")


st.markdown(
    """
This page summarizes [VHL saturation genome editing (SGE/MAVE) data from Findlay et al.](https://www.nature.com/articles/s41588-024-01800-z)
showing functional scores across the gene, colored by ClinVar-level annotation.
"""
)


# --- Load data ---
@st.cache_data
def load_mave_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Clean/standardize ClinVar annotation labels a bit
    df["clinvar_simple"] = df["clinvar_simple"].str.strip()
    return df


data_path = "data/Findlay_MAVE_processed.csv"  # update if needed
df = load_mave_data(data_path)


# --- Download data ---
st.download_button(
    label="Download processed MAVE CSV",
    data=df.to_csv(index=False),
    file_name="Findlay_MAVE_processed.csv",
    mime="text/csv",
)


# --- Basic preprocessing for plotting ---
# Use alt_pos (or ref_pos) as a genomic-position proxy along x-axis
df["genomic_position"] = df["alt_pos"]


# Map ClinVar categories to a reduced set if desired
clinvar_order = [
    "Absent",
    "Benign",
    "Likely benign",
    "Uncertain significance",
    "Conflicting interpretations of pathogenicity",
    "Likely pathogenic",
    "Pathogenic",
]
df["clinvar_simple"] = pd.Categorical(df["clinvar_simple"], categories=clinvar_order, ordered=True)


# --- Sidebar controls ---
with st.sidebar:
    st.header("MAVE / SGE filters")
    # Filter by consequence type
    consequence_options = sorted(df["consequence"].dropna().unique())
    selected_consequences = st.multiselect(
        "Consequence",
        options=consequence_options,
        default=consequence_options,
    )

    # Filter by SGE region / exon region if desired
    region_options = sorted(df["sge_region"].dropna().unique())
    selected_regions = st.multiselect(
        "SGE region",
        options=region_options,
        default=region_options,
    )

    # Option to restrict by ClinVar category
    clinvar_options = [c for c in clinvar_order if c in df["clinvar_simple"].unique()]
    selected_clinvar = st.multiselect(
        "ClinVar annotation",
        options=clinvar_options,
        default=clinvar_options,
    )

    # Toggle to use scaled vs raw genomic position (here simply rescaled 0–1)
    scale_genomic = st.checkbox("Genomic position at scale", value=False)


# --- Apply filters ---
mask = (
    df["consequence"].isin(selected_consequences)
    & df["sge_region"].isin(selected_regions)
    & df["clinvar_simple"].isin(selected_clinvar)
)
df_plot = df[mask].copy()


if scale_genomic:
    # Simple min-max scaling for display
    x_min, x_max = df_plot["genomic_position"].min(), df_plot["genomic_position"].max()
    df_plot["x_position"] = (df_plot["genomic_position"] - x_min) / (x_max - x_min)
    xaxis_title = "Scaled genomic position"
else:
    df_plot["x_position"] = df_plot["genomic_position"]
    xaxis_title = "Unscaled genomic position"


# --- Build scatter plot (similar layout to the screenshot) ---
fig = px.scatter(
    df_plot,
    x="x_position",
    y="function_score_final",
    color="clinvar_simple",
    hover_name="Preferred Variant Title",
    hover_data={
        "sge_region": True,
        "consequence": True,
        "tier_class": True,
        "function_score_final": True,
        "rna_score": True,
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
