import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path

st.title("MAVE / SGE summary")

st.markdown(
    """
This page summarizes [VHL saturation genome editing (SGE/MAVE) data from Findlay et al.](https://www.nature.com/articles/s41588-024-01800-z),
showing functional scores across the gene, colored by ClinVar-level annotation.
"""
)


# --- Load data ---
@st.cache_data
def load_mave_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["clinvar_simple"] = df["clinvar_simple"].astype(str).str.strip()
    return df


# CSV is in project root, one level above `pages/`
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
if "sge_region" not in df.columns:
    st.error("Column 'sge_region' is missing from the CSV. Needed for exon annotation.")
    st.stop()

# Sort by genomic position and create dense variant index
df = df.sort_values("alt_pos").reset_index(drop=True)
df["genomic_position"] = df["alt_pos"]
df["variant_index"] = range(len(df))


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
    region_options = sorted(df["sge_region"].dropna().unique())
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


# --- Build base scatter plot (x-axis = variant index) ---
fig = px.scatter(
    df_plot,
    x="variant_index",
    y="function_score_final",
    color="clinvar_simple",
    hover_name=(
        "Preferred Variant Title"
        if "Preferred Variant Title" in df_plot.columns
        else None
    ),
    hover_data={
        "sge_region": True,
        "consequence": "consequence" in df_plot.columns,
        "tier_class": "tier_class" in df_plot.columns,
        "function_score_final": True,
        "rna_score": "rna_score" in df_plot.columns,
        "genomic_position": True,
        "variant_index": True,
    },
    color_discrete_sequence=px.colors.qualitative.Set2,
)

# --- Add exon / SGE-region bands on x-axis ---
region_bounds = (
    df_plot.groupby("sge_region")["variant_index"]
    .agg(["min", "max"])
    .reset_index()
)

exon_colors = px.colors.qualitative.Pastel
for i, row in region_bounds.iterrows():
    color = exon_colors[i % len(exon_colors)]
    # translucent vertical rectangle
    fig.add_vrect(
        x0=row["min"],
        x1=row["max"],
        fillcolor=color,
        opacity=0.18,
        line_width=0,
        layer="below",
    )
    # exon label at top of band
    fig.add_annotation(
        x=(row["min"] + row["max"]) / 2,
        y=df_plot["function_score_final"].max() + 0.2,
        text=row["sge_region"],
        showarrow=False,
        font=dict(size=10, color="#CCCCCC"),
        yanchor="bottom",
    )


# --- Layout tweaks for a cleaner, wider look ---
fig.update_layout(
    template="plotly_dark",
    xaxis_title="Variant index",
    yaxis_title="Function score",
    legend_title="ClinVar annotation",
    height=600,
    margin=dict(l=40, r=40, t=80, b=40),
)

# Slight extra headroom for exon labels
fig.update_yaxes(range=[
    df_plot["function_score_final"].min() - 0.2,
    df_plot["function_score_final"].max() + 0.6,
])

st.plotly_chart(fig, use_container_width=True)
