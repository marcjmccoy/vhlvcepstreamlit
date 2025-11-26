import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path


st.title("MAVE / SGE summary")

st.markdown(
    """
This page summarizes [VHL saturation genome editing (SGE/MAVE) data from Findlay et al.](https://www.nature.com/articles/s41588-024-01800-z),
showing functional scores across the gene, with exon / intron structure annotated along the x-axis.
"""
)


# --- Load data ---
@st.cache_data
def load_mave_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["clinvar_simple"] = df["clinvar_simple"].astype(str).str.strip()
    return df


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
if "sge_region" not in df.columns:
    st.error("Column 'sge_region' is required to define exon / intron regions.")
    st.stop()

# Order regions in a biologically sensible way: exons then introns.
# Adjust this list to match your exact labels.
region_order = [
    "exon 1a", "exon 1b", "intron 1", "exon 2",
    "intron 2", "exon 3a", "exon 3b"
]
# Fallback: use whatever unique values exist if any are missing.
present_regions = [r for r in region_order if r in df["sge_region"].unique()]
if not present_regions:
    present_regions = list(df["sge_region"].dropna().unique())

df["sge_region"] = pd.Categorical(df["sge_region"], categories=present_regions, ordered=True)

# Sort by region, then by any positional proxy if available
pos_cols = [c for c in ["ref_pos", "cds_pos", "alt_pos"] if c in df.columns]
sort_cols = ["sge_region"] + pos_cols
df = df.sort_values(sort_cols).reset_index(drop=True)

# Dense positional index along gene; this is the x-axis
df["x_index"] = range(len(df))


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


# --- Sidebar filters ---
with st.sidebar:
    st.header("MAVE / SGE filters")

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
if selected_clinvar:
    mask &= df["clinvar_simple"].isin(selected_clinvar)

df_plot = df[mask].copy()
if df_plot.empty:
    st.warning("No variants match the current filters.")
    st.stop()


# --- Infer exon bands from x_index spans of each region ---
region_bounds = (
    df_plot.groupby("sge_region")["x_index"]
    .agg(["min", "max"])
    .reset_index()
    .rename(columns={"sge_region": "region"})
)

y_min = df_plot["function_score_final"].min()
y_max = df_plot["function_score_final"].max()


# --- Plot: x = x_index, color by sge_region, ClinVar in hover ---
fig = px.scatter(
    df_plot,
    x="x_index",
    y="function_score_final",
    color="sge_region",
    hover_name=(
        "Preferred Variant Title"
        if "Preferred Variant Title" in df_plot.columns
        else None
    ),
    hover_data={
        "sge_region": True,
        "clinvar_simple": True,
        "consequence": "consequence" in df_plot.columns,
        "tier_class": "tier_class" in df_plot.columns,
        "function_score_final": True,
        "rna_score": "rna_score" in df_plot.columns,
    },
    color_discrete_sequence=px.colors.qualitative.Set3,
)

# Add translucent bands for exons only (regions whose name contains 'exon')
band_colors = px.colors.qualitative.Pastel
for i, row in region_bounds.iterrows():
    if "exon" not in str(row["region"]).lower():
        continue
    color = band_colors[i % len(band_colors)]

    fig.add_vrect(
        x0=row["min"],
        x1=row["max"],
        fillcolor=color,
        opacity=0.18,
        line_width=0,
        layer="below",
    )

    fig.add_annotation(
        x=(row["min"] + row["max"]) / 2,
        y=y_max + 0.2,
        text=row["region"],
        showarrow=False,
        font=dict(size=10, color="#CCCCCC"),
        yanchor="bottom",
    )

fig.update_yaxes(range=[y_min - 0.2, y_max + 0.6])

fig.update_layout(
    template="plotly_dark",
    xaxis_title="Variant index (ordered by SGE region)",
    yaxis_title="Function score",
    legend_title="SGE / exon region",
    height=650,
    margin=dict(l=40, r=40, t=90, b=40),
)

st.plotly_chart(fig, use_container_width=True)
