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


BASE_DIR = Path(__file__).parent.parent
data_path = BASE_DIR / "Findlay_MAVE_processed.csv"
df = load_mave_data(data_path)


# --- Choose the coordinate column that runs along the gene ---
# CHANGE THIS if your coordinate is named differently (e.g. "alt_pos", "cds_pos", etc.).
coord_col = "ref_pos" if "ref_pos" in df.columns else "alt_pos"
if coord_col not in df.columns:
    st.error(
        "No suitable coordinate column found (expected 'ref_pos' or 'alt_pos'). "
        "Please update coord_col to match your processed CSV."
    )
    st.stop()

df = df.sort_values(coord_col).reset_index(drop=True)
df["genomic_position"] = df[coord_col]
df["variant_index"] = range(len(df))


# --- ClinVar categorical ordering ---
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


# --- Download data ---
st.download_button(
    label="Download processed MAVE CSV",
    data=df.to_csv(index=False),
    file_name="Findlay_MAVE_processed.csv",
    mime="text/csv",
)


# --- Sidebar filters ---
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
if selected_clinvar:
    mask &= df["clinvar_simple"].isin(selected_clinvar)

df_plot = df[mask].copy()
if df_plot.empty:
    st.warning("No variants match the current filters.")
    st.stop()


# --- Infer exon spans from sge_region on the same coordinate axis ---
exon_bounds = pd.DataFrame(columns=["region", "min", "max"])
if "sge_region" in df_plot.columns:
    tmp = (
        df_plot.groupby("sge_region")["genomic_position"]
        .agg(["min", "max"])
        .reset_index()
        .rename(columns={"sge_region": "region"})
    )
    # Keep only regions that look like exons (tweak condition to match your labels)
    exon_bounds = tmp[tmp["region"].str.contains("exon", case=False, na=False)]

y_min = df_plot["function_score_final"].min()
y_max = df_plot["function_score_final"].max()


# --- Build scatter: x = genomic_position, y = function score ---
fig = px.scatter(
    df_plot,
    x="genomic_position",
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
        "genomic_position": True,
        "variant_index": True,
    },
    color_discrete_sequence=px.colors.qualitative.Set2,
)


# --- Add exon bands based on inferred bounds ---
exon_colors = px.colors.qualitative.Pastel
for i, row in exon_bounds.reset_index(drop=True).iterrows():
    color = exon_colors[i % len(exon_colors)]

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

# Intronic variants automatically appear between these bands because they share the same x-axis.


# --- Layout: full-width, clean look ---
fig.update_yaxes(range=[y_min - 0.2, y_max + 0.6])

fig.update_layout(
    template="plotly_dark",
    xaxis_title=f"{coord_col} (position along VHL)",
    yaxis_title="Function score",
    legend_title="ClinVar annotation",
    height=650,
    margin=dict(l=40, r=40, t=90, b=40),
)

st.plotly_chart(fig, use_container_width=True)
