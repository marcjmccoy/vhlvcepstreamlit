# VHL/VCEP Classifier

VHL/VCEP Classifier is a Streamlit application that builds on the work of Dr. Raymond Kim and the ClinGen VHL Variant Curation Expert Panel (VCEP) to provide gene-wide, VHL-specific variant classification using ACMG/AMP criteria tailored for Von Hippel–Lindau disease and informed by multiplexed functional data. The app integrates VHL VCEP criteria specifications with saturation genome editing / MAVE-derived functional scores for nearly all single-nucleotide variants in VHL to support consistent, transparent variant interpretation across the gene.

## Overview

- Implements the VHL VCEP’s gene- and disease-specific ACMG/AMP specifications for VHL, including modified evidence codes (for example, PVS1, PS3, PS4, PM1, BS2, BS3, BS4, BP5) and unused codes, as published and versioned via ClinGen and related VHL curation manuscripts.
- Leverages high-throughput functional data from saturation genome editing / MAVE studies of VHL to transform variant effect scores into structured evidence that can be mapped to PS3/BS3 and related functional evidence categories.
- Aims to generate consistent provisional classifications across the entire coding region of VHL, while exposing the underlying evidence logic so expert users can review and override automated calls where appropriate.

## Scientific Background

Von Hippel–Lindau disease is a rare autosomal dominant cancer predisposition syndrome caused by pathogenic variants in VHL, with substantial phenotypic and allelic heterogeneity that complicates variant interpretation. The ClinGen VHL VCEP, led in part by Dr. Raymond Kim, has developed gene-specific ACMG/AMP specifications and a structured literature curation program to standardize VHL variant classification.

Recent saturation genome editing work has produced dense functional maps for VHL, measuring the impact of thousands of single-nucleotide variants on HIF-dependent cellular fitness and defining loss-of-function alleles with high accuracy. These VHL MAVE datasets enable systematic functional calibration of variants, providing quantitative scores that can be aligned with ACMG/AMP functional evidence strength levels and integrated alongside clinical, population, and computational data.

## Key Features

### VHL VCEP-aware ACMG/AMP engine

- Encodes the current VHL-specific ACMG/AMP criteria as published by the VCEP and exposed through the ClinGen Criteria Specifications Registry.
- Supports phenotype-driven and evidence-based criteria, including curated mutational hotspots, tumor-type associations, and gene-specific benign and pathogenic thresholds where defined.

### MAVE / SGE functional integration

- Ingests saturation genome editing / MAVE-derived functional scores for VHL variants, including quantitative measures of loss of function and mRNA dosage effects, and converts them into standardized functional evidence codes.
- Allows configurable score thresholds for assigning PS3 (supporting/moderate/strong) or BS3 evidence, enabling comparison of different calibration schemes or cut-points as the field evolves.

### Genome-wide VHL coverage

- Provides classifier outputs across the entire coding sequence of VHL for all assayed single-nucleotide variants, including missense, synonymous, and protein-truncating changes where data are available.
- Highlights regions with sparse editing or missing functional scores so users can contextualize classifier confidence.

### Transparent evidence view

- Presents, for each variant, the individual evidence codes triggered (functional, population, computational, segregation, phenotype, etc.), their strength, and the final ACMG/AMP classification.
- Designed to complement manual curation workflows by providing a structured starting point that expert curators can refine.

## Usage

1. **Input variant**
   - Enter a VHL variant using genomic or protein notation supported by the app (for example, HGVS-style cDNA or protein change) or select from a pre-loaded list of VHL variants if provided in the interface.

2. **Review evidence**
   - Inspect the automatically derived evidence codes, including functional scores, hotspot annotations, population data summaries, and any pre-integrated ClinGen VHL VCEP annotations.
   - Examine the functional MAVE/SGE score and its mapping to PS3/BS3 strength levels, along with any caveats about assay scope or performance.

3. **Interpret classification**
   - View the provisional ACMG/AMP classification (e.g., Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign) generated under the VHL VCEP specifications.
   - Optionally export or record the evidence summary for downstream use in curation systems, tumor boards, or research analyses, recognizing that this tool is decision-support only and not a stand-alone clinical reporting platform.

## Installation and Development

- Clone this repository and install Python dependencies (including Streamlit and the required scientific/ML packages) in a virtual environment.
- Run the app locally with the Streamlit CLI, for example:

## Intended Use and Limitations

This app is intended for use by clinical geneticists, molecular pathologists, variant curators, and researchers familiar with VHL and ACMG/AMP variant interpretation. It is a research and decision-support tool that reflects current public specifications and data sources, and any classifications should be reviewed in the context of full clinical, familial, and laboratory information before use in patient care.

The tool respects intellectual property and copyright of underlying publications and datasets, and users are expected to cite the original VHL VCEP specification papers and the VHL saturation genome editing/MAVE studies (e.g., Findlay et al. and subsequent VHL saturation genome editing work) when using classifier outputs in scientific or clinical work.