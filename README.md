## One Model To Classify Them All: A Semantically Enriched VHL/VCEP Classifier

## Authors
*   **Dr. Marc McCoy, MD/MSc** - Computational Geneticist, Resident Physician
*   **Veronica Andric**
*   **Sarah Ridd**
*   **Elif Tuzlali**
*   **Helia Purnaghshband**
*   **Dr. Raymond Kim, MD/PhD** - Principal Investigator

### Overview

VHL/VCEP Classifier is an ML/AI‑powered, gene‑centric decision‑support LLM application explicitly designed as a guideline‑forming intelligence system for Von Hippel–Lindau (VHL) variant interpretation. It operationalizes the ClinGen VHL Variant Curation Expert Panel (VCEP) gene‑ and disease‑specific ACMG/AMP specifications to deliver transparent, standards‑compliant classifications of VHL variants, combining a rules‑based implementation of VHL‑tailored criteria with multiplexed functional readouts from saturation genome editing / MAVE experiments to generate consistent, reproducible provisional classifications across the VHL coding region while exposing the full evidence stack for expert review and methods development.

### Academic Rationale and Ground Truth

The classifier serves both as a reference implementation of the written VHL VCEP guidelines and as an experimental platform for studying how those guidelines behave when applied at scale. It encodes the VHL-specific ACMG/AMP framework, including adapted pathogenic and benign evidence codes, as curated and versioned in the ClinGen Criteria Specification Registry. By integrating phenotype-driven and evidence-based features such as VHL mutational hotspots, tumor-type–specific associations, gene-level population thresholds, and a disease-specific annotation protocol, the system closely mirrors contemporary VCEP expert practice rather than relying on generic ACMG/AMP rules. It generates high-fidelity “ground truth” evidence profiles and labels suitable for benchmarking LLMs, training ML/AI models, and future refinement of guideline criteria.

### Functional Genomics Platform

The application also serves as a testbed for integrating high‑throughput functional genomics into clinical variant interpretation. It ingests quantitative data from VHL saturation genome editing / MAVE studies—capturing HIF‑dependent cellular fitness and mRNA dosage effects for thousands of single‑nucleotide variants—and converts these scores into structured functional evidence that can be mapped to PS3/BS3 at configurable strengths. This enables systematic exploration of functional evidence calibration (e.g., alternative PS3/BS3 thresholds), quantification of the impact of MAVE-derived evidence on reclassification (especially VUS resolution), and empirical assessment of how dense, gene-wide functional maps might be incorporated into future ACMG/AMP and ClinGen/VCEP recommendations.

### Operational Context

The VHL/VCEP Classifier provides genome-wide coverage for assayed VHL SNVs, including missense, synonymous, and protein-truncating variants, and flags regions with sparse or absent functional data to contextualize classifier confidence. For each variant, the interface presents (i) triggered evidence codes by category, (ii) their assigned strengths, and (iii) the resulting ACMG/AMP classification, yielding a structured, machine-readable evidence-and-label resource. This enables rigorous comparative analysis of human versus AI curation as well as academic studies of evidence use and conflict resolution.

### Impact

This work builds directly on the efforts of the VHL early detection and annotation teams and on the leadership of Dr. Raymond Kim, Medical Director of Cancer Early Detection at Princess Margaret Cancer Centre and chair of the ClinGen VHL Variant Curation Expert Panel, whose group has defined much of the present VHL-specific ACMG/AMP framework and curated many of the genotype–phenotype resources underlying this classifier. It is targeted to clinical geneticists, molecular pathologists, variant curators, and researchers working at the interface of clinical genomics, functional genomics, and computational method development, and is explicitly positioned as a research-grade tool: all outputs must be interpreted in the context of complete clinical, familial, and laboratory information, with appropriate citation of VHL VCEP criteria publications, disease-specific annotation work, and VHL saturation genome editing/MAVE studies in scholarly and clinical use.