## One Model To Classify Them All: A Semantically Enriched VHL/VCEP Classifier

A Streamlit application for ACMG/AMP classification of germline *VHL* variants, tuned to the ClinGen VHL VCEP criteria specification and enriched with structured case- and functional evidence from large-scale VHL biocuration efforts.

---

## Authors

- **Dr. Marc McCoy, MD/MSc** – Computational Geneticist, Resident Physician  
- **Veronica Andric**  
- **Sarah Ridd**  
- **Elif Tuzlali**  
- **Helia Purnaghshband**  
- **Kamalika Bhandari Deka**  
- **Kirsten Farncombe**
- **Dr. Raymond Kim, MD/PhD** – Principal Investigator  

---

## Overview

VHL/VCEP Classifier is an ML/AI‑powered, gene‑centric decision‑support LLM application explicitly designed as a guideline‑forming intelligence system for Von Hippel–Lindau (*VHL*) variant interpretation. It operationalizes the ClinGen VHL Variant Curation Expert Panel (VCEP) gene‑ and disease‑specific ACMG/AMP specifications to deliver transparent, standards‑compliant classifications of *VHL* variants. The system combines a rules‑based implementation of VHL‑tailored criteria with multiplexed functional readouts from saturation genome editing / MAVE experiments to generate consistent, reproducible provisional classifications across the *VHL* coding region while exposing the full evidence stack for expert review and methods development.

---

## Key Features

- **VHL VCEP–conformant rules engine**  
  Implements the ClinGen VHL VCEP ACMG/AMP criteria, including gene‑ and disease‑specific thresholds and code usage, rather than generic ACMG/AMP logic.

- **Integrated case‑level and phenotype evidence**  
  Supports incorporation of structured evidence from VHL‑focused biocuration efforts (e.g., mutational hotspots, tumor‑type–specific patterns, genotype–phenotype correlations, and disease‑specific annotation protocols).

- **Functional genomics integration**  
  Ingests quantitative MAVE / saturation genome editing data for *VHL* to derive PS3/BS3 strengths at configurable cut‑points, enabling systematic exploration of functional evidence calibration and VUS resolution.

- **Evidence‑and‑label export**  
  Produces machine‑readable evidence profiles (triggered codes, strengths, links to underlying data) and final ACMG/AMP labels suitable for benchmarking LLMs and training downstream ML models.

---

## Academic Rationale and Ground Truth

The classifier serves both as a reference implementation of the written VHL VCEP guidelines and as an experimental platform for studying how those guidelines behave when applied at scale. It encodes the VHL‑specific ACMG/AMP framework, including adapted pathogenic and benign evidence codes, as curated and versioned in the ClinGen Criteria Specification Registry. By integrating phenotype‑driven and evidence‑based features such as VHL mutational hotspots, tumor‑type–specific associations, gene‑level population thresholds, and a disease‑specific annotation protocol, the system closely mirrors contemporary VCEP expert practice rather than relying on generic ACMG/AMP rules.

The application is explicitly designed to generate high‑fidelity “ground truth” evidence profiles and labels, enabling:  
- benchmarking of LLMs and algorithmic classifiers against expert‑derived decisions  
- training of ML/AI models on richly structured, guideline‑aligned supervision  
- empirical refinement of gene‑ and disease‑specific guideline criteria over time.

---

## Functional Genomics Platform

The application serves as a testbed for integrating high‑throughput functional genomics into clinical variant interpretation. It ingests quantitative data from *VHL* saturation genome editing / MAVE studies—capturing HIF‑dependent cellular fitness and mRNA dosage effects for thousands of single‑nucleotide variants—and converts these scores into structured functional evidence that can be mapped to PS3/BS3 at configurable strengths.

This enables users to:  
- explore alternative PS3/BS3 thresholds and their effect on classification distributions  
- quantify the impact of MAVE‑derived evidence on reclassification, particularly VUS resolution  
- assess how dense, gene‑wide functional maps might be incorporated into future ACMG/AMP and ClinGen/VCEP recommendations.

---

## Operational Context

The VHL/VCEP Classifier provides genome‑wide coverage for assayed *VHL* SNVs, including missense, synonymous, and protein‑truncating variants, and flags regions with sparse or absent functional data to contextualize classifier confidence. For each variant, the interface presents:

1. Triggered evidence codes by ACMG/AMP category  
2. Their assigned strengths (supporting, moderate, strong, very strong, etc.)  
3. The resulting five‑tier ACMG/AMP classification

This yields a structured, machine‑readable evidence‑and‑label resource that supports:  
- rigorous comparative analysis of human versus AI curation  
- audits of evidence use and conflict resolution patterns  
- downstream data science on variant interpretation behavior.

---

## Installation and Usage

Once running, the web interface will guide you through:  
- entering or selecting a *VHL* variant  
- reviewing population, case‑level, and functional evidence  
- inspecting triggered criteria and the resulting ACMG/AMP classification  
- exporting evidence‑and‑label records for downstream analysis.

---

## Intended Users and Impact

This work builds directly on the efforts of the VHL early detection and annotation teams and on the leadership of **Dr. Raymond Kim**, Medical Director of Cancer Early Detection at Princess Margaret Cancer Centre and chair of the ClinGen VHL Variant Curation Expert Panel, whose group has defined much of the present VHL‑specific ACMG/AMP framework and curated many of the genotype–phenotype resources underlying this classifier.

The tool is targeted to:  
- clinical geneticists and molecular pathologists  
- variant curators and clinical laboratory scientists  
- researchers working at the interface of clinical genomics, functional genomics, and computational method development  

It is explicitly positioned as a **research‑grade** tool: all outputs must be interpreted in the context of complete clinical, familial, and laboratory information. Use in scholarly or clinical contexts should include appropriate citation of:  
- VHL VCEP criteria publications and ClinGen Criteria Specification Registry entries  
- disease‑specific annotation work and protocols  
- *VHL* saturation genome editing / MAVE studies and related functional genomics resources.

---

## Citation

If you use this classifier or its derived evidence‑and‑label datasets in your work, please cite the relevant VHL VCEP guideline publications, the VHL disease‑specific annotation protocol, and the underlying *VHL* functional genomics studies, along with this repository.