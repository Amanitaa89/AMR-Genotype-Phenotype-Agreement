# AMR-Genotype-Phenotype-Agreement
# MIC ↔ WGS (ARG) Class-Level Agreement Analysis

## Overview
This repository provides a **validated, reproducible framework** to assess **genotype–phenotype agreement** between antimicrobial susceptibility testing (MIC-based phenotype) and whole-genome sequencing (WGS)-detected antimicrobial resistance genes (ARGs), evaluated at the **antimicrobial class level**.

The workflow combines:
- **R** for rule-based, auditable data processing
- **Excel** for standardized calculation and documentation of agreement metrics

This approach is designed for **AMR surveillance, epidemiological analysis, and method development**, and is **not intended to replace phenotypic AST**.

---

## Methodological Concept

For each isolate and antimicrobial class:

- **Phenotype (class-level)**  
  `1` if **any antibiotic** in the class is resistant  
  `0` otherwise

- **Genotype (class-level)**  
  `1` if **any resistance gene** defined in the class rulebook is present  
  `0` otherwise

Agreement is evaluated using confusion-matrix–based metrics and chance-corrected statistics.

---

## Repository Structure

---

## Input Data Requirements

### Input file (CSV or Excel)
- One row per isolate
- Binary phenotype data (0 = susceptible, 1 = resistant)
- Binary genotype data (0 = gene absent, 1 = gene present)

Example:

| LIMS_ID | Ciprofloxacin | Levofloxacin | gyrA_S83Y | qnrS13 | tet(A) |
|--------|---------------|--------------|-----------|--------|--------|
| ID001  | 1             | 1            | 1         | 0      | 0      |

> MIC interpretation and breakpoint application must be completed **before** importing data into R.

---

## Role of R in the Workflow

R is used as the **primary analytical engine** to ensure:
- Reproducibility
- Consistent rule application
- Traceable generation of confusion matrix counts

R performs:
- Class-level phenotype and genotype construction
- Application of rulebook logic
- Generation of **TP, FP, TN, FN** per antimicrobial class

Excel is used **only** for calculation of performance metrics and reporting.

---

## Confusion Matrix Definitions

For each antimicrobial class:

- **TP (True Positive)**  
  Phenotype = 1 AND Genotype = 1

- **FP (False Positive)**  
  Phenotype = 0 AND Genotype = 1

- **TN (True Negative)**  
  Phenotype = 0 AND Genotype = 0

- **FN (False Negative)**  
  Phenotype = 1 AND Genotype = 0

---

## Performance Metrics Calculated

Metrics are calculated automatically in Excel using fixed formulas:

- Sensitivity
- Specificity
- Positive Predictive Value (PPV)
- Negative Predictive Value (NPV)
- Accuracy
- Cohen’s Kappa (chance-corrected agreement)
- Phi coefficient (Pearson correlation for binary data)

> High accuracy in classes with rare resistance does **not** imply strong agreement and must be interpreted alongside Kappa and event frequency.

---

## Workflow Summary

1. Prepare input dataset (phenotype + genotype)
2. Apply class-level rules using the R script
3. Export TP, FP, TN, FN per class
4. Enter counts into the Excel template (`Class_Raw_Counts`)
5. Review calculated metrics and archive results

---

## Interpretation Notes

- Sensitivity is interpreted only when sufficient resistant isolates are present
- Classes with near-zero resistant events should be reported as **not evaluable**
- Kappa and Phi quantify agreement beyond chance and should not be interpreted in isolation

---

## Reproducibility and Quality Control

- R scripts are version-controlled
- Excel formulas are fixed and auditable
- Input data, scripts, and outputs should be archived together
- Recommended for use as a controlled laboratory template

---

## Intended Use

- AMR surveillance programs
- Research and method validation
- Epidemiological genotype–phenotype analysis

---

## License

---

## Citation
This workflow, is an internal method "so far"

