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

