# METRN_IVDD_transcriptomic_analysis

Public companion code repository for the manuscript:

**"METRN (Meteorin) in Intervertebral Disc Degeneration: an exploratory public transcriptomic analysis with translational relevance to adolescent idiopathic scoliosis (AIS)"**
(target journal: *BMC Musculoskeletal Disorders*)

> This is a **hypothesis-generating, exploratory re-analysis** of publicly available transcriptomic datasets. It is **not** a definitive mechanistic proof. The translational AIS extension is a conceptual framework built on the IVDD evidence base and does not include AIS-specific transcriptomic data.

## Data sources

All input data are publicly available from the Gene Expression Omnibus (GEO):

| GEO ID | Description | Samples used |
|--------|-------------|--------------|
| GSE124272 | IVDD peripheral blood | n = 16 |
| GSE150408 | IVDD peripheral blood | n = 34 (25 treatment samples excluded; 17 IDD + 17 volunteers retained) |
| GSE23130 | Human surgical disc specimens, Thompson-graded | High-grade (Thompson IV–V, n = 8) vs Low-grade (I–III, n = 15) |
| GSE56081 | Same platform as GSE124272 | Excluded from meta-analysis (possible overlapping cohort) |
| GSE311180 | Supplementary descriptive context | Not included in the main meta-analysis |

## Repository structure

```
R/         R scripts: limma differential expression (FDR), random-effects meta-analysis,
           signature scoring (GSVA), GEO inspection, supplementary FDR tables
python/    Python equivalents for FDR calculation and signature scoring
results/   Processed outputs: FDR-corrected expression tables, meta-analysis results,
           co-expression tables, signature correlations
figures/   Publication figures (PNG + vector PDF): Workflow, Forest plot, Scatter, Bar chart
```

## Methods summary

1. **Differential expression** — `limma` with Benjamini–Hochberg FDR correction
   (`R/01_differential_expression_FDR.R`, `R/00_FDR_step1_all_datasets.R`).
2. **Meta-analysis** — DerSimonian–Laird random-effects model across GSE124272 and GSE150408
   (`R/02_meta_analysis.R`).
3. **Co-expression** — Spearman correlation of METRN with all genes in GSE124272
   (`R/FDR_supplementary.R`).
4. **Signature scoring** — Gene set variation analysis (GSVA) of seven curated pathological
   signatures; Spearman correlation with BH-FDR
   (`R/03_signature_scoring.R`, `python/03_signature_scoring.py`).

## Key results

- METRN showed **no significant differential expression after FDR correction** in any single IVDD dataset (all FDR > 0.05).
- Random-effects meta-analysis: pooled log₂FC = −0.063 (95% CI −2.25 to 2.12; *P* = 0.778; *I*² = 68.6%).
- Only the **ER-stress signature** was significantly correlated with METRN expression
  (Spearman *r* = 0.800, FDR = 0.001) in GSE124272.

## Requirements

- **R** ≥ 4.2 with `limma`, `metafor`, `GSVA`, `GEOquery`
- **Python** ≥ 3.9 with `pandas`, `numpy`, `scipy`

## Licence

Code is released under the **MIT Licence**. Input data remain under the licences of their original GEO submitters.

## Citation

See `CITATION.cff`. A permanent archived version with a DOI is available via Zenodo:
**10.5281/zenodo.21261996** (https://doi.org/10.5281/zenodo.21261996).
This DOI was minted from GitHub release v1.0.1 via the Zenodo GitHub integration.

## Competing interests / Funding

See the accompanying manuscript. The authors declare no competing interests; no specific
funding was received for this analysis.

## Data / Code availability

Processed data and analysis scripts are available in this repository (release v1.0.0) and
permanently archived on Zenodo (DOI: **10.5281/zenodo.21261996**,
https://doi.org/10.5281/zenodo.21261996).
