# Release v2.0.0 — METRN in IVDD: independent transcriptomic validation & mechanistic extension

This release accompanies the manuscript **"METRN in Intervertebral Disc Degeneration: Integrated Public Transcriptomic Evidence and Translational Relevance to Adolescent Idiopathic Scoliosis"** (BMC Musculoskeletal Disorders submission, manuscript v10).

## What changed since v1.x

- **Author list corrected to the true manuscript authorship** (8 authors, Wangjing Hospital of China Academy of Chinese Medical Sciences). This aligns the repository with the published single-cell anchor (Zhang P et al., *J Cell Physiol* 2025;240:e31506 — Ref. [16]) and the plasma-biomarker patent (ZL202310215676.1 — Ref. [17]).
- **Competing Interests added**: co-author P.Z. (Ping Zhang, Department of Pathology) is a named inventor on Chinese patent ZL202310215676.1 (plasma METRN-based IVDD diagnostic biomarker). All other authors declare no competing interests.
- **Funding statement added**:
  - China Academy of Chinese Medical Sciences Scientific and Technological Innovation Project — major grant **CI2026A02025**
  - Wangjing Hospital High-Level TCM Hospital Cultivation Team project — **WJYY-PYTD-2025-09**
- **Ethics approval referenced**: WJEC-KT-2025-027-P002.
- **Analysis code and processed data are unchanged** from v1.x (no methodological alterations; this release is a metadata/attribution update).

## Repository contents

| Path | Description |
|------|-------------|
| `R/` | Differential-expression (limma) and Benjamini–Hochberg FDR correction scripts |
| `python/` | Signature scoring and figure-generation scripts |
| `results/` | Processed CSVs: differential expression, random-effects meta-analysis, co-expression, signature correlation |
| `figures/` | Publication-ready figures (PNG / PDF / TIFF at 300 DPI) |
| `README.md` | Study overview, methods summary, data availability |
| `CITATION.cff` | Machine-readable citation (real authors) |
| `.zenodo.json` | Zenodo metadata (real creators) |
| `LICENSE` | MIT |

## Data availability

All primary transcriptomic data are publicly available from GEO:
- GSE124272 (peripheral blood, IVDD vs. control)
- GSE150408 (peripheral blood, sciatica/IVDD vs. control)
- GSE23130 (annulus fibrosus tissue, degenerative vs. non-degenerative)

Processed results and analysis code are archived in this repository and on Zenodo (concept DOI: **10.5281/zenodo.2126196**).

## How to cite

> Dong J, Guo J, Zhang P, Wang Y, Yu H, Zhang L, Zhu W, Kong D. METRN in Intervertebral Disc Degeneration: Integrated Public Transcriptomic Evidence and Translational Relevance to Adolescent Idiopathic Scoliosis. *BMC Musculoskeletal Disorders* (submitted). Code & data: https://github.com/kongdeyu0819-rgb/METRN_IVDD_transcriptomic_analysis (v2.0.0); Zenodo DOI: 10.5281/zenodo.2126196.

## Notes for maintainers

- A new GitHub Release triggers the Zenodo GitHub integration to mint a **new Zenodo version** with the corrected author list. After publishing, verify the Zenodo record (https://zenodo.org/record/21261996) shows the 8 named creators and update the manuscript's DOI reference to the concept DOI (10.5281/zenodo.2126196) so it always resolves to the latest version.
