# DV90 scATAC-seq CNA Caller Benchmarking

Repository for the benchmarking of three CNA inference tools on scATAC-seq data from the DV90 lung cancer cell line, as part of the **onAir-CNA** project (WP3).

## Overview

Copy Number Alterations (CNAs) can be inferred from single-cell ATAC-seq data by exploiting the relationship between fragment depth and DNA copy number. This repository benchmarks three published tools:

| Tool | Strategy | Reference |
|---|---|---|
| [AtaCNV](https://github.com/XiDsLab/AtaCNV) | BICseq2 segmentation, GC correction | XiDsLab |
| [epiAneufinder](https://github.com/colomemaria/epiAneufinder) | Binary segmentation, GC correction | Marini et al. |
| [CONGASp](https://github.com/caravagnalab/CONGASp) | Bayesian mixture model | Patruno et al. 2023 |

All tools were run with **1 Mb genomic bins** and evaluated against the **CCLE ABSOLUTE ground truth** for DV90.

---

## Dataset

- **Cell line**: DV90 (lung cancer)
- **Platform**: 10x Genomics scMultiomics
- **Cells**: 3,472 (filtered barcodes)
- **Fragments**: 183,873,732 unique fragments (PCR-deduplicated)
- **Median fragments/cell**: 35,947

Fragment file and barcodes are **not included** in this repository due to size.

---

## Repository Structure

\`\`\`
toolEval/
├── Script/
│   ├── build_matrix.py          # Builds 1 Mb count matrix from fragment file
│   ├── 1_AtaCNV_DV90.ipynb      # AtaCNV pipeline (R kernel)
│   ├── 2_epiAneufinder.r        # epiAneufinder run script
│   ├── 3_congas.py              # CONGASp run script
│   ├── 4_Comparison.ipynb       # Benchmarking vs CCLE ground truth
│   ├── 5_Report_v4.ipynb        # PDF report generation
│   ├── 6_coverageAnalysis.r     # Coverage subsampling (N=100 x 4 levels)
│   └── 7_coverageStat.ipynb     # Coverage performance plot
├── Results/                     # Outputs, figures, CSV, PDF report
├── docker_image_parts/          # Docker image split in 90MB parts
├── Dockerfile                   # Original Dockerfile (see note)
├── reassemble.sh                # Reassembles and loads the Docker image
└── runMe.sh                     # Launches the JupyterLab container
\`\`\`

---

## Results Summary

| Tool | Sensitivity (gain) | Specificity (gain) | Sensitivity (loss) | Specificity (loss) |
|---|---|---|---|---|
| **AtaCNV** | **0.870** | **0.791** | 0.241 | 0.714 |
| epiAneufinder | 0.870 | 0.581 | 0.586 | 0.784 |
| CONGASp | 0.314 | 0.666 | 0.310 | 0.707 |

**→ AtaCNV is the recommended tool.**

### Most reliable CNA calls
- **Gain**: chr8, chr19 (AtaCNV + epiAneufinder, 2/3 tools)
- **Loss**: chr9, chr18, chr21 (all 3 tools, 3/3)

---

## Coverage Analysis

AtaCNV evaluated at 25%, 33%, 50%, 67% of full depth — N=100 random subsamples per level.

- **Gain sensitivity** stable across all coverage levels (~0.87)
- **Minimum recommended coverage**: 33% (~12,000 median fragments/cell)

---

## Docker Environment

\`\`\`bash
# Reassemble image and launch JupyterLab
bash runMe.sh
\`\`\`

Before running, set your password hash in \`runMe.sh\`:
\`\`\`bash
-e JUPYTER_PASSWORD_HASH='<insert_your_password_hash_here>'
\`\`\`

> ⚠️ The \`Dockerfile\` reflects the original build but does not reproduce the exact environment — additional packages were installed interactively. Use the pre-built image via \`reassemble.sh\` for full reproducibility.

---

## Note on PCR Duplicate Correction

The 10x fragment file column 5 (\`readSupport\`) counts how many times each unique fragment was sequenced. An early version mistakenly used this as a weight. All scripts now count each unique fragment exactly once.

---

## Authors

- **Luca Stormreig** — analysis, benchmarking, pipeline
- **Chiara Castelli** — Docker, reproducibility, generalization
- **Raffaele Calogero** (PI) — project design and supervision

## Project

[onAir-CNA](https://github.com/onAir-CNA) — WP3: CNA calling from scATAC-seq
