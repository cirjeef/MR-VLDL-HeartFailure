# R Analysis Code for Manuscript Submission

## Files

1. **MR_HF_LipidsHealthDisease_final.R** — Main analysis (Modules A–G, Figures 2–6)
2. **MR_HF_Supplementary_Analyses.R** — Supplementary analyses (Modules H–L, Figures S1–S4)

## Running Order

```r
# Step 1: Main analysis (~60-90 min, requires OpenGWAS API)
source("MR_HF_LipidsHealthDisease_final.R")

# Step 2: Supplementary analysis (~60-90 min, requires OpenGWAS API)
# Optional: download LDSC files for full MRlap correction
source("MR_HF_Supplementary_Analyses.R")
```

## Requirements

- R >= 4.3
- Packages: TwoSampleMR, MRPRESSO, ieugwasr, dplyr, tidyr, ggplot2, ggrepel, cowplot, writexl, data.table, patchwork
- Optional: MRlap, GenomicSEM (for full sample overlap correction)
- OpenGWAS JWT token in ~/.Renviron

## Outputs

| Directory | Contents |
|-----------|----------|
| figures/ | Main text figures (PDF + TIFF + PNG) |
| figures_supp/ | Supplementary figures |
| tables/ | Main supplementary tables |
| tables_supp/ | Additional supplementary tables (S8–S11) |
| rds_new/ | R data objects for reproducibility |
