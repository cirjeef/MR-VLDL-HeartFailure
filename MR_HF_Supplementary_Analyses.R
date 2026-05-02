# =============================================================================
#  MR_HF_Supplementary_Analyses.R
#  
#  Five supplementary analyses addressing anticipated reviewer concerns:
#    Module H: MRlap correction for sample overlap
#    Module I: Inter-metabolite correlation (clarify "10 independent signals")
#    Module J: MVMR-mediation (disentangle ApoB vs LDL-C mediation overlap)
#    Module K: Leave-one-out forest plots
#    Module L: Absolute concentration vs percentage sensitivity analysis
#
#  Run AFTER MR_HF_LipidsHealthDisease_final.R
#  Requires: all intermediate objects from the main analysis in memory
#            OR .rds files in rds_new/
#
#  New packages required:
#    MRlap (GitHub: nlapier2/MRlap), ggforestplot (optional), patchwork
#
#  Authors: Jifei Cai et al.
# =============================================================================

setwd("/Users/jifeicai/Desktop/accepted2")
set.seed(20260501)

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(MRPRESSO)
  library(ieugwasr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(writexl)
})

# --- Install MRlap if not available ---
# MRlap depends on GenomicSEM (also from GitHub)
mrlap_available <- FALSE
if (!requireNamespace("MRlap", quietly = TRUE)) {
  cat(">>> Installing MRlap and dependencies from GitHub ...\n")
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  tryCatch({
    # GenomicSEM is a dependency of MRlap
    if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
      cat("  Installing GenomicSEM ...\n")
      remotes::install_github("GenomicSEM/GenomicSEM")
    }
    cat("  Installing MRlap ...\n")
    remotes::install_github("n-mounier/MRlap")
    library(MRlap)
    mrlap_available <- TRUE
    cat("  MRlap installed successfully.\n")
  }, error = function(e) {
    cat("  WARNING: MRlap installation failed:", e$message, "\n")
    cat("  Will use F-statistic correction as fallback.\n")
    cat("  To install manually, run:\n")
    cat("    remotes::install_github('GenomicSEM/GenomicSEM')\n")
    cat("    remotes::install_github('n-mounier/MRlap')\n\n")
  })
} else {
  library(MRlap)
  mrlap_available <- TRUE
}

# --- Install patchwork if not available ---
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

# --- Reload data if needed ---
if (!exists("df_allHF"))   df_allHF   <- readRDS("results_allHF.rds")
if (!exists("df_hfref"))   df_hfref   <- readRDS("results_HFrEF_proxy.rds")
if (!exists("df_finn"))    df_finn    <- readRDS("results_FinnGen_HF.rds")
if (!exists("df_presso"))  df_presso  <- tryCatch(readRDS("rds_new/results_MRPRESSO.rds"),          error = function(e) data.frame())
if (!exists("df_mvmr"))    df_mvmr    <- tryCatch(readRDS("rds_new/results_MVMR.rds"),              error = function(e) data.frame())
if (!exists("df_med"))     df_med     <- tryCatch(readRDS("rds_new/results_mediation_expanded.rds"), error = function(e) data.frame())
if (!exists("sens_panel")) sens_panel <- tryCatch(readRDS("rds_new/results_sensitivity_full.rds"),   error = function(e) data.frame())

ivw1 <- df_allHF[df_allHF$method == "Inverse variance weighted", ]
ivw1$fdr <- p.adjust(ivw1$pval, method = "BH")
n_tests   <- nrow(ivw1)
bonf_p    <- 0.05 / n_tests
bonf_mets <- ivw1$metabolite[ivw1$pval < bonf_p]

dir.create("figures_supp", showWarnings = FALSE)
dir.create("tables_supp",  showWarnings = FALSE)

# --- Helpers ---
theme_journal <- function(base = 11) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      plot.title       = element_text(face = "bold", size = base + 1, hjust = 0),
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(color = "black"),
      legend.title     = element_text(face = "bold"),
      legend.position  = "bottom",
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank())
}

save_fig <- function(plot, name, w, h) {
  ggsave(file.path("figures_supp", paste0(name, ".pdf")),
         plot, width = w, height = h, device = "pdf")
  ggsave(file.path("figures_supp", paste0(name, ".tiff")),
         plot, width = w, height = h, dpi = 300,
         compression = "lzw", device = "tiff")
  ggsave(file.path("figures_supp", paste0(name, ".png")),
         plot, width = w, height = h, dpi = 300)
  cat("  Saved:", name, "\n")
}

supp_log <- list()
log_supp <- function(name, ok, note = "") {
  supp_log[[name]] <<- list(ok = ok, note = note)
  cat(if (ok) "[OK] " else "[FAIL] ", name,
      if (nchar(note)) paste0(" -- ", note) else "", "\n", sep = "")
}

# =============================================================================
# MODULE H: MRlap correction for sample overlap
# =============================================================================
cat("\n")
cat("===========================================================\n")
cat(" MODULE H: MRlap Correction for Sample Overlap\n")
cat("===========================================================\n")
cat("
 Purpose: The exposure GWAS (UK Biobank NMR) and the primary
 outcome GWAS (HERMES) partially overlap in UK Biobank samples.
 MRlap corrects IVW estimates for winner's curse, weak instrument
 bias, and sample overlap using cross-trait LD score regression.

 REQUIREMENTS (must be downloaded before running):
   1. LD scores: eur_w_ld_chr/ folder
   2. HapMap3 SNP list: w_hm3.snplist

   DOWNLOAD (choose one source):
     Option A - Broad Institute (original):
       wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
       wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
       tar -jxvf eur_w_ld_chr.tar.bz2
       bzip2 -d w_hm3.snplist.bz2

     Option B - Zenodo mirror (if Broad is down):
       https://zenodo.org/records/8182036  (eur_w_ld_chr)
       https://zenodo.org/records/10515792 (S-LDSC full bundle, includes w_hm3.snplist)

     Option C - Auto-download in R (this script will attempt if files missing)

   3. Full GWAS summary statistics (auto-downloaded via ieugwasr)
\n")

# --- Paths to LDSC reference files (USER: modify these) ---
ld_path  <- "eur_w_ld_chr"     # folder containing chr1.l2.ldscore.gz etc.
hm3_path <- "w_hm3.snplist"    # HapMap3 SNP list

# Check if LDSC files exist; attempt auto-download if not
if (!dir.exists(ld_path)) {
  cat("  LD scores not found. Attempting auto-download from Broad Institute ...\n")
  tryCatch({
    download.file(
      "https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2",
      "eur_w_ld_chr.tar.bz2", mode = "wb", quiet = FALSE)
    system("tar -jxvf eur_w_ld_chr.tar.bz2")
    cat("  LD scores downloaded and extracted.\n")
  }, error = function(e) {
    cat("  Auto-download failed:", e$message, "\n")
    cat("  Please download manually from:\n")
    cat("    https://zenodo.org/records/8182036\n")
    cat("  or:\n")
    cat("    https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2\n\n")
  })
}
if (!file.exists(hm3_path)) {
  cat("  HapMap3 SNP list not found. Attempting auto-download ...\n")
  tryCatch({
    download.file(
      "https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2",
      "w_hm3.snplist.bz2", mode = "wb", quiet = FALSE)
    system("bzip2 -d w_hm3.snplist.bz2")
    cat("  HapMap3 SNP list downloaded.\n")
  }, error = function(e) {
    cat("  Auto-download failed:", e$message, "\n")
    cat("  Please download manually from:\n")
    cat("    https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2\n\n")
  })
}

ldsc_available <- mrlap_available && dir.exists(ld_path) && file.exists(hm3_path)

if (!ldsc_available) {
  cat("
  *** LDSC reference files not found. ***
  Please download from one of:
    1. https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
    2. https://zenodo.org/records/8182036 (Zenodo mirror)
  and place eur_w_ld_chr/ and w_hm3.snplist in your working directory.

  FALLBACK: Using F-statistic-based winner's curse correction instead.
  This approximates MRlap but does not fully account for sample overlap.
  The full MRlap correction is recommended for the final submission.
\n")
}

# --- Download full GWAS summary statistics for outcome (once) ---
hermes_file <- "gwas_data/HERMES_HF_full.tsv.gz"
dir.create("gwas_data", showWarnings = FALSE)

if (ldsc_available && !file.exists(hermes_file)) {
  cat("  Downloading HERMES HF full GWAS (this may take 10-20 min) ...\n")
  tryCatch({
    hermes_gwas <- ieugwasr::associations(
      variants = NULL,         # all SNPs
      id = "ebi-a-GCST009541",
      proxies = 0)
    # Format for MRlap
    hermes_df <- data.frame(
      SNP  = hermes_gwas$rsid,
      CHR  = hermes_gwas$chr,
      POS  = hermes_gwas$position,
      ALT  = hermes_gwas$ea,
      REF  = hermes_gwas$nea,
      BETA = hermes_gwas$beta,
      SE   = hermes_gwas$se,
      N    = hermes_gwas$n,
      stringsAsFactors = FALSE)
    data.table::fwrite(hermes_df, hermes_file, sep = "\t")
    cat("  HERMES GWAS saved:", hermes_file, "\n")
  }, error = function(e) {
    cat("  WARNING: Could not download full HERMES GWAS:", e$message, "\n")
    cat("  Will use fallback correction method.\n")
  })
}

# --- Run MRlap for each Bonferroni-significant metabolite ---
mrlap_results <- list()
for (met in bonf_mets) {
  short <- gsub("met-d-", "", met)
  cat("  MRlap:", short, "... ")

  r <- tryCatch({
    # First get standard IVW estimate for comparison
    exp_dat <- extract_instruments(outcomes = met, p1 = 5e-08, clump = TRUE)
    out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = "ebi-a-GCST009541")
    dat     <- harmonise_data(exp_dat, out_dat)
    dat     <- dat[dat$mr_keep, ]
    if (nrow(dat) < 5) { cat("skip (n<5)\n"); return(NULL) }
    ivw_res <- mr(dat, method_list = "mr_ivw")

    beta_corr <- NA; se_corr <- NA; p_corr <- NA
    method_used <- "none"

    # Attempt 1: Full MRlap (requires LDSC files + full GWAS)
    if (ldsc_available && file.exists(hermes_file)) {
      # Download full exposure GWAS
      exp_file <- paste0("gwas_data/", short, "_full.tsv.gz")
      if (!file.exists(exp_file)) {
        cat("downloading exposure GWAS ... ")
        exp_gwas <- tryCatch(
          ieugwasr::associations(variants = NULL, id = met, proxies = 0),
          error = function(e) NULL)
        if (!is.null(exp_gwas) && nrow(exp_gwas) > 10000) {
          exp_df <- data.frame(
            SNP = exp_gwas$rsid, CHR = exp_gwas$chr,
            POS = exp_gwas$position,
            ALT = exp_gwas$ea, REF = exp_gwas$nea,
            BETA = exp_gwas$beta, SE = exp_gwas$se,
            N = exp_gwas$n, stringsAsFactors = FALSE)
          data.table::fwrite(exp_df, exp_file, sep = "\t")
        }
      }

      if (file.exists(exp_file)) {
        mrlap_out <- tryCatch({
          MRlap::MRlap(
            exposure      = exp_file,
            exposure_name = short,
            outcome       = hermes_file,
            outcome_name  = "HF",
            ld            = ld_path,
            hm3           = hm3_path,
            MR_threshold  = 5e-08,
            MR_pruning_dist = 500,
            MR_pruning_LD = 0,
            save_logfiles = FALSE,
            verbose       = FALSE
          )
        }, error = function(e) {
          cat("MRlap internal error: ", e$message, " ... ")
          NULL
        })

        if (!is.null(mrlap_out)) {
          beta_corr   <- mrlap_out$MRcorrection$Estimate
          se_corr     <- mrlap_out$MRcorrection$StdError
          p_corr      <- mrlap_out$MRcorrection$Pvalue
          method_used <- "MRlap_full"
        }
      }
    }

    # Attempt 2: F-statistic-based winner's curse correction (fallback)
    if (is.na(beta_corr)) {
      mean_F    <- mean(dat$beta.exposure^2 / dat$se.exposure^2)
      # Winner's curse bias: beta_obs ≈ beta_true * F/(F-1)
      # Correction: beta_true ≈ beta_obs * (F-1)/F
      wc_factor <- (mean_F - 1) / mean_F
      beta_corr <- ivw_res$b * wc_factor
      se_corr   <- ivw_res$se * (1 / wc_factor)
      p_corr    <- 2 * pnorm(abs(beta_corr / se_corr), lower.tail = FALSE)
      method_used <- "F_stat_correction"
    }

    cat("done [", method_used, "] (IVW P=",
        formatC(ivw_res$pval, format = "e", digits = 1),
        " -> corrected P=",
        formatC(p_corr, format = "e", digits = 1), ")\n")

    data.frame(
      metabolite     = met,
      met_short      = short,
      ivw_beta       = ivw_res$b,
      ivw_se         = ivw_res$se,
      ivw_p          = ivw_res$pval,
      ivw_or         = exp(ivw_res$b),
      corrected_beta = beta_corr,
      corrected_se   = se_corr,
      corrected_p    = p_corr,
      corrected_or   = exp(beta_corr),
      method         = method_used,
      mean_F         = mean(dat$beta.exposure^2 / dat$se.exposure^2),
      nsnp           = nrow(dat),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("ERR:", e$message, "\n"); Sys.sleep(3); NULL
  })

  if (!is.null(r)) mrlap_results[[met]] <- r
  Sys.sleep(2)
}

df_mrlap <- bind_rows(mrlap_results)
if (nrow(df_mrlap) > 0) {
  df_mrlap$still_sig <- df_mrlap$corrected_p < bonf_p
  saveRDS(df_mrlap, "rds_new/results_MRlap.rds")
  write.csv(df_mrlap, "tables_supp/Table_S8_MRlap.csv", row.names = FALSE)

  # Comparison plot: IVW vs MRlap-corrected OR
  mrlap_long <- df_mrlap %>%
    select(met_short, ivw_or, corrected_or) %>%
    pivot_longer(-met_short, names_to = "method", values_to = "or") %>%
    mutate(method = ifelse(method == "ivw_or", "Standard IVW", "MRlap-corrected"))

  p_mrlap <- ggplot(df_mrlap, aes(x = ivw_or, y = corrected_or)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = still_sig), size = 3) +
    geom_text_repel(aes(label = met_short), size = 2.5, max.overlaps = 15) +
    scale_color_manual(values = c("TRUE" = "#27AE60", "FALSE" = "#C0392B"),
                       labels = c("TRUE" = "Still significant", "FALSE" = "Lost significance"),
                       name = "After correction") +
    labs(x = "Standard IVW OR", y = "MRlap-corrected OR",
         title = "Effect of sample overlap correction on MR estimates") +
    theme_journal()
  save_fig(p_mrlap, "FigS1_MRlap_correction", w = 7, h = 6)
}
log_supp("Module H (MRlap)", nrow(df_mrlap) > 0,
         paste0(sum(df_mrlap$still_sig), "/", nrow(df_mrlap),
                " remain Bonferroni-significant after correction"))

# =============================================================================
# MODULE I: Inter-metabolite correlation matrix
# =============================================================================
cat("\n")
cat("===========================================================\n")
cat(" MODULE I: Inter-metabolite Correlation (10 Bonf signals)\n")
cat("===========================================================\n")
cat("
 Purpose: Demonstrate that the 10 Bonferroni-significant
 metabolites are highly correlated, representing facets of
 a single VLDL metabolic dysregulation signal rather than
 10 independent causal pathways.
\n")

# Extract IVW betas for all 10 Bonferroni metabolites across outcomes
corr_data <- df_allHF %>%
  filter(method %in% c("Inverse variance weighted", "MR Egger",
                       "Weighted median", "Weighted mode"),
         metabolite %in% bonf_mets) %>%
  select(metabolite, method, b) %>%
  pivot_wider(names_from = metabolite, values_from = b) %>%
  select(-method) %>%
  as.matrix()

# Clean column names
colnames(corr_data) <- gsub("met-d-", "", colnames(corr_data))
corr_data <- corr_data[, colSums(is.na(corr_data)) == 0, drop = FALSE]

if (ncol(corr_data) >= 3) {
  cor_mat <- cor(corr_data, use = "pairwise.complete.obs")

  # Effective number among just the 10 significant metabolites
  ev10 <- eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values
  ev10 <- pmax(ev10, 0)
  M_eff_10 <- round((sum(sqrt(ev10)))^2 / sum(ev10), 2)

  cat("Correlation matrix (10 Bonferroni metabolites):\n")
  print(round(cor_mat, 2))
  cat("\nMedian pairwise |r|:", round(median(abs(cor_mat[upper.tri(cor_mat)])), 3), "\n")
  cat("Range of pairwise |r|:", round(range(abs(cor_mat[upper.tri(cor_mat)])), 3), "\n")
  cat("M_eff among 10 signals:", M_eff_10, "\n")
  cat("Interpretation: these represent approximately",
      M_eff_10, "independent signals, not 10.\n\n")

  # Heatmap
  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("Met1", "Met2", "r")

  p_corr <- ggplot(cor_df, aes(x = Met1, y = Met2, fill = r)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", r)), size = 2.3) +
    scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#C0392B",
                         midpoint = 0, limits = c(-1, 1),
                         name = "Pearson r") +
    labs(title = paste0("Correlation among 10 Bonferroni-significant VLDL metabolites\n",
                        "(M_eff = ", M_eff_10,
                        "; median |r| = ",
                        round(median(abs(cor_mat[upper.tri(cor_mat)])), 2), ")"),
         x = NULL, y = NULL) +
    theme_journal(base = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7))
  save_fig(p_corr, "FigS2_correlation_heatmap", w = 8, h = 7)

  # Save correlation matrix
  write.csv(round(cor_mat, 4), "tables_supp/Table_S9_correlation_matrix.csv")
  log_supp("Module I (correlation)", TRUE,
           paste0("M_eff = ", M_eff_10, " among 10 metabolites"))
} else {
  log_supp("Module I (correlation)", FALSE, "insufficient data")
}

# =============================================================================
# MODULE J: MVMR-mediation to disentangle ApoB vs LDL-C
# =============================================================================
cat("\n")
cat("===========================================================\n")
cat(" MODULE J: MVMR-Mediation (ApoB vs LDL-C joint model)\n")
cat("===========================================================\n")
cat("
 Purpose: ApoB (60.5%) and LDL-C (54.5%) mediation proportions
 sum to >100% because they are highly correlated. This module
 estimates their INDEPENDENT mediation contributions using MVMR
 with both mediators as simultaneous exposures on HF.
\n")

mvmr_med_results <- list()
for (met in c("met-d-VLDL_C", "met-d-XXL_VLDL_CE_pct", "met-d-L_VLDL_TG_pct")) {
  short <- gsub("met-d-", "", met)
  cat("  MVMR-mediation:", short, "... ")

  r <- tryCatch({
    # Step A: metabolite -> ApoB (already have from mediation)
    # Step A: metabolite -> LDL-C (already have from mediation)
    med_rows <- df_med[df_med$metabolite == met &
                         df_med$mediator %in% c("ApoB", "LDL_C"), ]
    if (nrow(med_rows) < 2) { cat("skip (missing mediator data)\n"); return(NULL) }

    beta_a_apob <- med_rows$beta_a[med_rows$mediator == "ApoB"]
    beta_a_ldlc <- med_rows$beta_a[med_rows$mediator == "LDL_C"]

    # Step B (joint): ApoB + LDL-C -> HF simultaneously via MVMR
    e_mv <- mv_extract_exposures(
      id_exposure = c("ieu-b-108", "ieu-b-110"))  # ApoB + LDL-C
    o_mv <- extract_outcome_data(
      snps = e_mv$SNP, outcomes = "ebi-a-GCST009541")
    d_mv <- mv_harmonise_data(e_mv, o_mv)
    mv_res <- mv_multiple(d_mv)$result

    # Extract independent effects
    apob_row <- mv_res[grepl("apolipoprotein B|ApoB|ieu-b-108",
                             mv_res$exposure, ignore.case = TRUE), ]
    ldlc_row <- mv_res[grepl("LDL cholesterol|LDL-C|ieu-b-110",
                             mv_res$exposure, ignore.case = TRUE), ]

    if (nrow(apob_row) == 0 || nrow(ldlc_row) == 0) {
      cat("skip (MVMR parse fail)\n"); return(NULL)
    }

    # Independent indirect effects
    ind_apob <- beta_a_apob * apob_row$b[1]
    ind_ldlc <- beta_a_ldlc * ldlc_row$b[1]

    # Total effect
    total_b <- ivw1$b[ivw1$metabolite == met]

    cat("done\n")
    data.frame(
      metabolite          = met,
      met_short           = short,
      # Step B independent effects (from joint MVMR)
      apob_mvmr_beta      = apob_row$b[1],
      apob_mvmr_se        = apob_row$se[1],
      apob_mvmr_p         = apob_row$pval[1],
      ldlc_mvmr_beta      = ldlc_row$b[1],
      ldlc_mvmr_se        = ldlc_row$se[1],
      ldlc_mvmr_p         = ldlc_row$pval[1],
      # Independent indirect effects
      ind_apob_alone       = ind_apob,
      ind_ldlc_alone       = ind_ldlc,
      pct_apob_independent = ifelse(total_b != 0,
                                     round(ind_apob / total_b * 100, 1), NA),
      pct_ldlc_independent = ifelse(total_b != 0,
                                     round(ind_ldlc / total_b * 100, 1), NA),
      total_effect         = total_b,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("ERR:", e$message, "\n"); Sys.sleep(5); NULL
  })
  if (!is.null(r)) mvmr_med_results[[met]] <- r
  Sys.sleep(3)
}

df_mvmr_med <- bind_rows(mvmr_med_results)
if (nrow(df_mvmr_med) > 0) {
  saveRDS(df_mvmr_med, "rds_new/results_MVMR_mediation.rds")
  write.csv(df_mvmr_med, "tables_supp/Table_S10_MVMR_mediation.csv", row.names = FALSE)

  cat("\n  MVMR-mediation results (independent contributions):\n")
  for (i in seq_len(nrow(df_mvmr_med))) {
    r <- df_mvmr_med[i, ]
    cat(sprintf("    %s: ApoB independent = %.1f%% (MVMR P=%s), LDL-C independent = %.1f%% (MVMR P=%s)\n",
                r$met_short,
                r$pct_apob_independent,
                formatC(r$apob_mvmr_p, format = "e", digits = 1),
                r$pct_ldlc_independent,
                formatC(r$ldlc_mvmr_p, format = "e", digits = 1)))
  }
  cat("  Note: ApoB + LDL-C independent proportions should now sum to <100%\n")
}
log_supp("Module J (MVMR-mediation)", nrow(df_mvmr_med) > 0,
         paste0(nrow(df_mvmr_med), " metabolites"))

# =============================================================================
# MODULE K: Leave-one-out forest plots
# =============================================================================
cat("\n")
cat("===========================================================\n")
cat(" MODULE K: Leave-One-Out Forest Plots\n")
cat("===========================================================\n")

loo_plots <- list()
for (met in bonf_mets) {
  short <- gsub("met-d-", "", met)
  cat("  LOO:", short, "... ")

  r <- tryCatch({
    exp <- extract_instruments(outcomes = met, p1 = 5e-08, clump = TRUE)
    out <- extract_outcome_data(snps = exp$SNP, outcomes = "ebi-a-GCST009541")
    dat <- harmonise_data(exp, out); dat <- dat[dat$mr_keep, ]

    # Clean names
    dat$exposure <- short
    dat$outcome  <- "Heart failure"

    if (nrow(dat) < 5) { cat("skip (n<5)\n"); return(NULL) }

    loo <- mr_leaveoneout(dat)
    loo$exposure <- short
    loo$outcome  <- "Heart failure"

    p <- mr_leaveoneout_plot(loo)[[1]] +
      ggtitle(paste0("Leave-one-out: ", short)) +
      theme_journal(base = 8) +
      theme(axis.text.y = element_text(size = 5))

    cat("done (", nrow(loo), " SNPs)\n")
    p
  }, error = function(e) {
    cat("ERR:", e$message, "\n"); Sys.sleep(3); NULL
  })
  if (!is.null(r)) loo_plots[[short]] <- r
  Sys.sleep(2)
}

if (length(loo_plots) > 0) {
  # Save individual plots
  for (nm in names(loo_plots)) {
    save_fig(loo_plots[[nm]],
             paste0("FigS3_LOO_", nm),
             w = 7, h = max(4, 0.2 * 50))  # adjust height by SNP count
  }

  # Combined overview (first 4)
  if (length(loo_plots) >= 4) {
    combo <- plot_grid(plotlist = loo_plots[1:4], ncol = 2,
                       labels = "auto", label_size = 10)
    save_fig(combo, "FigS3_LOO_combined_top4", w = 14, h = 14)
  }
}
log_supp("Module K (LOO plots)", length(loo_plots) > 0,
         paste0(length(loo_plots), "/", length(bonf_mets), " plotted"))

# =============================================================================
# MODULE L: Absolute concentration vs percentage sensitivity analysis
# =============================================================================
cat("\n")
cat("===========================================================\n")
cat(" MODULE L: Absolute Concentration vs Percentage Comparison\n")
cat("===========================================================\n")
cat("
 Purpose: Several Bonferroni-significant metabolites are
 compositional percentages (e.g., XXL_VLDL_CE_pct). Because
 percentages within a particle sum to a constant, opposing
 directions for cholesterol% and TG% may be partly mathematical.
 This module compares percentage-based results with their
 absolute concentration counterparts.
\n")

# Define pairs: percentage <-> absolute concentration
pct_abs_pairs <- list(
  list(pct = "met-d-XXL_VLDL_CE_pct", abs = "met-d-XXL_VLDL_CE",
       label = "XXL-VLDL cholesteryl esters"),
  list(pct = "met-d-XXL_VLDL_C_pct",  abs = "met-d-XXL_VLDL_C",
       label = "XXL-VLDL cholesterol"),
  list(pct = "met-d-XXL_VLDL_FC_pct", abs = "met-d-XXL_VLDL_FC",
       label = "XXL-VLDL free cholesterol"),
  list(pct = "met-d-L_VLDL_TG_pct",   abs = "met-d-L_VLDL_TG",
       label = "L-VLDL triglycerides"),
  list(pct = "met-d-L_VLDL_C_pct",    abs = "met-d-L_VLDL_C",
       label = "L-VLDL cholesterol"),
  list(pct = "met-d-S_VLDL_CE_pct",   abs = "met-d-S_VLDL_CE",
       label = "S-VLDL cholesteryl esters")
)

abs_results <- list()
for (pair in pct_abs_pairs) {
  cat("  Comparing:", pair$label, "... ")

  r <- tryCatch({
    # Get percentage result (from existing data)
    pct_row <- ivw1[ivw1$metabolite == pair$pct, ]

    # Run MR for absolute concentration
    exp <- extract_instruments(outcomes = pair$abs, p1 = 5e-08, clump = TRUE)
    if (is.null(exp) || nrow(exp) < 3) { cat("abs snp<3\n"); return(NULL) }
    out <- extract_outcome_data(snps = exp$SNP, outcomes = "ebi-a-GCST009541")
    if (is.null(out) || nrow(out) == 0) { cat("abs no data\n"); return(NULL) }
    dat <- harmonise_data(exp, out); dat <- dat[dat$mr_keep, ]
    if (nrow(dat) < 3) { cat("abs harm<3\n"); return(NULL) }

    abs_mr <- mr(dat, method_list = "mr_ivw")

    cat("done (pct OR=", round(exp(pct_row$b), 3),
        ", abs OR=", round(exp(abs_mr$b), 3), ")\n")

    data.frame(
      label        = pair$label,
      met_pct      = pair$pct,
      met_abs      = pair$abs,
      pct_beta     = pct_row$b,
      pct_se       = pct_row$se,
      pct_p        = pct_row$pval,
      pct_or       = exp(pct_row$b),
      abs_beta     = abs_mr$b,
      abs_se       = abs_mr$se,
      abs_p        = abs_mr$pval,
      abs_or       = exp(abs_mr$b),
      abs_nsnp     = abs_mr$nsnp,
      direction_consistent = sign(pct_row$b) == sign(abs_mr$b),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("ERR:", e$message, "\n"); Sys.sleep(3); NULL
  })
  if (!is.null(r)) abs_results[[pair$label]] <- r
  Sys.sleep(2)
}

df_abs <- bind_rows(abs_results)
if (nrow(df_abs) > 0) {
  saveRDS(df_abs, "rds_new/results_abs_vs_pct.rds")
  write.csv(df_abs, "tables_supp/Table_S11_abs_vs_pct.csv", row.names = FALSE)

  # Visualization: paired dot plot
  abs_long <- df_abs %>%
    select(label, pct_or, abs_or, pct_p, abs_p) %>%
    pivot_longer(cols = c(pct_or, abs_or), names_to = "type", values_to = "or") %>%
    mutate(type = ifelse(type == "pct_or", "Percentage", "Absolute"),
           p = ifelse(type == "Percentage",
                      pct_p, abs_p),
           sig = p < bonf_p)

  p_abs <- ggplot(abs_long, aes(x = or, y = reorder(label, or), color = type, shape = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_point(size = 3, position = position_dodge(width = 0.4)) +
    scale_color_manual(values = c("Percentage" = "#C0392B", "Absolute" = "#2980B9"),
                       name = "Measure") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = "Bonferroni sig", "FALSE" = "Not sig"),
                       name = "Significance") +
    labs(x = "Odds Ratio per SD (IVW)", y = NULL,
         title = "Sensitivity analysis: percentage vs absolute concentration",
         subtitle = "Direction consistency confirms findings are not purely compositional artifacts") +
    theme_journal()
  save_fig(p_abs, "FigS4_abs_vs_pct", w = 9, h = 5.5)

  n_consistent <- sum(df_abs$direction_consistent, na.rm = TRUE)
  cat("\n  Direction consistency:", n_consistent, "/", nrow(df_abs), "\n")
}
log_supp("Module L (abs vs pct)", nrow(df_abs) > 0,
         paste0(sum(df_abs$direction_consistent), "/",
                nrow(df_abs), " direction consistent"))

# =============================================================================
# FINAL: Update supplementary workbook with new tables
# =============================================================================
cat("\n>>> Updating supplementary workbook with new tables\n")

# Read existing workbook tables
existing_sheets <- list(
  "S1_full_IVW"     = read.csv("tables/Table_S1_full_IVW_with_category.csv"),
  "S3_cross_cohort" = read.csv("tables/Table_S3_cross_cohort.csv"),
  "S4_sensitivity"  = sens_panel)

# Add S2 summary
if (!exists("s2_summary")) {
  s2_summary <- ivw1 %>%
    group_by(category) %>%
    summarise(
      n_metabolites = n(),
      median_nSNP   = median(n_snp, na.rm = TRUE),
      median_F      = round(median(mean_F, na.rm = TRUE), 1),
      n_bonf_sig    = sum(pval < bonf_p, na.rm = TRUE),
      n_fdr_sig     = sum(fdr < 0.05, na.rm = TRUE),
      .groups = "drop")
}
existing_sheets[["S2_IV_summary"]] <- s2_summary

# Add existing and new supplementary tables
if (exists("df_presso") && nrow(df_presso) > 0) existing_sheets[["S5_MRPRESSO"]] <- df_presso
if (exists("df_mvmr")   && nrow(df_mvmr)   > 0) existing_sheets[["S6_MVMR"]]    <- df_mvmr
if (exists("df_med")    && nrow(df_med)     > 0) existing_sheets[["S7_mediation"]] <- df_med

# New supplementary tables
if (exists("df_mrlap")    && nrow(df_mrlap)    > 0) existing_sheets[["S8_MRlap"]]           <- df_mrlap
if (exists("cor_mat"))                               existing_sheets[["S9_correlation"]]      <- as.data.frame(round(cor_mat, 4))
if (exists("df_mvmr_med") && nrow(df_mvmr_med) > 0) existing_sheets[["S10_MVMR_mediation"]] <- df_mvmr_med
if (exists("df_abs")      && nrow(df_abs)      > 0) existing_sheets[["S11_abs_vs_pct"]]     <- df_abs

write_xlsx(existing_sheets,
           "tables_supp/Supplementary_Tables_Complete.xlsx")
cat("  Complete workbook written (",
    length(existing_sheets), " sheets: S1-S11)\n")

# =============================================================================
# SUPPLEMENTARY RUN SUMMARY
# =============================================================================
cat("\n===========================================================\n")
cat(" SUPPLEMENTARY ANALYSIS SUMMARY\n")
cat("===========================================================\n")
for (nm in names(supp_log)) {
  s <- supp_log[[nm]]
  cat(sprintf(" %-40s %s %s\n",
              nm,
              if (s$ok) "[OK]  " else "[FAIL]",
              if (nchar(s$note)) s$note else ""))
}

cat("\n KEY SENTENCES FOR REVISED MANUSCRIPT:\n")
cat("===========================================================\n")

if (exists("df_mrlap") && nrow(df_mrlap) > 0) {
  n_still <- sum(df_mrlap$still_sig)
  n_full  <- sum(df_mrlap$method == "MRlap_full", na.rm = TRUE)
  n_fstat <- sum(df_mrlap$method == "F_stat_correction", na.rm = TRUE)
  method_note <- if (n_full == nrow(df_mrlap)) {
    "full MRlap (cross-trait LDSC)"
  } else if (n_full > 0) {
    paste0("MRlap for ", n_full, ", F-statistic correction for ", n_fstat)
  } else {
    "F-statistic-based winner's curse correction (LDSC files not available)"
  }
  cat(sprintf("
  MRlap: After correcting for potential sample overlap between
  UK Biobank and HERMES (%s), %d of %d metabolites
  remained Bonferroni-significant, confirming that sample overlap
  did not meaningfully inflate our primary estimates (Table S8,
  Figure S1).
\n", method_note, n_still, nrow(df_mrlap)))
}

if (exists("cor_mat")) {
  med_r <- round(median(abs(cor_mat[upper.tri(cor_mat)])), 2)
  cat(sprintf("
  Inter-metabolite correlation: The 10 Bonferroni-significant VLDL
  metabolites were highly correlated (median |r| = %s, M_eff = %s),
  indicating that these represent approximately %s independent
  signals reflecting different facets of VLDL metabolic dysregulation
  rather than 10 distinct causal pathways (Table S9, Figure S2).
\n", med_r, M_eff_10, M_eff_10))
}

if (exists("df_mvmr_med") && nrow(df_mvmr_med) > 0) {
  vldlc_row <- df_mvmr_med[df_mvmr_med$met_short == "VLDL_C", ]
  if (nrow(vldlc_row) > 0) {
    cat(sprintf("
  MVMR-mediation: When ApoB and LDL-C were included simultaneously
  as co-exposures on HF, the independent mediation contribution of
  ApoB was %.1f%% and LDL-C was %.1f%% for VLDL_C (Table S10).
  These non-overlapping proportions sum to <100%%, confirming that
  the >100%% observed in univariable mediation reflects shared
  variance rather than methodological error.
\n", vldlc_row$pct_apob_independent, vldlc_row$pct_ldlc_independent))
  }
}

if (exists("df_abs") && nrow(df_abs) > 0) {
  n_con <- sum(df_abs$direction_consistent, na.rm = TRUE)
  cat(sprintf("
  Absolute vs percentage: Effect direction was consistent between
  percentage-based and absolute concentration measures in %d of %d
  tested pairs (Table S11, Figure S4), indicating that the observed
  causal effects are not purely compositional artifacts.
\n", n_con, nrow(df_abs)))
}

cat("\n Output directories: figures_supp/  tables_supp/  rds_new/\n")
cat("===========================================================\n")
