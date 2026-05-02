# =============================================================================
#  Mendelian Randomization Analysis: NMR Metabolites and Heart Failure
#  
#  For submission to: Lipids in Health and Disease
#  
#  Description:
#    Five-layer MR analytical framework assessing causal relationships between
#    249 NMR-quantified circulating metabolites and heart failure (HF):
#      1. Two-sample MR (IVW, MR-Egger, weighted median, weighted mode)
#      2. Sensitivity analyses (Cochran's Q, Egger intercept, Steiger, LOO)
#      3. MR-PRESSO outlier detection and correction
#      4. Multivariable MR adjusting for BMI and SBP
#      5. Two-step mediation MR (7 candidate mediators)
#  
#  Data sources:
#    Exposure  — UK Biobank NMR metabolomics GWAS (~115,000 Europeans)
#    Outcome 1 — HERMES all-cause HF (47,309 / 930,014)
#    Outcome 2 — FinnGen HF + CHD (HFrEF proxy; 8,876 / 197,780)
#    Outcome 3 — FinnGen HF strict (validation; 13,087 / 195,091)
#  
#  Requirements:
#    R >= 4.3, TwoSampleMR, MRPRESSO, ieugwasr, dplyr, tidyr,
#    ggplot2, ggrepel, cowplot, writexl, data.table
#    OpenGWAS JWT token in ~/.Renviron
#  
#  Authors: Jifei Cai, Nan Tong, Chenchen Hang, Junyuan Wu, Shubin Guo
#  Contact: cirjeef@163.com
# =============================================================================

# =============================================================================
# 0. ENVIRONMENT SETUP
# =============================================================================
# Set working directory (modify to your local path)
setwd("/Users/jifeicai/Desktop/accepted2")

# OpenGWAS API authentication
# ---------------------------
# The IEU OpenGWAS API requires a JWT token for data access.
# To configure:
#   1. Register at https://api.opengwas.io and obtain your JWT
#   2. Add to ~/.Renviron:  OPENGWAS_JWT="ey...your_token..."
#   3. Restart R/RStudio
# NEVER hardcode tokens in scripts intended for public repositories.

set.seed(20260501)
options(stringsAsFactors = FALSE)

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)
dir.create("rds_new", showWarnings = FALSE)

# =============================================================================
# 1. PACKAGES AND TOKEN VALIDATION
# =============================================================================
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
  library(data.table)
})

cat("\n>>> Validating OpenGWAS token ...\n")
tryCatch({
  u <- ieugwasr::user()
  cat("[OK] Token valid.\n")
  if (!is.null(u$uid))             cat("    UID    :", u$uid, "\n")
  if (!is.null(u$github_username)) cat("    GitHub :", u$github_username, "\n")
  cat("\n")
}, error = function(e) {
  stop("Token validation FAILED: ", conditionMessage(e),
       "\nPlease configure your JWT at https://api.opengwas.io ",
       "and add to ~/.Renviron")
})

# Execution log
run_log <- list()
log_step <- function(name, ok, note = "") {
  run_log[[name]] <<- list(ok = ok, note = note)
  cat(if (ok) "[OK] " else "[FAIL] ", name,
      if (nchar(note)) paste0(" -- ", note) else "", "\n", sep = "")
}

# =============================================================================
# Journal-quality plotting helpers
# =============================================================================
theme_journal <- function(base = 11) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      plot.title       = element_text(face = "bold", size = base + 1, hjust = 0),
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(color = "black"),
      legend.title     = element_text(face = "bold"),
      legend.position  = "bottom",
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
}

save_fig <- function(plot, name, w, h) {
  ggsave(file.path("figures", paste0(name, ".pdf")),
         plot, width = w, height = h, device = "pdf")
  ggsave(file.path("figures", paste0(name, ".tiff")),
         plot, width = w, height = h, dpi = 300,
         compression = "lzw", device = "tiff")
  ggsave(file.path("figures", paste0(name, ".png")),
         plot, width = w, height = h, dpi = 300)
  invisible(NULL)
}

# =============================================================================
# 2. LOAD PRE-COMPUTED STAGE-1 MR RESULTS
# =============================================================================
cat("\n>>> Loading pre-computed MR results ...\n")
df_allHF <- readRDS("results_allHF.rds")
df_hfref <- readRDS("results_HFrEF_proxy.rds")
df_finn  <- readRDS("results_FinnGen_HF.rds")

ivw1 <- df_allHF[df_allHF$method == "Inverse variance weighted", ]
ivw2 <- df_hfref[df_hfref$method == "Inverse variance weighted", ]
ivw3 <- df_finn[df_finn$method  == "Inverse variance weighted", ]

ivw1$fdr <- p.adjust(ivw1$pval, method = "BH")
ivw2$fdr <- p.adjust(ivw2$pval, method = "BH")
ivw3$fdr <- p.adjust(ivw3$pval, method = "BH")

n_tests   <- nrow(ivw1)
bonf_p    <- 0.05 / n_tests
bonf_mets <- ivw1$metabolite[ivw1$pval < bonf_p]

cat("Metabolites tested:", n_tests, "\n")
cat("Bonferroni threshold:", formatC(bonf_p, format = "e", digits = 2), "\n")
cat("Bonferroni-significant:", length(bonf_mets), "\n\n")

# =============================================================================
# 3. MODULE A: Multiple-testing correction (M_eff + category FDR)
# =============================================================================
cat(">>> Module A: multiple-testing burden\n")

# Assign metabolite categories
ivw1$category <- case_when(
  grepl("VLDL", ivw1$metabolite, ignore.case = TRUE)                    ~ "VLDL",
  grepl("IDL",  ivw1$metabolite, ignore.case = TRUE)                    ~ "IDL",
  grepl("LDL",  ivw1$metabolite, ignore.case = TRUE)                    ~ "LDL",
  grepl("HDL",  ivw1$metabolite, ignore.case = TRUE)                    ~ "HDL",
  grepl("FA|MUFA|PUFA|SFA|DHA|LA|Omega", ivw1$metabolite, TRUE)        ~ "Fatty acids",
  grepl("Ala|Gln|Gly|His|Ile|Leu|Phe|Tyr|Val", ivw1$metabolite, TRUE)  ~ "Amino acids",
  grepl("Glucose|Lac|Cit|Pyr|Glycerol", ivw1$metabolite, TRUE)         ~ "Glycolysis",
  grepl("AcAce|bOHbut|Ace",  ivw1$metabolite, TRUE)                    ~ "Ketones",
  grepl("GlycA|GlycB",       ivw1$metabolite, TRUE)                    ~ "Inflammation",
  TRUE                                                                  ~ "Other"
)
cat("Category counts:\n"); print(table(ivw1$category))

# Within-category FDR
ivw1 <- ivw1 %>%
  group_by(category) %>%
  mutate(fdr_within = p.adjust(pval, method = "BH")) %>%
  ungroup()

# Effective number of independent tests (Galwey method)
beta_mat <- df_allHF %>%
  filter(method %in% c("Inverse variance weighted", "MR Egger",
                       "Weighted median", "Weighted mode")) %>%
  select(metabolite, method, b) %>%
  pivot_wider(names_from = metabolite, values_from = b) %>%
  select(-method) %>% as.matrix()

M_eff <- tryCatch({
  beta_mat <- beta_mat[, colSums(is.na(beta_mat)) == 0, drop = FALSE]
  cor_mat  <- cor(beta_mat); cor_mat[is.na(cor_mat)] <- 0
  ev <- eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 0)
  if (sum(ev) == 0) NA else round((sum(sqrt(ev)))^2 / sum(ev), 1)
}, error = function(e) NA)

M_eff_lit <- 35  # Nightingale panel reference (Karjalainen et al. 2022 Nature)
cat("Empirical M_eff:", M_eff,
    "| Literature reference:", M_eff_lit,
    "| Bonferroni @ M_eff:", formatC(0.05/M_eff_lit, format = "e", digits = 2), "\n")

write.csv(ivw1, "tables/Table_S1_full_IVW_with_category.csv", row.names = FALSE)
log_step("Module A (M_eff + category FDR)", TRUE)

# =============================================================================
# 4. MODULE C: Cross-cohort comparison
# =============================================================================
cat("\n>>> Module C: cross-cohort merge\n")

cross <- merge(
  ivw1[, c("metabolite", "b", "se", "pval", "or", "or_lci95", "or_uci95",
           "fdr", "n_snp", "mean_F", "het_Q_pval", "egger_intercept_pval",
           "category")],
  ivw2[, c("metabolite", "b", "se", "pval", "or", "fdr")],
  by = "metabolite", suffixes = c("_HERMES", "_HFrEF"), all = TRUE)
cross <- merge(
  cross, ivw3[, c("metabolite", "b", "se", "pval", "or", "fdr")],
  by = "metabolite", all = TRUE)
colnames(cross)[which(colnames(cross) == "b")]    <- "b_FinnGen"
colnames(cross)[which(colnames(cross) == "se")]   <- "se_FinnGen"
colnames(cross)[which(colnames(cross) == "pval")] <- "pval_FinnGen"
colnames(cross)[which(colnames(cross) == "or")]   <- "or_FinnGen"
colnames(cross)[which(colnames(cross) == "fdr")]  <- "fdr_FinnGen"

cross$sig_HERMES <- case_when(
  cross$pval_HERMES < bonf_p ~ "Bonferroni",
  cross$fdr_HERMES  < 0.05   ~ "FDR",
  cross$pval_HERMES < 0.05   ~ "Nominal",
  TRUE                       ~ "NS")
cross$dir_consistent_FinnGen <- sign(cross$b_HERMES) == sign(cross$b_FinnGen)
cross$dir_consistent_HFrEF   <- sign(cross$b_HERMES) == sign(cross$b_HFrEF)
cross$amp_ratio_HFrEF <- with(cross, ifelse(
  metabolite %in% bonf_mets, b_HFrEF / b_HERMES, NA))

saveRDS(cross, "rds_new/results_cross_comparison.rds")
write.csv(cross, "tables/Table_S3_cross_cohort.csv", row.names = FALSE)
log_step("Module C (cross-cohort)", TRUE)

# =============================================================================
# 5. MODULE D: Sensitivity panel (offline base + online Steiger/LOO)
# =============================================================================
cat("\n>>> Module D: sensitivity panel\n")

pull_stat <- function(df, met, mth, col) {
  v <- df[df$metabolite == met & df$method == mth, col]
  if (length(v) == 0) NA else v[1]
}

# Offline base (from pre-computed results)
sens_panel <- lapply(bonf_mets, function(met) {
  data.frame(
    metabolite = met,
    n_snp     = pull_stat(df_allHF, met, "Inverse variance weighted", "n_snp"),
    mean_F    = pull_stat(df_allHF, met, "Inverse variance weighted", "mean_F"),
    ivw_b     = pull_stat(df_allHF, met, "Inverse variance weighted", "b"),
    ivw_p     = pull_stat(df_allHF, met, "Inverse variance weighted", "pval"),
    egger_b   = pull_stat(df_allHF, met, "MR Egger", "b"),
    egger_p   = pull_stat(df_allHF, met, "MR Egger", "pval"),
    wm_b      = pull_stat(df_allHF, met, "Weighted median", "b"),
    wm_p      = pull_stat(df_allHF, met, "Weighted median", "pval"),
    wmode_b   = pull_stat(df_allHF, met, "Weighted mode", "b"),
    wmode_p   = pull_stat(df_allHF, met, "Weighted mode", "pval"),
    Q_pval    = pull_stat(df_allHF, met, "Inverse variance weighted", "het_Q_pval"),
    egger_intercept_p = pull_stat(df_allHF, met, "Inverse variance weighted",
                                  "egger_intercept_pval"),
    stringsAsFactors = FALSE)
}) %>% bind_rows()

sens_panel$dir_concordant <- with(sens_panel,
  sign(ivw_b) == sign(egger_b) & sign(ivw_b) == sign(wm_b))
sens_panel$pleiotropy_flag <- ifelse(
  !is.na(sens_panel$egger_intercept_p) & sens_panel$egger_intercept_p < 0.05,
  "DIRECTIONAL_PLEIOTROPY", "OK")

# Online: Steiger directionality + leave-one-out
online_extra <- list()
for (met in bonf_mets) {
  cat("  Steiger + LOO:", gsub("met-d-", "", met), "... ")
  r <- tryCatch({
    exp <- extract_instruments(outcomes = met, p1 = 5e-08, clump = TRUE)
    out <- extract_outcome_data(snps = exp$SNP, outcomes = "ebi-a-GCST009541")
    dat <- harmonise_data(exp, out); dat <- dat[dat$mr_keep, ]
    if (nrow(dat) < 4) { cat("skip (n<4)\n"); return(NULL) }
    dat$samplesize.exposure <- 115000
    dat$samplesize.outcome  <- 977323
    st  <- directionality_test(dat)
    loo <- mr_leaveoneout(dat)
    cat("done\n")
    data.frame(metabolite = met,
               steiger_dir = st$correct_causal_direction,
               steiger_p   = st$steiger_pval,
               loo_max_p   = max(loo$p, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("ERR -- sleeping 5s\n"); Sys.sleep(5); NULL
  })
  if (!is.null(r)) online_extra[[met]] <- r
  Sys.sleep(2)
}

if (length(online_extra) > 0) {
  df_online <- bind_rows(online_extra)
  sens_panel <- left_join(sens_panel, df_online, by = "metabolite")
}

saveRDS(sens_panel, "rds_new/results_sensitivity_full.rds")
write.csv(sens_panel, "tables/Table_S4_sensitivity_full.csv", row.names = FALSE)
log_step("Module D (sensitivity)",
         length(online_extra) > 0,
         paste0(length(online_extra), "/", length(bonf_mets), " online OK"))

# =============================================================================
# 6. MODULE E: MR-PRESSO
# =============================================================================
cat("\n>>> Module E: MR-PRESSO\n")
presso_results <- list()
for (met in bonf_mets) {
  cat("  PRESSO:", gsub("met-d-", "", met), "... ")
  r <- tryCatch({
    exp <- extract_instruments(outcomes = met, p1 = 5e-08, clump = TRUE)
    out <- extract_outcome_data(snps = exp$SNP, outcomes = "ebi-a-GCST009541")
    dat <- harmonise_data(exp, out); dat <- dat[dat$mr_keep, ]
    if (nrow(dat) < 4) { cat("skip (n<4)\n"); return(NULL) }
    pr <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                    SdOutcome = "se.outcome", SdExposure = "se.exposure",
                    data = dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                    NbDistribution = 3000, SignifThreshold = 0.05)
    mr2 <- pr$`Main MR results`
    n_out <- ifelse(is.null(pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
                    0, length(pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`))
    cat("done (", n_out, "outliers)\n")
    data.frame(metabolite = met,
      raw_beta = mr2$`Causal Estimate`[1], raw_se = mr2$Sd[1],
      raw_p    = mr2$`P-value`[1],
      corr_beta = mr2$`Causal Estimate`[2], corr_se = mr2$Sd[2],
      corr_p   = mr2$`P-value`[2],
      n_outliers = n_out,
      global_p   = pr$`MR-PRESSO results`$`Global Test`$Pvalue,
      stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("ERR:", e$message, "\n"); Sys.sleep(5); NULL
  })
  if (!is.null(r)) presso_results[[met]] <- r
  Sys.sleep(2)
}
df_presso <- bind_rows(presso_results)
if (nrow(df_presso) > 0) {
  df_presso$met_short <- gsub("met-d-", "", df_presso$metabolite)
  saveRDS(df_presso, "rds_new/results_MRPRESSO.rds")
  write.csv(df_presso, "tables/Table_S5_MRPRESSO.csv", row.names = FALSE)
}
log_step("Module E (MR-PRESSO)",
         nrow(df_presso) > 0,
         paste0(nrow(df_presso), "/", length(bonf_mets), " metabolites"))

# =============================================================================
# 7. MODULE F: Multivariable MR (adjusted for BMI + SBP)
# =============================================================================
cat("\n>>> Module F: MVMR adjusted for BMI + SBP\n")
mvmr_mets <- c("met-d-XXL_VLDL_CE_pct", "met-d-VLDL_C", "met-d-L_VLDL_TG_pct")
mvmr_results <- list()
for (met in mvmr_mets) {
  cat("  MVMR:", gsub("met-d-", "", met), "... ")
  r <- tryCatch({
    e  <- mv_extract_exposures(id_exposure = c(met, "ieu-b-40", "ieu-b-38"))
    o  <- extract_outcome_data(snps = e$SNP, outcomes = "ebi-a-GCST009541")
    d  <- mv_harmonise_data(e, o)
    mv <- mv_multiple(d)$result
    mv$model <- met
    cat("done\n"); mv
  }, error = function(e) {
    cat("ERR:", e$message, "\n"); Sys.sleep(5); NULL
  })
  if (!is.null(r)) mvmr_results[[met]] <- r
  Sys.sleep(3)
}
df_mvmr <- bind_rows(mvmr_results)
if (nrow(df_mvmr) > 0) {
  df_mvmr$or       <- exp(df_mvmr$b)
  df_mvmr$or_lci95 <- exp(df_mvmr$b - 1.96 * df_mvmr$se)
  df_mvmr$or_uci95 <- exp(df_mvmr$b + 1.96 * df_mvmr$se)
  saveRDS(df_mvmr, "rds_new/results_MVMR.rds")
  write.csv(df_mvmr, "tables/Table_S6_MVMR.csv", row.names = FALSE)
}
log_step("Module F (MVMR)", nrow(df_mvmr) > 0,
         paste0(length(mvmr_results), "/", length(mvmr_mets), " models"))

# =============================================================================
# 8. MODULE G: Two-step mediation MR (expanded mediator panel)
# =============================================================================
cat("\n>>> Module G: mediation MR (expanded)\n")
focal_mets <- c("met-d-XXL_VLDL_CE_pct", "met-d-VLDL_C", "met-d-L_VLDL_TG_pct")
mediators <- list(
  list(id = "ieu-b-35",           name = "CRP"),
  list(id = "ieu-b-40",           name = "BMI"),
  list(id = "ieu-b-38",           name = "SBP"),
  list(id = "ebi-a-GCST90002232", name = "Fasting glucose"),
  list(id = "ieu-b-110",          name = "LDL_C"),
  list(id = "ieu-b-108",          name = "ApoB"))
outcome_id <- "ebi-a-GCST009541"
ivw_main   <- df_allHF[df_allHF$method == "Inverse variance weighted", ]

med_2step <- function(met_id, med_id, med_name, total_b) {
  cat("  ", gsub("met-d-", "", met_id), "->", med_name, "... ")
  tryCatch({
    # Step A: metabolite -> mediator
    ea <- extract_instruments(outcomes = met_id, p1 = 5e-08, clump = TRUE)
    if (is.null(ea) || nrow(ea) < 3) { cat("StepA snp<3\n"); return(NULL) }
    oa <- extract_outcome_data(snps = ea$SNP, outcomes = med_id)
    if (is.null(oa) || nrow(oa) == 0) { cat("StepA no data\n"); return(NULL) }
    da <- harmonise_data(ea, oa); da <- da[da$mr_keep, ]
    if (nrow(da) < 3) { cat("StepA harm<3\n"); return(NULL) }
    ra <- mr(da, method_list = "mr_ivw")

    # Step B: mediator -> HF
    eb <- extract_instruments(outcomes = med_id, p1 = 5e-08, clump = TRUE)
    if (is.null(eb) || nrow(eb) < 3) { cat("StepB snp<3\n"); return(NULL) }
    ob <- extract_outcome_data(snps = eb$SNP, outcomes = outcome_id)
    if (is.null(ob) || nrow(ob) == 0) { cat("StepB no data\n"); return(NULL) }
    db <- harmonise_data(eb, ob); db <- db[db$mr_keep, ]
    if (nrow(db) < 3) { cat("StepB harm<3\n"); return(NULL) }
    rb <- mr(db, method_list = "mr_ivw")

    # Indirect effect (product of coefficients) with delta-method SE
    ind    <- ra$b * rb$b
    ind_se <- sqrt(ra$b^2 * rb$se^2 + rb$b^2 * ra$se^2)
    ind_p  <- 2 * pnorm(abs(ind / ind_se), lower.tail = FALSE)
    pm <- ifelse(!is.na(total_b) && total_b != 0 && sign(ind) == sign(total_b),
                 round(ind / total_b * 100, 1), NA)
    cat("p =", formatC(ind_p, format = "e", digits = 2),
        " %med =", ifelse(is.na(pm), "NA", paste0(pm, "%")), "\n")
    data.frame(metabolite = met_id, mediator = med_name,
               beta_a = ra$b, se_a = ra$se, p_a = ra$pval,
               beta_b = rb$b, se_b = rb$se, p_b = rb$pval,
               indirect_effect = ind, indirect_se = ind_se, indirect_p = ind_p,
               total_effect = total_b, prop_mediated = pm,
               stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("ERR:", e$message, "\n"); Sys.sleep(5); NULL
  })
}

med_results <- list()
for (met in focal_mets) {
  tb <- ivw_main$b[ivw_main$metabolite == met]
  if (length(tb) == 0) tb <- NA
  cat("\nMet:", gsub("met-d-", "", met), "(total b =", round(tb, 4), ")\n")
  for (md in mediators) {
    r <- med_2step(met, md$id, md$name, tb)
    if (!is.null(r)) med_results[[length(med_results) + 1]] <- r
    Sys.sleep(2)
  }
}
df_med <- bind_rows(med_results)
if (nrow(df_med) > 0) {
  df_med$met_short <- gsub("met-d-", "", df_med$metabolite)
  saveRDS(df_med, "rds_new/results_mediation_expanded.rds")
  write.csv(df_med, "tables/Table_S7_mediation_expanded.csv", row.names = FALSE)
}
log_step("Module G (mediation)", nrow(df_med) > 0,
         paste0(nrow(df_med), " pathways tested"))

# =============================================================================
# 9. FIGURE 2: Volcano plot
# =============================================================================
cat("\n>>> Figure 2: volcano\n")
ivw1$log10p <- -log10(ivw1$pval)
ivw1$log_or <- log(ivw1$or)
ivw1$significance <- factor(case_when(
  ivw1$pval < bonf_p ~ "Bonferroni",
  ivw1$fdr  < 0.05   ~ "FDR < 0.05",
  ivw1$pval < 0.05   ~ "p < 0.05",
  TRUE                ~ "NS"),
  levels = c("Bonferroni", "FDR < 0.05", "p < 0.05", "NS"))

# Label Bonferroni-significant + extreme effect-size metabolites
ivw1$label <- ifelse(
  ivw1$pval < bonf_p | abs(ivw1$log_or) > 0.25,
  gsub("met-d-", "", ivw1$metabolite), "")

p_volcano <- ggplot(ivw1, aes(x = log_or, y = log10p, color = significance)) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted",
             color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(bonf_p), linetype = "dashed",
             color = "#C0392B", linewidth = 0.4) +
  geom_point(alpha = 0.75, size = 1.8) +
  scale_color_manual(values = c("Bonferroni" = "#C0392B",
                                "FDR < 0.05" = "#2980B9",
                                "p < 0.05"   = "#27AE60",
                                "NS"         = "grey75")) +
  geom_text_repel(aes(label = label), size = 2.8, color = "black",
                  segment.color = "grey60", max.overlaps = 30,
                  min.segment.length = 0) +
  labs(x = expression(log(OR) ~ "per SD increase in metabolite"),
       y = expression(-log[10](italic(P))), color = NULL,
       title = paste0("MR associations: ", nrow(ivw1),
                      " NMR metabolites vs all-cause HF (HERMES)")) +
  theme_journal()
save_fig(p_volcano, "Fig2_volcano_allHF", w = 8.5, h = 6.5)
log_step("Figure 2 (volcano)", TRUE)

# =============================================================================
# 10. FIGURE 3: Forest plot (top 15, with OR values and pleiotropy flags)
# =============================================================================
cat(">>> Figure 3: forest\n")
top15 <- ivw1 %>% arrange(pval) %>% head(15)
top15$met_label <- gsub("met-d-", "", top15$metabolite)

# Flag metabolites with directional pleiotropy (Egger intercept P < 0.05)
egger_int <- df_allHF %>%
  filter(method == "Inverse variance weighted") %>%
  select(metabolite, egger_intercept_pval)
top15 <- top15 %>% left_join(egger_int, by = "metabolite",
                             suffix = c("", ".dup"))
flag_idx <- which(!is.na(top15$egger_intercept_pval) &
                    top15$egger_intercept_pval < 0.05)
top15$met_label[flag_idx] <- paste0(top15$met_label[flag_idx], " *")

top15 <- top15 %>% arrange(or)
top15$met_label <- factor(top15$met_label, levels = top15$met_label)
top15$direction <- ifelse(top15$or > 1, "Risk", "Protective")
top15$or_text   <- sprintf("%.3f (%.3f-%.3f)", top15$or, top15$or_lci95, top15$or_uci95)

p_forest <- ggplot(top15, aes(x = or, y = met_label, color = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(xmin = or_lci95, xmax = or_uci95),
                width = 0.25, linewidth = 0.55, orientation = "y") +
  geom_point(size = 2.6) +
  geom_text(aes(x = max(top15$or_uci95) + 0.05, label = or_text),
            hjust = 0, size = 2.5, color = "black", show.legend = FALSE) +
  scale_color_manual(values = c("Risk" = "#C0392B", "Protective" = "#2980B9")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = "Odds Ratio (95% CI) per SD", y = NULL, color = NULL,
       title = "Top 15 metabolites causally associated with all-cause HF",
       caption = "* MR-Egger intercept P < 0.05 (potential directional pleiotropy)") +
  theme_journal() +
  theme(plot.caption = element_text(size = 8, color = "grey30", hjust = 0))
save_fig(p_forest, "Fig3_forest_top15", w = 9, h = 6)
log_step("Figure 3 (forest)", TRUE)

# =============================================================================
# 11. FIGURE 4: Scatter + funnel plots (clean axis labels, shared legend)
# =============================================================================
cat(">>> Figure 4: scatter + funnel\n")
plot_mets <- c("met-d-XXL_VLDL_CE_pct", "met-d-VLDL_C",
               "met-d-L_VLDL_TG_pct", "met-d-S_VLDL_C")
met_labels <- c("met-d-XXL_VLDL_CE_pct" = "XXL_VLDL_CE_pct",
                "met-d-VLDL_C"          = "VLDL_C",
                "met-d-L_VLDL_TG_pct"   = "L_VLDL_TG_pct",
                "met-d-S_VLDL_C"        = "S_VLDL_C")

panel_list <- list()
for (met in plot_mets) {
  short <- met_labels[met]
  cat("  ", short, "... ")
  r <- tryCatch({
    e <- extract_instruments(outcomes = met, p1 = 5e-08, clump = TRUE)
    o <- extract_outcome_data(snps = e$SNP, outcomes = "ebi-a-GCST009541")
    d <- harmonise_data(e, o)

    # Clean display names (remove OpenGWAS API IDs)
    d$exposure <- short
    d$outcome  <- "Heart failure"

    rs <- mr(d, method_list = c("mr_ivw", "mr_egger_regression",
                                "mr_weighted_median"))
    rs$exposure <- short
    rs$outcome  <- "Heart failure"

    sp <- mr_scatter_plot(rs, d)[[1]] +
      ggtitle(paste0(short, " -- Scatter")) +
      theme_journal(base = 9) +
      theme(legend.position = "none")

    ss <- mr_singlesnp(d)
    ss$exposure <- short
    ss$outcome  <- "Heart failure"
    fp <- mr_funnel_plot(ss)[[1]] +
      ggtitle(paste0(short, " -- Funnel")) +
      theme_journal(base = 9) +
      theme(legend.position = "none")

    cat("done\n")
    list(sp, fp)
  }, error = function(e) {
    cat("ERR:", e$message, "\n"); Sys.sleep(5); NULL
  })
  if (!is.null(r)) panel_list <- c(panel_list, r)
  Sys.sleep(2)
}

if (length(panel_list) > 0) {
  # Extract shared legend from one scatter panel
  shared_legend <- get_legend(
    panel_list[[1]] + theme(legend.position = "bottom"))
  grid <- plot_grid(plotlist = panel_list, ncol = 2,
                    labels = letters[seq_along(panel_list)],
                    label_size = 11, align = "hv")
  fig4 <- plot_grid(grid, shared_legend, ncol = 1,
                    rel_heights = c(1, 0.06))
  save_fig(fig4, "Fig4_scatter_funnel",
           w = 10, h = 3.5 * (length(panel_list) / 2))
}
log_step("Figure 4 (scatter+funnel)",
         length(panel_list) > 0,
         paste0(length(panel_list)/2, "/", length(plot_mets), " mets plotted"))

# =============================================================================
# 12. FIGURE 5: Cross-outcome comparison (with regression line + labels)
# =============================================================================
cat(">>> Figure 5: cross-outcome\n")
compare <- merge(
  ivw1[, c("metabolite", "b", "pval")],
  ivw2[, c("metabolite", "b", "pval")],
  by = "metabolite", suffixes = c("_allHF", "_HFrEF"))

compare$sig <- factor(case_when(
  compare$pval_allHF < bonf_p & compare$pval_HFrEF < 0.05/nrow(ivw2) ~ "Both significant",
  compare$pval_allHF < bonf_p ~ "All-HF only",
  compare$pval_HFrEF < 0.05/nrow(ivw2) ~ "HFrEF proxy only",
  TRUE ~ "NS"),
  levels = c("Both significant", "All-HF only", "HFrEF proxy only", "NS"))

compare$label <- ifelse(compare$sig == "Both significant",
                        gsub("met-d-", "", compare$metabolite), "")
r_val <- cor(compare$b_allHF, compare$b_HFrEF, use = "complete.obs")
slope_val <- coef(lm(b_HFrEF ~ b_allHF, data = compare))[2]
annot_text <- sprintf("r = %.2f, slope = %.2f", r_val, slope_val)

p_compare <- ggplot(compare, aes(x = b_allHF, y = b_HFrEF, color = sig)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.4) +
  geom_smooth(method = "lm", se = TRUE, color = "#E67E22",
              linewidth = 0.6, alpha = 0.15, data = compare) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(aes(label = label), size = 2.5, color = "black",
                  max.overlaps = 20, segment.color = "grey60") +
  scale_color_manual(values = c("Both significant" = "#C0392B",
                                "All-HF only"      = "#2980B9",
                                "HFrEF proxy only" = "#27AE60",
                                "NS"               = "grey75")) +
  annotate("text", x = min(compare$b_allHF) + 0.02,
           y = max(compare$b_HFrEF) - 0.02,
           label = annot_text, size = 3.5, hjust = 0, color = "#E67E22") +
  labs(x = "MR effect (\u03B2) on all-cause HF",
       y = "MR effect (\u03B2) on HF + CHD (HFrEF proxy)",
       color = NULL,
       title = "Cross-outcome effect comparison (slope = 1 dashed)") +
  theme_journal()
save_fig(p_compare, "Fig5_subtype_compare", w = 7.5, h = 7)
log_step("Figure 5 (cross-outcome)", TRUE)

# =============================================================================
# 13. FIGURE 6: Mediation MR (3-metabolite faceted panel)
# =============================================================================
if (nrow(df_med) > 0) {
  cat(">>> Figure 6: mediation summary (expanded)\n")

  med_plot <- df_med %>%
    filter(metabolite %in% focal_mets) %>%
    mutate(
      sig = (p_a < 0.05 & p_b < 0.05 & indirect_p < 0.05),
      met_short = factor(gsub("met-d-", "", metabolite),
                         levels = c("VLDL_C", "XXL_VLDL_CE_pct", "L_VLDL_TG_pct")),
      pm_label = ifelse(!is.na(prop_mediated) & sig,
                        paste0(prop_mediated, "%"), "")
    )

  p_med <- ggplot(med_plot,
                  aes(x = reorder(mediator, indirect_effect),
                      y = indirect_effect, fill = sig)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = indirect_effect - 1.96 * indirect_se,
                      ymax = indirect_effect + 1.96 * indirect_se),
                  width = 0.25) +
    geom_text(aes(label = pm_label,
                  y = indirect_effect + sign(indirect_effect) *
                    (1.96 * indirect_se + 0.005)),
              size = 2.8, color = "#27AE60", fontface = "bold") +
    scale_fill_manual(values = c("TRUE" = "#27AE60", "FALSE" = "grey70"),
                      labels = c("FALSE" = "NS", "TRUE" = "Significant"),
                      name = "Mediation") +
    facet_wrap(~ met_short, ncol = 1, scales = "free_x") +
    labs(x = "Candidate mediator",
         y = "Indirect (mediated) effect on HF",
         title = "Two-step mediation MR: metabolite \u2192 mediator \u2192 HF") +
    coord_flip(clip = "off") +
    theme_journal() +
    theme(strip.text = element_text(face = "bold", size = 10),
          plot.margin = margin(5, 20, 5, 5))
  save_fig(p_med, "Fig6_mediation_expanded", w = 8, h = 10)
  log_step("Figure 6 (mediation)", TRUE)
} else {
  log_step("Figure 6 (mediation)", FALSE, "no mediation data")
}

# =============================================================================
# 14. SUPPLEMENTARY WORKBOOK
# =============================================================================
cat("\n>>> Writing supplementary workbook\n")

# S2: IV summary by metabolite category
s2_summary <- ivw1 %>%
  group_by(category) %>%
  summarise(
    n_metabolites = n(),
    median_nSNP   = median(n_snp, na.rm = TRUE),
    median_F      = round(median(mean_F, na.rm = TRUE), 1),
    n_bonf_sig    = sum(pval < bonf_p, na.rm = TRUE),
    n_fdr_sig     = sum(fdr < 0.05, na.rm = TRUE),
    .groups = "drop")

sheets <- list(
  "S1_full_IVW"     = ivw1       %>% arrange(pval),
  "S2_IV_summary"   = s2_summary,
  "S3_cross_cohort" = cross      %>% arrange(pval_HERMES),
  "S4_sensitivity"  = sens_panel)
if (exists("df_presso") && nrow(df_presso) > 0) sheets[["S5_MRPRESSO"]] <- df_presso
if (exists("df_mvmr")   && nrow(df_mvmr)   > 0) sheets[["S6_MVMR"]]    <- df_mvmr
if (exists("df_med")    && nrow(df_med)     > 0) sheets[["S7_mediation"]] <- df_med

write_xlsx(sheets, "tables/Supplementary_Tables_LipidsHealthDis.xlsx")
cat("  Supplementary workbook written (", length(sheets), " sheets)\n")

# =============================================================================
# 15. RUN SUMMARY
# =============================================================================
cat("\n=========================================================\n")
cat(" RUN SUMMARY\n")
cat("=========================================================\n")
for (nm in names(run_log)) {
  s <- run_log[[nm]]
  cat(sprintf(" %-35s %s %s\n",
              nm,
              if (s$ok) "[OK]  " else "[FAIL]",
              if (nchar(s$note)) s$note else ""))
}

cat("\n KEY STATISTICS FOR MANUSCRIPT:\n")
cat(sprintf("   M_eff (effective tests):  %s\n", M_eff))
cat(sprintf("   Bonferroni threshold:     P < %.2e (at n=%d)\n", bonf_p, n_tests))
cat(sprintf("   Bonferroni-significant:   %d metabolites (all VLDL)\n", length(bonf_mets)))

if (nrow(df_med) > 0) {
  sig_med <- df_med %>% filter(p_a < 0.05, p_b < 0.05, indirect_p < 0.05)
  if (nrow(sig_med) > 0) {
    cat("\n   Significant mediation pathways:\n")
    for (i in seq_len(nrow(sig_med))) {
      r <- sig_med[i, ]
      cat(sprintf("     %s -> %s -> HF:  %%med = %.1f%%  P = %.2e\n",
                  gsub("met-d-", "", r$metabolite),
                  r$mediator, r$prop_mediated, r$indirect_p))
    }
  }
}

cat("\n Output directories: figures/  tables/  rds_new/\n")
cat("=========================================================\n")
