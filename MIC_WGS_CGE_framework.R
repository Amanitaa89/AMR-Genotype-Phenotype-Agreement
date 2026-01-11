#!/usr/bin/env Rscript
# ============================================================
# MIC ↔ WGS (ARG) ↔ CGE Class-Level Agreement Analysis
# Author: (Your Name / Lab)
# ============================================================

suppressPackageStartupMessages({
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
  if (!requireNamespace("tools", quietly = TRUE)) install.packages("tools")
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
  library(dplyr)
  library(readr)
  library(tools)
  library(stringr)
})

# -----------------------------
# 0) User settings (EDIT HERE)
# -----------------------------
# Input file: MIC + gene presence/absence (0/1) per isolate
INPUT_FILE <- "/Users/aat/Documents/Enteritidis/New_MIC_WGS_Data.csv"

# ID column name in MIC/WGS file
ID_COL <- "LIMS ID"

# CGE report (LIMS ID, Genotype, Predicted Phenotype, CGE Predicted Phenotype)
CGE_FILE <- "/Users/aat/Documents/Enteritidis/WGS_CGE_Report.csv"

# Output parent directory
OUT_PARENT <- "/Users/aat/Documents/MIC_WGS_Kappa_Output"

# Definitions used:
PHENO_DEF <- "Phenotype_class = 1 if ANY antibiotic in class == 1, else 0"
GENO_DEF  <- "Genotype_class  = 1 if ANY gene in class rule present == 1, else 0"
CGE_DEF   <- "CGE_class       = 1 if CGE predicted phenotype lists ANY antibiotic in class, else 0"

# -----------------------------
# 1) Class rulebook (EDIT HERE)
# -----------------------------
# Use EXACT column names from your MIC/WGS input file.
class_rules <- list(
  FQ_Quinolones = list(
    antibiotics = c("Ciprofloxacin","Levofloxacin","NalidixicAcid"),
    genes = c("gyrA_S83Y","gyrA_S83F","gyrA_D87G","gyrA_D87N","gyrA_D87Y","parC_S80I","qnrS13","qnrB19")
  ),
  ESBL_Cephalosporins = list(
    antibiotics = c("Cefotaxime","Ceftriaxone","Ceftazidime","Cefepime","Ceftiofur"),
    genes = c("blaCTX-M-8","blaCTX-M-65","blaSHV-12")
  ),
  Tetracyclines = list(
    antibiotics = c("Tetracycline","Doxycycline","Minocycline"),
    genes = c("tet(A)","tet(B)","tet(G)")
  ),
  SXT_SulfaTrim = list(
    antibiotics = c("Sulfamethoxazole","Sulfizoxazole","Trimethoprim_Sulfamethoxazole","Trimethoprim"),
    genes = c("sul1","sul2","sul3","dfrA5","dfrA14","dfrA17")
  ),
  Phenicols = list(
    antibiotics = c("Chloramphenicol"),
    genes = c("floR","cmlA1")
  ),
  Aminoglycosides = list(
    antibiotics = c("Streptomycin","Gentamicin","Tobramycin","Amikacin"),
    genes = c("aadA1","aadA2","aadA5","aadA22",
              "aac(3)-IId","aac(3)-IVa",
              "aph(6)-Id","aph(3'')-Ib","aph(3')-Ia","aph(4)-Ia")
  ),
  Penicillins = list(
    antibiotics = c("Ampicillin"),
    genes = c("blaTEM","blaTEM-1","blaCARB-2")
  )
)

# -----------------------------
# 2) Helper functions
# -----------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
}

read_input <- function(path) {
  ext <- tolower(file_ext(path))
  if (ext == "csv") {
    df <- readr::read_csv(path, show_col_types = FALSE) |> as.data.frame()
  } else if (ext %in% c("xlsx","xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
    df <- readxl::read_excel(path) |> as.data.frame()
  } else {
    stop("Unsupported file type: ", ext, " (use CSV or Excel)")
  }
  names(df) <- trimws(names(df))
  df
}

any1 <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  as.integer(any(x == 1, na.rm = TRUE))
}

as01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  suppressWarnings({
    y <- as.integer(as.character(x))
  })
  y
}

# Cohen's kappa for two binary vectors
kappa_binary <- function(y_true, y_pred) {
  ok <- stats::complete.cases(y_true, y_pred)
  y_true <- y_true[ok]; y_pred <- y_pred[ok]
  if (length(y_true) == 0) return(NA_real_)

  tab <- table(factor(y_true, levels=c(0,1)),
               factor(y_pred, levels=c(0,1)))

  n <- sum(tab)
  if (n == 0) return(NA_real_)

  po <- (tab[1,1] + tab[2,2]) / n
  pe <- ((sum(tab[1,]) / n) * (sum(tab[,1]) / n)) +
        ((sum(tab[2,]) / n) * (sum(tab[,2]) / n))

  if (isTRUE(all.equal(1 - pe, 0))) return(NA_real_)
  (po - pe) / (1 - pe)
}

calc_metrics <- function(y_true, y_pred){
  y_true <- as.integer(y_true); y_pred <- as.integer(y_pred)
  ok <- stats::complete.cases(y_true, y_pred)
  y_true <- y_true[ok]; y_pred <- y_pred[ok]

  TP <- sum(y_true==1 & y_pred==1)
  FP <- sum(y_true==0 & y_pred==1)
  TN <- sum(y_true==0 & y_pred==0)
  FN <- sum(y_true==1 & y_pred==0)
  n  <- TP + FP + TN + FN

  Sens <- ifelse((TP+FN)==0, NA_real_, TP/(TP+FN))
  Spec <- ifelse((TN+FP)==0, NA_real_, TN/(TN+FP))
  PPV  <- ifelse((TP+FP)==0, NA_real_, TP/(TP+FP))
  NPV  <- ifelse((TN+FN)==0, NA_real_, TN/(TN+FN))
  Acc  <- ifelse(n==0, NA_real_, (TP+TN)/n)

  Kap <- kappa_binary(y_true, y_pred)
  Phi <- suppressWarnings(stats::cor(y_true, y_pred, method="pearson"))

  data.frame(
    TP=TP, FP=FP, TN=TN, FN=FN, n=n,
    Sensitivity=Sens, Specificity=Spec,
    PPV=PPV, NPV=NPV,
    Accuracy=Acc,
    Kappa=Kap,
    Phi=Phi
  )
}

make_rulebook_df <- function(rules){
  do.call(rbind, lapply(names(rules), function(cls){
    data.frame(
      Class = cls,
      Phenotype_definition = "Any antibiotic in class == 1",
      Genotype_definition  = "Any gene in class rule present == 1",
      Antibiotics = paste(rules[[cls]]$antibiotics, collapse=", "),
      Genes = paste(rules[[cls]]$genes, collapse=", "),
      stringsAsFactors = FALSE
    )
  }))
}

# -----------------------------
# 3) Load MIC/WGS data
# -----------------------------
stop_if_missing(INPUT_FILE)

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(OUT_PARENT, paste0("Run_", run_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df_mic <- read_input(INPUT_FILE)

if (!ID_COL %in% names(df_mic)) {
  stop("ID column not found: '", ID_COL, "'. Available columns include: ",
       paste(head(names(df_mic), 20), collapse=", "), " ...")
}

df_mic <- df_mic |> dplyr::rename(isolate_id = !!ID_COL)

# Save rulebook
rulebook_df <- make_rulebook_df(class_rules)
write.csv(rulebook_df, file.path(out_dir, "rulebook.csv"), row.names = FALSE)

# -----------------------------
# 4) Optional: parse CGE predicted phenotype
# -----------------------------
cge_bin <- NULL

if (!is.null(CGE_FILE) && nzchar(CGE_FILE)) {
  stop_if_missing(CGE_FILE)

  cge_raw <- readr::read_csv(CGE_FILE, show_col_types = FALSE)

  # حاول نلقط عمود LIMS ID و CGE Predicted Phenotype
  lims_col <- grep("LIMS", names(cge_raw), ignore.case = TRUE, value = TRUE)[1]
  cge_pheno_col <- grep("CGE", names(cge_raw), ignore.case = TRUE, value = TRUE)[1]

  if (is.na(lims_col) || is.na(cge_pheno_col)) {
    stop("Could not find 'LIMS' or 'CGE Predicted Phenotype' columns in CGE file.")
  }

  cge2 <- cge_raw %>%
    dplyr::rename(
      isolate_id = !!lims_col,
      CGE_Pheno_raw = !!cge_pheno_col
    ) %>%
    mutate(
      CGE_Pheno = tolower(ifelse(is.na(CGE_Pheno_raw), "", CGE_Pheno_raw))
    )

  # Map MIC drug column names → text tokens used in CGE_Pheno
  drug_map <- c(
    "Ciprofloxacin"                 = "ciprofloxacin",
    "Levofloxacin"                  = "levofloxacin",
    "NalidixicAcid"                 = "nalidixic acid",
    "Tetracycline"                  = "tetracycline",
    "Doxycycline"                   = "doxycycline",
    "Minocycline"                   = "minocycline",
    "Chloramphenicol"               = "chloramphenicol",
    "Ampicillin"                    = "ampicillin",
    "Amoxicillin_ClavulanicAcid"    = "amoxicillin",
    "Piperacillin_Tazobactam"       = "piperacillin",
    "Ticarcillin/Clavulanic Acid"   = "ticarcillin",
    "Streptomycin"                  = "streptomycin",
    "Gentamicin"                    = "gentamicin",
    "Tobramycin"                    = "tobramycin",
    "Amikacin"                      = "amikacin",
    "Sulfamethoxazole"              = "sulfamethoxazole",
    "Sulfizoxazole"                 = "sulfizoxazole",
    "Trimethoprim_Sulfamethoxazole" = "trimethoprim-sulfamethoxazole",
    "Trimethoprim"                  = "trimethoprim",
    "Cefotaxime"                    = "cefotaxime",
    "Ceftriaxone"                   = "ceftriaxone",
    "Ceftazidime"                   = "ceftazidime",
    "Cefepime"                      = "cefepime",
    "Ceftiofur"                     = "ceftiofur"
  )

  all_ab <- unique(unlist(lapply(class_rules, `[[`, "antibiotics")))

  for (ab in all_ab) {
    token <- drug_map[[ab]]
    if (is.null(token) || is.na(token)) {
      token <- tolower(gsub("_", " ", ab))
    }
    colname <- paste0(ab, "_CGE")

    cge2[[colname]] <- dplyr::case_when(
      cge2$CGE_Pheno %in% c("", "susceptible", "-") ~ 0L,
      stringr::str_detect(cge2$CGE_Pheno, stringr::fixed(token, ignore_case = TRUE)) ~ 1L,
      TRUE ~ 0L
    )
  }

  # Keep only ID + CGE binary columns
  keep_cols <- c("isolate_id", grep("_CGE$", names(cge2), value = TRUE))
  cge_bin <- cge2[, keep_cols, drop = FALSE]

  # Save expanded CGE file for audit
  write.csv(cge_bin, file.path(out_dir, "CGE_binary_calls_per_drug.csv"), row.names = FALSE)
}

# -----------------------------
# 5) Merge MIC/WGS with CGE (if available)
# -----------------------------
if (!is.null(cge_bin)) {
  df_full <- df_mic %>% left_join(cge_bin, by = "isolate_id")
} else {
  df_full <- df_mic
}

# -----------------------------
# 6) Main loop per class
# -----------------------------
results_mic_vs_wgs  <- list()
results_with_cge    <- list()

# Log file
log_path <- file.path(out_dir, "run_log.txt")
sink(log_path)
cat("Run ID:", run_id, "\n")
cat("Input MIC/WGS file:", INPUT_FILE, "\n")
cat("Input CGE file    :", ifelse(is.null(CGE_FILE),"NONE",CGE_FILE), "\n")
cat("Rows (MIC file):", nrow(df_mic), "Columns:", ncol(df_mic), "\n")
cat("Phenotype definition:", PHENO_DEF, "\n")
cat("Genotype definition :", GENO_DEF, "\n")
cat("CGE definition      :", CGE_DEF, "\n\n")
sink()

for (cls in names(class_rules)) {

  ab <- class_rules[[cls]]$antibiotics
  genes <- class_rules[[cls]]$genes

  ab_present   <- intersect(ab, names(df_full))
  genes_present <- intersect(genes, names(df_full))

  # CGE columns corresponding to this class
  ab_cge_cols <- paste0(ab, "_CGE")
  ab_cge_present <- intersect(ab_cge_cols, names(df_full))

  if (length(ab_present) == 0) {
    results_mic_vs_wgs[[cls]] <- data.frame(
      Class = cls,
      note = "No phenotype (MIC) columns found for this class",
      stringsAsFactors = FALSE
    )
    next
  }

  tmp <- df_full %>%
    dplyr::select(isolate_id,
                  dplyr::all_of(ab_present),
                  dplyr::all_of(genes_present),
                  dplyr::all_of(ab_cge_present)) %>%
    dplyr::mutate(across(-isolate_id, as01)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pheno_class = any1(c_across(all_of(ab_present))),
      geno_class  = if (length(genes_present)==0) 0L else any1(c_across(all_of(genes_present))),
      cge_class   = if (length(ab_cge_present)==0) NA_integer_
                    else any1(c_across(all_of(ab_cge_present)))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(isolate_id, pheno_class, geno_class, cge_class)

  # Save per-isolate calls (MIC/WGS/CGE)
  write.csv(tmp, file.path(out_dir, paste0(cls, "_per_isolate_calls_MIC_WGS_CGE.csv")),
            row.names = FALSE)

  # Confusion MIC vs WGS genes (2x2)
  tab_pg <- table(
    Phenotype = factor(tmp$pheno_class, levels=c(0,1)),
    Genotype  = factor(tmp$geno_class,  levels=c(0,1))
  )
  write.csv(as.data.frame.matrix(tab_pg),
            file.path(out_dir, paste0(cls, "_confusion_MIC_vs_WGS_2x2.csv")))

  # Metrics MIC vs WGS genes
  met_pg <- calc_metrics(tmp$pheno_class, tmp$geno_class)

  results_mic_vs_wgs[[cls]] <- cbind(
    data.frame(
      Class = cls,
      Antibiotics_used = paste(ab_present, collapse=", "),
      Genes_used       = paste(genes_present, collapse=", "),
      stringsAsFactors = FALSE
    ),
    met_pg
  )

  # ---- If CGE data available, compute extra metrics ----
  if (!all(is.na(tmp$cge_class))) {

    # MIC vs CGE
    met_pc <- calc_metrics(tmp$pheno_class, tmp$cge_class)
    names(met_pc) <- paste0(names(met_pc), "_MIC_vs_CGE")

    # WGS genes vs CGE
    met_gc <- calc_metrics(tmp$geno_class, tmp$cge_class)
    names(met_gc) <- paste0(names(met_gc), "_GENE_vs_CGE")

    results_with_cge[[cls]] <- cbind(
      data.frame(
        Class = cls,
        Antibiotics_used = paste(ab_present, collapse=", "),
        Genes_used       = paste(genes_present, collapse=", "),
        stringsAsFactors = FALSE
      ),
      met_pg,      # الأساس: MIC vs WGS
      met_pc,      # MIC vs CGE
      met_gc       # GENE vs CGE
    )
  } else {
    results_with_cge[[cls]] <- cbind(
      data.frame(
        Class = cls,
        Antibiotics_used = paste(ab_present, collapse=", "),
        Genes_used       = paste(genes_present, collapse=", "),
        stringsAsFactors = FALSE
      ),
      met_pg
    )
  }
}

# -----------------------------
# 7) Save summaries
# -----------------------------
final_pg   <- dplyr::bind_rows(results_mic_vs_wgs)
final_cge  <- dplyr::bind_rows(results_with_cge)

# الأساس القديم: MIC vs WGS فقط
write.csv(final_pg,
          file.path(out_dir, "class_metrics_MIC_vs_WGS.csv"),
          row.names = FALSE)

# الملف الموسع الذي يشمل CGE أيضاً
write.csv(final_cge,
          file.path(out_dir, "class_metrics_with_CGE.csv"),
          row.names = FALSE)

# Append session info to log
sink(log_path, append = TRUE)
cat("\nR sessionInfo():\n")
print(sessionInfo())
sink()

message("DONE ✅  Outputs saved to: ", out_dir)
print(final_cge)
