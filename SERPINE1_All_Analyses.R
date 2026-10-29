#!/usr/bin/env Rscript

# =============================================================================
# normalize_GSE145989_GSE18995.R
#
# Purpose:
#   Normalize the two GEO-derived, GPL570-annotated gene-expression matrices
#   using the workflow in the supplied geo05.normalize.R, then select the final
#   donor-paired samples and write:
#     inputs/GSE145989_normalized_gene_expression.tsv
#     inputs/GSE18995_normalized_gene_expression.tsv
#
# Input files (edit paths below if needed):
#   GPL570_annotated_output/GSE145989_gene_symbol_expression.tsv
#   GPL570_annotated_output/GSE18995_gene_symbol_expression.tsv
#
# Processing:
#   - Reads annotated gene-symbol matrices; no CEL files are used.
#   - Averages any duplicated gene rows using limma::avereps().
#   - Uses the same scale heuristic as geo05.normalize.R to decide whether to
#     apply log2(x + 1).
#   - Applies limma::normalizeBetweenArrays() to all samples in each dataset
#     separately, before selecting the final paired samples.
#   - Does not merge cohorts or perform cross-cohort batch correction.
#   - Selects samples from the final Table S1 pair map: 53 pairs for GSE145989
#     and 16 pairs for GSE18995.
#
# Important GSE18995 note:
#   The supplied s1(2).txt and s2(2).txt contain 18 Control and 17 Treat GSMs.
#   The final Table S1 pair map retains 16 matched pairs. This script follows
#   the final pair map and therefore excludes GSM469903, GSM469906, and GSM469910
#   from the final GSE18995 paired output. They may still be present in the
#   input matrix and are included in the within-dataset normalization step.
#
# Run from the project root:
#   Rscript normalize_GSE145989_GSE18995.R
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Package 'limma' is required. Install it with BiocManager::install('limma').")
}
library(limma)

# ------------------------------- USER SETTINGS -------------------------------

INPUT_FILES <- c(
  GSE145989 = "GPL570_annotated_output/GSE145989_gene_symbol_expression.tsv",
  GSE18995  = "GPL570_annotated_output/GSE18995_gene_symbol_expression.tsv"
)

OUTPUT_DIR <- "inputs"
REPORT_DIR <- "normalization_reports"
NORMALIZE_ALL_MATRIX_SAMPLES_FIRST <- TRUE

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- FINAL PAIR MAPS ---------------------------------
PAIR_MAPS <- list(
  GSE145989 = do.call(rbind, list(
    data.frame(Donor_ID="Toronto_571", CIT_GSM="GSM4343791", HR_GSM="GSM4343761", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_575", CIT_GSM="GSM4343769", HR_GSM="GSM4343770", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_578", CIT_GSM="GSM4343775", HR_GSM="GSM4343776", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_589", CIT_GSM="GSM4343779", HR_GSM="GSM4343781", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_602", CIT_GSM="GSM4343744", HR_GSM="GSM4343746", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_605", CIT_GSM="GSM4343784", HR_GSM="GSM4343787", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_606", CIT_GSM="GSM4343732", HR_GSM="GSM4343733", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_607", CIT_GSM="GSM4343795", HR_GSM="GSM4343796", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_614", CIT_GSM="GSM4343800", HR_GSM="GSM4343801", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_618", CIT_GSM="GSM4343676", HR_GSM="GSM4343677", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_626", CIT_GSM="GSM4343693", HR_GSM="GSM4343694", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_627", CIT_GSM="GSM4343777", HR_GSM="GSM4343778", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_635", CIT_GSM="GSM4343797", HR_GSM="GSM4343798", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_639", CIT_GSM="GSM4343807", HR_GSM="GSM4343808", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_645", CIT_GSM="GSM4343762", HR_GSM="GSM4343763", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_652", CIT_GSM="GSM4343666", HR_GSM="GSM4343668", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_654", CIT_GSM="GSM4343702", HR_GSM="GSM4343704", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_659", CIT_GSM="GSM4343771", HR_GSM="GSM4343774", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_667", CIT_GSM="GSM4343809", HR_GSM="GSM4343810", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_679", CIT_GSM="GSM4343678", HR_GSM="GSM4343679", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_684", CIT_GSM="GSM4343683", HR_GSM="GSM4343684", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_686", CIT_GSM="GSM4343685", HR_GSM="GSM4343686", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_690", CIT_GSM="GSM4343816", HR_GSM="GSM4343817", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_692", CIT_GSM="GSM4343742", HR_GSM="GSM4343743", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_699", CIT_GSM="GSM4343673", HR_GSM="GSM4343675", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_704", CIT_GSM="GSM4343741", HR_GSM="GSM4343749", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_705", CIT_GSM="GSM4343714", HR_GSM="GSM4343715", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_717", CIT_GSM="GSM4343669", HR_GSM="GSM4343670", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_719", CIT_GSM="GSM4343719", HR_GSM="GSM4343720", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_728", CIT_GSM="GSM4343811", HR_GSM="GSM4343812", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_740", CIT_GSM="GSM4343707", HR_GSM="GSM4343708", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_744", CIT_GSM="GSM4343736", HR_GSM="GSM4343737", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_747", CIT_GSM="GSM4343691", HR_GSM="GSM4343692", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_757", CIT_GSM="GSM4343680", HR_GSM="GSM4343682", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_758", CIT_GSM="GSM4343696", HR_GSM="GSM4343697", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_763", CIT_GSM="GSM4343700", HR_GSM="GSM4343701", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_772", CIT_GSM="GSM4343705", HR_GSM="GSM4343706", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_780", CIT_GSM="GSM4343671", HR_GSM="GSM4343672", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_785", CIT_GSM="GSM4343757", HR_GSM="GSM4343758", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_786", CIT_GSM="GSM4343750", HR_GSM="GSM4343751", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_793", CIT_GSM="GSM4343709", HR_GSM="GSM4343711", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_796", CIT_GSM="GSM4343734", HR_GSM="GSM4343735", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_824", CIT_GSM="GSM4343712", HR_GSM="GSM4343713", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_828", CIT_GSM="GSM4343729", HR_GSM="GSM4343730", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_832", CIT_GSM="GSM4343753", HR_GSM="GSM4343754", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_833", CIT_GSM="GSM4343782", HR_GSM="GSM4343783", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_834", CIT_GSM="GSM4343739", HR_GSM="GSM4343740", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_853", CIT_GSM="GSM4343802", HR_GSM="GSM4343803", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_855", CIT_GSM="GSM4343721", HR_GSM="GSM4343690", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_857", CIT_GSM="GSM4343766", HR_GSM="GSM4343768", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_858", CIT_GSM="GSM4343727", HR_GSM="GSM4343728", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_864", CIT_GSM="GSM4343793", HR_GSM="GSM4343794", stringsAsFactors=FALSE),
    data.frame(Donor_ID="Toronto_898", CIT_GSM="GSM4343725", HR_GSM="GSM4343726", stringsAsFactors=FALSE)
  )),
  GSE18995 = do.call(rbind, list(
    data.frame(Donor_ID="DBD_patient_01", CIT_GSM="GSM469912", HR_GSM="GSM469918", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DBD_patient_02", CIT_GSM="GSM469913", HR_GSM="GSM469919", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DBD_patient_03", CIT_GSM="GSM469914", HR_GSM="GSM469920", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DBD_patient_04", CIT_GSM="GSM469915", HR_GSM="GSM469921", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DBD_patient_05", CIT_GSM="GSM469916", HR_GSM="GSM469922", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DBD_patient_06", CIT_GSM="GSM469917", HR_GSM="GSM469923", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DBD_patient_07", CIT_GSM="GSM469924", HR_GSM="GSM469930", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DBD_patient_08", CIT_GSM="GSM469925", HR_GSM="GSM469931", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DBD_patient_09", CIT_GSM="GSM469926", HR_GSM="GSM469932", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DBD_patient_10", CIT_GSM="GSM469927", HR_GSM="GSM469933", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DBD_patient_11", CIT_GSM="GSM469928", HR_GSM="GSM469934", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DBD_patient_12", CIT_GSM="GSM469929", HR_GSM="GSM469935", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DCD_patient_01", CIT_GSM="GSM469901", HR_GSM="GSM469907", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DCD_patient_02", CIT_GSM="GSM469902", HR_GSM="GSM469908", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DCD_patient_04", CIT_GSM="GSM469904", HR_GSM="GSM469909", stringsAsFactors=FALSE),
    data.frame(Donor_ID="DCD_patient_06", CIT_GSM="GSM469905", HR_GSM="GSM469911", stringsAsFactors=FALSE)
  ))
)

for (gse in names(PAIR_MAPS)) {
  row.names(PAIR_MAPS[[gse]]) <- NULL
}
stopifnot(nrow(PAIR_MAPS$GSE145989) == 53L)
stopifnot(nrow(PAIR_MAPS$GSE18995) == 16L)

# ------------------------------- FUNCTIONS ------------------------------------

read_gene_matrix <- function(path, gse) {
  if (!file.exists(path)) {
    stop(
      "Input matrix not found for ", gse, ": ", path,
      "\nRun the GPL570 matrix annotation script first, or edit INPUT_FILES."
    )
  }
  
  dat <- read.delim(
    path, header = TRUE, sep = "\t", check.names = FALSE,
    quote = "", comment.char = "", stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", "null")
  )
  if (ncol(dat) < 2L) stop(gse, ": input must have a gene-ID column plus sample columns.")
  
  gene_id <- trimws(as.character(dat[[1]]))
  if (anyNA(gene_id) || any(!nzchar(gene_id))) {
    stop(gse, ": blank or missing gene IDs found in the first column.")
  }
  
  value_df <- dat[, -1, drop = FALSE]
  value_df[] <- lapply(value_df, function(x) suppressWarnings(as.numeric(as.character(x))))
  expr <- as.matrix(value_df)
  rownames(expr) <- gene_id
  
  if (all(is.na(expr))) stop(gse, ": no numeric expression values found.")
  if (anyDuplicated(colnames(expr))) stop(gse, ": duplicated sample-column names found.")
  
  # Same behavior as geo05.normalize.R: average duplicate row IDs if any remain.
  expr <- limma::avereps(expr)
  expr
}

normalize_matrix <- function(expr, gse) {
  qx <- as.numeric(stats::quantile(
    expr, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1.0),
    na.rm = TRUE, names = FALSE
  ))
  log_needed <- (qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0)
  
  if (log_needed) {
    # Match the supplied geo05.normalize.R behavior.
    expr[expr < 0] <- 0
    expr <- log2(expr + 1)
    message(gse, ": applied log2(x + 1) based on the expression-scale heuristic.")
  } else {
    message(gse, ": no additional log2 transformation applied.")
  }
  
  if (NORMALIZE_ALL_MATRIX_SAMPLES_FIRST) {
    # Match geo05.normalize.R: normalize all available matrix columns first,
    # then subset the final donor-paired samples.
    normalized <- limma::normalizeBetweenArrays(expr)
  } else {
    normalized <- expr
  }
  
  list(matrix = normalized, log_transform_applied = log_needed, quantiles = qx)
}

make_paired_output <- function(norm_expr, pair_map, gse) {
  cit <- pair_map$CIT_GSM
  hr  <- pair_map$HR_GSM
  required <- c(cit, hr)
  missing <- setdiff(required, colnames(norm_expr))
  if (length(missing)) {
    stop(
      gse, ": input matrix is missing required GSM(s): ",
      paste(missing, collapse = ", "),
      "\nCheck that the correct GEO matrix was annotated and that its columns use GSM accessions."
    )
  }
  
  if (NORMALIZE_ALL_MATRIX_SAMPLES_FIRST) {
    data_norm <- norm_expr
  } else {
    # If the optional setting is FALSE, normalize only the final selected samples.
    data_norm <- limma::normalizeBetweenArrays(norm_expr[, required, drop = FALSE])
  }
  
  # Preserve the downstream input layout: all CIT/Control columns first,
  # followed by all HR/Treat columns.
  paired <- cbind(
    data_norm[, cit, drop = FALSE],
    data_norm[, hr, drop = FALSE]
  )
  colnames(paired) <- c(paste0(cit, "_Control"), paste0(hr, "_Treat"))
  paired
}

# ---------------------------------- RUN ---------------------------------------

summary_rows <- list()

for (gse in names(INPUT_FILES)) {
  message("\n========== ", gse, " ==========")
  expr <- read_gene_matrix(INPUT_FILES[[gse]], gse)
  pair_map <- PAIR_MAPS[[gse]]
  
  message("Input matrix: ", nrow(expr), " genes x ", ncol(expr), " samples")
  norm_result <- normalize_matrix(expr, gse)
  norm_expr <- norm_result$matrix
  
  paired <- make_paired_output(norm_expr, pair_map, gse)
  out_file <- file.path(OUTPUT_DIR, paste0(gse, "_normalized_gene_expression.tsv"))
  
  # Write a conventional tab-delimited matrix with gene IDs in column 1 and
  # sample IDs in the header; compatible with downstream read.delim(...).
  out_df <- data.frame(id = rownames(paired), paired, check.names = FALSE)
  write.table(
    out_df, file = out_file, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
  
  sample_map <- data.frame(
    Donor_ID = pair_map$Donor_ID,
    CIT_GSM = pair_map$CIT_GSM,
    HR_GSM = pair_map$HR_GSM,
    CIT_output_column = paste0(pair_map$CIT_GSM, "_Control"),
    HR_output_column = paste0(pair_map$HR_GSM, "_Treat"),
    stringsAsFactors = FALSE
  )
  write.csv(
    sample_map,
    file = file.path(REPORT_DIR, paste0(gse, "_sample_map.csv")),
    row.names = FALSE
  )
  
  report <- c(
    paste0("Dataset: ", gse),
    paste0("Input file: ", INPUT_FILES[[gse]]),
    paste0("Input genes after duplicate averaging: ", nrow(expr)),
    paste0("Input matrix samples: ", ncol(expr)),
    paste0("Final donor pairs: ", nrow(pair_map)),
    paste0("Final output samples: ", ncol(paired)),
    paste0("Scale heuristic applied log2(x + 1): ", norm_result$log_transform_applied),
    "Between-array normalization: limma::normalizeBetweenArrays",
    paste0("Normalize all matrix samples before selecting final pairs: ",
           NORMALIZE_ALL_MATRIX_SAMPLES_FIRST),
    "CEL files used for expression calculation: no",
    paste0("R version: ", R.version.string),
    paste0("limma version: ", as.character(utils::packageVersion("limma"))),
    paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  )
  writeLines(report, file.path(REPORT_DIR, paste0(gse, "_normalization_report.txt")))
  
  summary_rows[[gse]] <- data.frame(
    Dataset = gse,
    Input_genes = nrow(expr),
    Input_samples = ncol(expr),
    Final_pairs = nrow(pair_map),
    Final_samples = ncol(paired),
    Log2_transform_applied = norm_result$log_transform_applied,
    Output_file = out_file,
    stringsAsFactors = FALSE
  )
  
  message("Wrote: ", out_file)
  message("Output dimensions: ", nrow(paired), " genes x ", ncol(paired), " samples.")
}

summary_table <- do.call(rbind, summary_rows)
write.csv(
  summary_table,
  file = file.path(REPORT_DIR, "normalization_summary.csv"),
  row.names = FALSE
)

message("\nCompleted both datasets independently.")
message("GSE145989: 53 Control + 53 Treat samples.")
message("GSE18995: 16 Control + 16 Treat samples.")
message("No CEL files were used to calculate expression values.")

# =============================================================================
# SERPINE1 / PAI-1 LIRI MANUSCRIPT: SINGLE-FILE ENGLISH ANALYSIS WORKFLOW
# Figure order: 1-8, followed by supplementary integration diagnostics.
# Final Figure 5: A-F trajectory; G-J CellChat; K paired expression proxies.
# No standalone sample-level Figure 5 program-comparison test is included.
# =============================================================================
# HOW TO USE
# 1. Save this file and edit PROJECT_ROOT and CONFIG input paths below.
# 2. Source the complete file once. This defines all stages without running them.
# 3. list_stages() lists stage names. Run a stage with run_stage("figure1_wgcna").
# 4. To run several stages, set STAGES_TO_RUN to their names, in order.
# Command line: Rscript SERPINE1_All_Analyses_No_Original_5G.R --list
# Command line: Rscript SERPINE1_All_Analyses_No_Original_5G.R --check
# Command line: Rscript SERPINE1_All_Analyses_No_Original_5G.R figure1_wgcna
#
# INPUTS MUST BE THE ACTUAL FINAL DATA AND METADATA, NOT PUBLISHED P VALUES.
# This is code only: no input matrices, computed results, or figures are bundled.
# The actual 113-model list, test matrix, atlas/endothelial checkpoints,
# Monocle session, vendor report and in-vitro measurements must be supplied.
# The original final annotations and trajectory construction settings cannot be
# reconstructed from missing source files. Saved final objects are required by
# the corresponding continuation stages. Figure 8 uses a schema-driven summary
# module; its real measurements and prespecified contrast table are required.
# Source-derived Figure 2 filenames retain historical subpanel names; match their
# contents to the final manuscript rather than interpreting filenames as evidence.
#
# VALIDATION: static delimiter/string checks only in the preparation environment.
# R was unavailable there. Parse and execute this file in the original analysis
# environment before claiming complete reproduction or uploading a final release.
# ChatGPT (OpenAI) assisted with code editing and consolidation. Record the
# authors' actual execution/review procedures and package versions on release.
# =============================================================================

PROJECT_ROOT <- "."  # Replace with your analysis-project directory.
STAGES_TO_RUN <- character(0)  # Example: c("figure1_wgcna","figure1_blue_enrichment")
REPO_ROOT <- normalizePath(PROJECT_ROOT,winslash="/",mustWork=TRUE)

# Find this single source file for input-version hashing; no helper R files
# are needed. If pasting into a console, set SCRIPT_FILE to its saved path.
SCRIPT_FILE <- ""
.file_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(length(.file_arg)) SCRIPT_FILE <- normalizePath(sub("^--file=","",.file_arg[1]),mustWork=FALSE)
if(!nzchar(SCRIPT_FILE)) {
 .frames <- sys.frames()
 for(.frame in rev(.frames)) if(!is.null(.frame$ofile)) {
  SCRIPT_FILE <- normalizePath(.frame$ofile,mustWork=FALSE);break
 }
}

# Required sample crosswalk schemas:
# metadata/GSE145989_pairs.csv: 53 rows, donor_id,CIT_GSM,R2H_GSM,
#   CIT_sample_ID,R2H_sample_ID; optional stage_review_required,Source_note.
# metadata/GSE18995_pairs.csv: 16 rows with the same columns.
# metadata/single_cell_sample_map.csv: Sample,Donor,Condition (12 samples).
# ML matrices: sample rows, 14 candidate-gene columns and Type; CIT/HR labels.
# ml_methods_113.tsv: the actual 113 unique fitted combinations in a Model column.
# Figure 8: Endpoint,Group,Replicate,Value, optional Field; contrast file:
# inputs/figure8_contrasts.csv with Endpoint,Group1,Group2,Test.

# Edit paths only after locating the actual final source inputs.
CONFIG <- list(
 output_directory="outputs",
 endothelial_state_column="Endothelial_Annotation",
 endothelial_label="Endothelial cell",
 stress_state="Stress-responsive capillary EC",
 nk_label="NK cell",monocyte_label="Monocyte",macrophage_label="Macrophage",
 inputs=list(
  ml_train="inputs/ml_train.tsv",
  ml_test="inputs/ml_test.tsv",
  ml_methods="inputs/ml_methods_113.tsv",
  ml_run_pointer="outputs/Figure1/ML/selected_run_pointer.rds",
  bulk_train_expression="inputs/GSE145989_normalized_gene_expression.tsv",
  bulk_paired_results="inputs/bulk_paired_limma.tsv",
  bulk_paired_metadata="inputs/bulk_paired_metadata.tsv",
  wgcna_checkpoint="outputs/Figure1/WGCNA/wgcna_checkpoint.rds",
  sce_rdata="inputs/final_annotated_atlas.RData",
  tenx_directory="inputs/10X",
  figure3_paired_tables="outputs/Figure3/Paired/Tables",
  figure3_signed_rank="outputs/Figure3/Paired/Tables/06_paired_signed_rank_for_GSEA.csv",
  single_cell_qc_directory="inputs/single_cell_QC",
  endothelial_bundle="inputs/Figure5_analysis_objects_PREMIUM.rds",
  network_master="outputs/Shared/Figure6_MASTER_DATA_ONLY.rds",
  trajectory_session="inputs/trajectory_session.rds",
  trajectory_data="inputs/trajectory_data.rds",
  perturbation_bundle="outputs/Figure6/Compute/05_OBJECTS/SERPINE1_virtualKO_OPTIMIZED_ANALYSIS_BUNDLE.rds",
  metabolic_bundle="outputs/Figure7/Metabolic/03_objects/Figure7_STAGE1_SC_METABOLIC_DISCOVERY_BUNDLE.rds",
  metabolic_panels="outputs/Figure7/Metabolic/Figure7_STAGE1_A_C_FINAL_LOGIC_BUNDLE.rds",
  lipidomics_directory="inputs/vendor_lipidomics",
  in_vitro_measurements="inputs/figure8_measurements.csv"
 ))

# Shared input/output checks
# Shared repository utilities. No data or statistical results are fabricated.
require_packages <- function(pkgs) {
 missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
 if(length(missing)) stop("Install required packages first: ",paste(missing,collapse=", "))
}
repo_input <- function(key) {
 if(!key %in% names(CONFIG$inputs)) stop("Undefined input key: ",key)
 p <- CONFIG$inputs[[key]]
 if(!grepl("^(/|[A-Za-z]:)",p)) p <- file.path(REPO_ROOT,p)
 if(!file.exists(p)) stop("Required input is missing: ",key," -> ",p,
 ". See docs/INPUTS.md; supply the actual final input, not a placeholder result.")
 p
}
repo_output <- function(relative) {
 p <- file.path(REPO_ROOT,CONFIG$output_directory,relative)
 dir.create(if(grepl("\\.[A-Za-z0-9]+$",basename(p))) dirname(p) else p,
            recursive=TRUE,showWarnings=FALSE)
 p
}
restore_named_list <- function(path) {
 x<-readRDS(path)
 if(!is.list(x)||is.null(names(x))) stop("Expected a named RDS list: ",path)
 list2env(x,envir=.GlobalEnv);invisible(x)
}
load_sce <- function() {
 p<-repo_input("sce_rdata")
 if(grepl("\\.rds$",p,ignore.case=TRUE)) return(readRDS(p))
 e<-new.env();load(p,envir=e)
 hits<-ls(e)[vapply(as.list(e),inherits,logical(1),what="Seurat")]
 if(length(hits)!=1) stop("The saved atlas must contain exactly one Seurat object.")
 e[[hits]]
}
save_publication_plot <- function(p,stem,width=8,height=6) {
 require_packages("ggplot2")
 ggplot2::ggsave(paste0(stem,".pdf"),p,width=width,height=height)
 ggplot2::ggsave(paste0(stem,".png"),p,width=width,height=height,dpi=300)
}
repo_ml_run <- function() {
 p<-repo_input("ml_run_pointer");x<-readRDS(p)
 if(!dir.exists(x$RUN_DIR)) stop("Model-run pointer is stale; rerun the model stage or configure a valid run.")
 x$RUN_DIR
}
validate_ml_inputs <- function(train,test,methods) {
 a<-read.delim(train,row.names=1,check.names=FALSE);b<-read.delim(test,row.names=1,check.names=FALSE)
 ga<-setdiff(names(a),"Type");gb<-setdiff(names(b),"Type")
 if(length(ga)!=14L||!setequal(ga,gb)) stop("Both ML matrices must contain the same 14 candidate genes and Type.")
 m<-read.delim(methods,check.names=FALSE)
 if(!"Model" %in% names(m)||anyDuplicated(m$Model)||nrow(m)!=113L)
   stop("Supply the actual 113 unique model combinations in a tab-separated Model column.")
 invisible(TRUE)
}

# Embedded donor-safe ML helpers
load_ml_helpers <- function() {
 eval(quote({
# R/ml_helpers.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

# Enhanced 2026-10-21: donor-safe RF importance selection and fixed-nt
# plsRglm (experimental), SVM oriented decision scores without row-CV Platt fit.
# 113 proposed methods can be attempted; failures are explicitly reported.
# NOTE: RF importance selection is not identical to the original var.select.
# =============================================================================
# Donor-aware model-selection and classification backends.
# ALL ROC selection metrics are computed on held-out DONORS.
# For inner tuning of glmnet/glmBoost/GBM/XGBoost, held-out folds are DONOR groups.
# All supervised selectors are refit INSIDE each outer fold; no test labels used.
# Unsupported algorithm families are SKIPPED, never silently relabelled.
# =============================================================================

make_donor_folds <- function(donor, k = 5L, seed = 123456L) {
  donor <- as.character(donor)
  if (anyNA(donor) || any(!nzchar(donor))) stop("Invalid donor IDs")
  u <- sort(unique(donor)); k <- as.integer(k)
  if (length(u) < k || k < 2L) stop("Number of donors must be >= 2 and >= fold count")
  set.seed(seed)
  u <- sample(u, length(u))
  mapping <- stats::setNames(rep(seq_len(k), length.out = length(u)), u)
  f <- unname(mapping[donor])
  if (anyNA(f) || any(vapply(split(f, donor), function(x) length(unique(x)) != 1L, logical(1)))) {
    stop("Donor-fold assignment failed")
  }
  as.integer(f)
}

# Independently supplied verified GEO GSM<->donor crosswalk.
read_verified_pairs <- function(path, sample_ids, y, expected_pairs) {
  if (!file.exists(path)) stop("Missing verified donor mapping: ", path,
    "\nProvide the already verified pair table from your earlier CEL/GSM audit; ",
    "do not construct donor pairs from GSM order.")
  tab <- if (grepl("\\.csv$", path, ignore.case=TRUE)) {
    utils::read.csv(path, check.names=FALSE, stringsAsFactors=FALSE,
                    fileEncoding="UTF-8-BOM")
  } else {
    utils::read.delim(path, check.names=FALSE, stringsAsFactors=FALSE,
                      fileEncoding="UTF-8-BOM")
  }
  if (all(c("donor_id", "CIT_GSM", "R2H_GSM") %in% names(tab))) {
    long <- data.frame(
      GSM = c(as.character(tab$CIT_GSM), as.character(tab$R2H_GSM)),
      Donor = rep(as.character(tab$donor_id), 2L),
      Type = rep(c(0L, 1L), each=nrow(tab)), stringsAsFactors=FALSE)
  } else if (all(c("sample", "donor", "condition") %in% names(tab))) {
    condition <- tolower(as.character(tab$condition))
    type <- ifelse(condition %in% c("cit","control","0"), 0L,
                   ifelse(condition %in% c("r2h","treat","hr","2h","1"), 1L, NA_integer_))
    long <- data.frame(
      GSM = sub(".*(GSM[0-9]+).*", "\\1", as.character(tab$sample)),
      Donor = as.character(tab$donor), Type=type, stringsAsFactors=FALSE)
  } else stop("Invalid crosswalk columns. Require donor_id,CIT_GSM,R2H_GSM ",
              "or sample,donor,condition")
  if (anyNA(long$Type) || anyNA(long$Donor) || any(!nzchar(long$Donor)) ||
      anyDuplicated(long$GSM) || any(!grepl("^GSM[0-9]+$", long$GSM))) {
    stop("Crosswalk has invalid, duplicate or missing GSM/donor/condition entries")
  }
  gsm <- sub(".*(GSM[0-9]+).*", "\\1", sample_ids)
  if (any(!grepl("^GSM[0-9]+$", gsm))) stop("Sample IDs do not encode GSM numbers")
  matched <- match(gsm, long$GSM)
  if (anyNA(matched)) {
    stop("Unmapped GSM: ", paste(unique(gsm[is.na(matched)]), collapse=", "))
  }
  donor <- long$Donor[matched]
  if (any(as.integer(y) != long$Type[matched])) stop("Type disagrees with verified donor mapping")
  tab2 <- table(donor, factor(as.integer(y), levels=0:1))
  if (!all(tab2==1L) || nrow(tab2)!=expected_pairs || length(y)!=2L*expected_pairs) {
    stop("Expected ",expected_pairs," complete donors; observed ",nrow(tab2))
  }
  # Extra title checks for unexpected label reversals when present in sample ID.
  labels <- ifelse(grepl("_(Control|CIT)$", sample_ids), 0L,
             ifelse(grepl("_(Treat|R2H|HR)$", sample_ids), 1L, NA_integer_))
  if (anyNA(labels) || any(labels != y)) stop("Sample suffix and Type disagree")
  donor
}

# Within EACH outer fold, fit the center/scale on that fold's training donors.
fit_group_scaler <- function(X) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  if (!ncol(X) || any(!is.finite(X))) stop("Invalid expression matrix")
  mu <- colMeans(X); sdv <- apply(X, 2L, stats::sd)
  keep <- is.finite(sdv) & sdv > 1e-10
  if (!any(keep)) stop("No nonconstant predictors")
  list(center=mu[keep], scale=sdv[keep], kept=names(mu)[keep], dropped=names(mu)[!keep])
}
transform_by_training <- function(X, scaler) {
  X <- as.matrix(X[,scaler$kept,drop=FALSE]); storage.mode(X) <- "double"
  z <- sweep(sweep(X,2L,scaler$center,"-"),2L,scaler$scale,"/")
  if (any(!is.finite(z))) stop("Non-finite values after applying TRAINING-only scaler")
  as.data.frame(z, check.names=FALSE)
}

parse_ml_method <- function(s) {
  s <- gsub("[[:space:]]", "", as.character(s))
  ss <- strsplit(s,"\\+",fixed=FALSE)[[1]]
  if (length(ss)==1L) ss <- c("simple",ss)
  if (length(ss)!=2L) stop("Unsupported combination string: ",s)
  parse_piece <- function(x) {
    a <- sub("\\[.*$", "", x)
    param <- if (grepl("\\[",x)) sub("^[^\\[]*\\[([^]]*)\\]$","\\1",x) else ""
    list(method=a, param=param)
  }
  list(label=s,selector=parse_piece(ss[1]),learner=parse_piece(ss[2]))
}

# Algorithm availability: no groups can be preserved inside standard caret CV,
# naive plsRglm CV, or RF var.select OOB. Do not pretend those are paired.
selector_supported <- function(nm) nm %in% c(
  "simple","Lasso","Ridge","Enet","Stepglm","glmBoost","SVM","LDA",
  "NaiveBayes","GBM","XGBoost","RF")
learner_supported <- function(nm) nm %in% c(
  "Lasso","Ridge","Enet","Stepglm","glmBoost","SVM","LDA",
  "NaiveBayes","RF","GBM","XGBoost","plsRglm")

# No artificial feature-number threshold; report actual number after refit.
# RF-selector is implemented as a FIXED within-training-fold importance selector,
# not the original var.select. plsRglm fixed nt<=2 and may fail in some R versions.
parse_alpha <- function(name,param) {
  if (name=="Lasso") return(1)
  if (name=="Ridge") return(0)
  a <- suppressWarnings(as.numeric(sub("^alpha=", "",param)))
  if (is.na(a) || a<0 || a>1) stop("Bad Enet alpha: ",param)
  a
}

make_inner_folds <- function(donor,k=5L,seed=123456L) {
  make_donor_folds(donor,k=min(as.integer(k),length(unique(donor))),seed=seed)
}

fit_boost_grouped <- function(X,y,donor,seed,max_mstop=200L,k=5L) {
  if (!requireNamespace("mboost",quietly=TRUE)) stop("Package mboost missing")
  dat <- data.frame(Type=factor(ifelse(y==1,"R2H","CIT"),levels=c("CIT","R2H")),
                    X,check.names=FALSE)
  set.seed(seed)
  fit <- mboost::glmboost(Type~.,data=dat,family=mboost::Binomial(),
                         control=mboost::boost_control(mstop=max_mstop,nu=0.1))
  ids <- make_inner_folds(donor,k=k,seed=seed+101L)
  foldweights <- vapply(seq_len(max(ids)),function(v) as.numeric(ids!=v),numeric(length(ids)))
  cv <- mboost::cvrisk(fit,folds=foldweights,papply=lapply)
  best <- as.integer(mboost::mstop(cv))[1]
  if (!is.finite(best) || best<1) stop("glmBoost has invalid tuning result")
  list(fit=fit[best],mstop=best)
}

fit_glmnet_grouped <- function(X,y,donor,alpha,seed,k=5L) {
  if (!requireNamespace("glmnet",quietly=TRUE)) stop("Package glmnet missing")
  set.seed(seed)
  glmnet::cv.glmnet(as.matrix(X),as.integer(y),family="binomial",alpha=alpha,
                    foldid=make_inner_folds(donor,k=k,seed=seed+21L),
                    nfolds=min(k,length(unique(donor))),standardize=FALSE,
                    type.measure="deviance")
}

stepglm_fit <- function(X,y,direction="both") {
  # FIX: a dot upper scope in stats::step() expands relative to the
  # intercept-only starting formula, leaving no candidate variables.
  # Explicitly enumerate all candidate genes as the upper model scope.
  dat <- data.frame(Type=as.integer(y), X, check.names=FALSE)
  if (!direction %in% c("both","forward","backward"))
    stop("Unsupported Stepglm direction: ", direction)
  genes <- colnames(X)
  if (!length(genes)) stop("SKIP: Stepglm has no candidate genes")
  # Use backticks to make the upper model safe for non-syntactic gene names.
  term_names <- paste0("`", gsub("`", "", genes, fixed=TRUE), "`")
  full_form <- stats::reformulate(term_names, response="Type")
  if (direction=="forward") {
    start <- stats::glm(Type ~ 1, data=dat, family=stats::binomial())
    fit <- stats::step(start,
                       scope=list(lower=~1, upper=full_form),
                       direction="forward", trace=0)
  } else {
    start <- stats::glm(full_form, data=dat, family=stats::binomial())
    fit <- stats::step(start, direction=direction, trace=0)
  }
  kept <- setdiff(names(stats::coef(fit))[is.finite(stats::coef(fit))], "(Intercept)")
  if (length(kept)==0L)
    stop("SKIP: Stepglm[", direction, "] selected intercept only; no valid gene-based classifier")
  fit
}

fit_gbm_grouped <- function(X,y,donor,seed,max_trees=700L,
                            grid=c(50L,100L,200L,400L,700L),k=5L) {
  if (!requireNamespace("gbm",quietly=TRUE)) stop("Package gbm missing")
  k <- min(as.integer(k),length(unique(donor)))
  fold <- make_inner_folds(donor,k=k,seed=seed+73L)
  grid <- sort(unique(as.integer(grid[grid<=max_trees])))
  if (length(grid)<1) stop("GBM tree grid is empty")
  losses <- matrix(NA_real_, nrow=k,ncol=length(grid))
  for (i in seq_len(k)) {
    a <- which(fold!=i); b <- which(fold==i)
    set.seed(seed+i+2000L)
    mod <- gbm::gbm(Type~.,data=data.frame(Type=as.integer(y[a]),X[a,,drop=FALSE],check.names=FALSE),
                    distribution="bernoulli",n.trees=max_trees,
                    interaction.depth=min(3L,ncol(X)),n.minobsinnode=5L,
                    shrinkage=0.01,bag.fraction=0.75,cv.folds=0,n.cores=1,verbose=FALSE)
    for (j in seq_along(grid)) {
      pr <- as.numeric(stats::predict(mod,X[b,,drop=FALSE],n.trees=grid[j],type="response"))
      pr <- pmax(1e-8,pmin(1-1e-8,pr))
      losses[i,j] <- mean(-y[b]*log(pr)-(1-y[b])*log1p(-pr))
    }
  }
  risks <- colMeans(losses)
  if (any(!is.finite(risks))) stop("Non-finite grouped GBM tuning risks")
  best <- grid[which.min(risks)]
  set.seed(seed+3001L)
  mod <- gbm::gbm(Type~.,data=data.frame(Type=as.integer(y),X,check.names=FALSE),
                  distribution="bernoulli",n.trees=best,
                  interaction.depth=min(3L,ncol(X)),n.minobsinnode=5L,
                  shrinkage=0.01,bag.fraction=0.75,cv.folds=0,n.cores=1,verbose=FALSE)
  list(fit=mod,ntrees=best)
}

fit_xgb_grouped <- function(X,y,donor,seed,max_rounds=100L,k=5L) {
  if (!requireNamespace("xgboost",quietly=TRUE)) stop("Package xgboost missing")
  folds <- make_inner_folds(donor,k=k,seed=seed+31L)
  # Explicit donor-grouped fold indices; NEVER xgboost's default row folds.
  val_indices <- lapply(seq_len(max(folds)),function(i) which(folds==i))
  dtrain <- xgboost::xgb.DMatrix(as.matrix(X),label=as.numeric(y))
  set.seed(seed)
  cv <- xgboost::xgb.cv(params=list(objective="binary:logistic",eval_metric="logloss",
                                max_depth=2L,eta=0.1,nthread=1),
                         data=dtrain,nrounds=max_rounds,folds=val_indices,
                         verbose=0,early_stopping_rounds=10)
  best <- if (!is.null(cv$best_iteration)) as.integer(cv$best_iteration) else
            as.integer(which.min(cv$evaluation_log[["test_logloss_mean"]]))
  if (!is.finite(best) || best<1) stop("XGBoost grouped CV failed")
  set.seed(seed+1L)
  xgboost::xgb.train(params=list(objective="binary:logistic",max_depth=2L,
                                 eta=0.1,nthread=1),data=dtrain,nrounds=best,verbose=0)
}

select_in_train <- function(spec,X,y,donor,seed) {
  nm <- spec$method; param <- spec$param
  if (!selector_supported(nm)) stop("Unsupported selector (not group-audited): ",nm)
  if (nm %in% c("simple","SVM","LDA","NaiveBayes","XGBoost")) return(colnames(X))
  if (nm %in% c("Lasso","Ridge","Enet")) {
    alpha <- parse_alpha(nm,param)
    mod <- fit_glmnet_grouped(X,y,donor,alpha=alpha,seed=seed)
    cmat <- as.matrix(stats::coef(mod,s="lambda.min"))
    return(intersect(colnames(X), rownames(cmat)[cmat[,1]!=0]))
  }
  if (nm=="Stepglm") {
    mod <- stepglm_fit(X,y,direction=param)
    return(intersect(colnames(X),names(stats::coef(mod))[is.finite(stats::coef(mod))]))
  }
  if (nm=="glmBoost") {
    mod <- fit_boost_grouped(X,y,donor,seed=seed)
    b <- unlist(stats::coef(mod$fit),use.names=TRUE)
    names(b) <- gsub("^`|`$","",names(b))
    return(intersect(colnames(X),names(b)[is.finite(b)&abs(b)>1e-10]))
  }
  if (nm=="RF") {
    # RF importance calculated ONLY from the current outer-training donors.
    # This is a fixed, prospective selector distinct from rfsrc::var.select.
    if (!requireNamespace("randomForestSRC",quietly=TRUE)) stop("Package randomForestSRC missing")
    set.seed(seed)
    mod <- randomForestSRC::rfsrc(
      Type ~ ., data=data.frame(Type=factor(y,levels=0:1),X,check.names=FALSE),
      ntree=500L, nodesize=5L, importance="permute", proximity=FALSE, forest=TRUE)
    imp <- mod$importance
    if (is.matrix(imp)) imp <- imp[,ncol(imp)]
    if (is.null(names(imp)) && length(imp)==ncol(X)) names(imp) <- colnames(X)
    imp <- as.numeric(imp[colnames(X)])
    names(imp) <- colnames(X)
    good <- which(is.finite(imp) & imp>0)
    # Fixed minimum top-five fallback, documented in output. Not based on holdout labels.
    if (length(good)<2L) good <- order(imp,decreasing=TRUE,na.last=NA)[seq_len(min(5L,sum(is.finite(imp))))]
    if (!length(good)) stop("RF importance selector has no finite importances")
    return(names(imp)[good])
  }
  if (nm=="GBM") {
    mod <- fit_gbm_grouped(X,y,donor,seed=seed,max_trees=300L,
                           grid=c(50L,100L,200L,300L),k=4L)
    imp <- gbm::summary.gbm(mod$fit,n.trees=mod$ntrees,plotit=FALSE)
    return(intersect(colnames(X),as.character(imp$var[imp$rel.inf>0])))
  }
  stop("Selector is not implemented: ",nm)
}

fit_classifier <- function(spec,X,y,donor,seed) {
  nm <- spec$method; param <- spec$param
  if (!learner_supported(nm)) stop("Unsupported learner (not group-audited): ",nm)
  if (nm %in% c("Lasso","Ridge","Enet")) {
    m <- fit_glmnet_grouped(X,y,donor,alpha=parse_alpha(nm,param),seed=seed)
    return(list(method=nm,object=m))
  }
  if (nm=="Stepglm") return(list(method=nm,object=stepglm_fit(X,y,param)))
  if (nm=="glmBoost") {
    m <- fit_boost_grouped(X,y,donor,seed=seed)
    return(list(method=nm,object=m$fit,mstop=m$mstop))
  }
  if (nm=="SVM") {
    if (!requireNamespace("e1071",quietly=TRUE)) stop("Package e1071 missing")
    set.seed(seed)
    m <- e1071::svm(x=as.matrix(X),y=factor(y,levels=0:1),
                    probability=FALSE,scale=FALSE)
    trp <- stats::predict(m,newdata=as.matrix(X),decision.values=TRUE)
    mar <- as.numeric(attr(trp,"decision.values"))
    if (length(mar)!=length(y) || any(!is.finite(mar))) stop("SVM decision values invalid")
    direction <- if (mean((mar>0)==(as.character(trp)=="1")) >= 0.5) 1 else -1
    return(list(method=nm,object=m,margin_sign=direction))
  }
  if (nm=="LDA") {
    if (!requireNamespace("MASS",quietly=TRUE)) stop("Package MASS missing")
    m <- MASS::lda(x=X,grouping=factor(y,levels=0:1))
    return(list(method=nm,object=m))
  }
  if (nm=="NaiveBayes") {
    if (!requireNamespace("e1071",quietly=TRUE)) stop("Package e1071 missing")
    m <- e1071::naiveBayes(x=X,y=factor(y,levels=0:1))
    return(list(method=nm,object=m))
  }
  if (nm=="RF") {
    # Fixed RF parameters (not row-level OOB tuning).
    if (!requireNamespace("randomForestSRC",quietly=TRUE)) stop("Package randomForestSRC missing")
    set.seed(seed)
    dat <- data.frame(Type=factor(y,levels=0:1),X,check.names=FALSE)
    m <- randomForestSRC::rfsrc(Type~.,data=dat,ntree=500L,nodesize=5L,
                                importance=FALSE,proximity=FALSE,forest=TRUE)
    return(list(method=nm,object=m))
  }
  if (nm=="GBM") {
    m <- fit_gbm_grouped(X,y,donor,seed=seed,max_trees=700L)
    return(list(method=nm,object=m$fit,ntrees=m$ntrees))
  }
  if (nm=="XGBoost") {
    m <- fit_xgb_grouped(X,y,donor,seed=seed)
    return(list(method=nm,object=m))
  }
  if (nm=="plsRglm") {
    # No row-level tuning: fixed two latent components, learned on training donors only.
    if (!requireNamespace("plsRglm",quietly=TRUE)) stop("Package plsRglm missing")
    nt_use <- min(2L,ncol(X))
    set.seed(seed)
    m <- plsRglm::plsRglm(as.numeric(y),as.data.frame(X),nt=nt_use,
                          modele="pls-glm-logistic",verbose=FALSE,sparse=TRUE)
    return(list(method=nm,object=m,nt=nt_use))
  }
  stop("No audited fit for ",nm)
}

predict_class1_prob <- function(m,X) {
  X <- as.data.frame(X,check.names=FALSE)
  nm <- m$method
  if (nm %in% c("Lasso","Ridge","Enet")) {
    p <- stats::predict(m$object,newx=as.matrix(X),s="lambda.min",type="response")
  } else if (nm=="Stepglm") {
    p <- stats::predict(m$object,newdata=X,type="response")
  } else if (nm=="glmBoost") {
    p <- stats::predict(m$object,newdata=X,type="response")
  } else if (nm=="SVM") {
    # Logistic transform of oriented decision margin: ROC score, NOT calibrated probability.
    pred <- stats::predict(m$object,newdata=as.matrix(X),decision.values=TRUE)
    mar <- as.numeric(attr(pred,"decision.values"))
    if (length(mar)!=nrow(X) || any(!is.finite(mar))) stop("SVM margin invalid")
    p <- stats::plogis(m$margin_sign * mar)
  } else if (nm=="LDA") {
    p <- stats::predict(m$object,newdata=X)$posterior[,"1"]
  } else if (nm=="NaiveBayes") {
    p <- stats::predict(m$object,newdata=X,type="raw")[,"1"]
  } else if (nm=="RF") {
    p <- stats::predict(m$object,newdata=X)$predicted[,"1"]
  } else if (nm=="GBM") {
    p <- stats::predict(m$object,newdata=X,n.trees=m$ntrees,type="response")
  } else if (nm=="XGBoost") {
    p <- stats::predict(m$object,newdata=as.matrix(X))
  } else if (nm=="plsRglm") {
    raw <- stats::predict(m$object,newdata=X,type="response")
    # plsRglm may return predictions at each component count; use the last.
    if (is.matrix(raw)) raw <- raw[,ncol(raw)]
    if (is.array(raw) && length(dim(raw))>2L) {
      if (dim(raw)[1]!=nrow(X)) stop("Unexpected plsRglm prediction layout")
      raw <- matrix(raw,nrow=nrow(X))[,ncol(matrix(raw,nrow=nrow(X)))]
    }
    p <- raw
  } else stop("No score predictor for ",nm)
  p <- as.numeric(p)
  if (length(p)!=nrow(X) || any(!is.finite(p)) || any(p<0|p>1))
    stop("Invalid predictions from ", nm)
  p
}

fit_ml_pipeline <- function(method,X,y,donor,seed=123456L,min_features=2L) {
  parts <- parse_ml_method(method)
  if (!selector_supported(parts$selector$method))
    stop("SKIP: unavailable audited selector: ",parts$selector$method)
  if (!learner_supported(parts$learner$method))
    stop("SKIP: unavailable audited learner: ",parts$learner$method)
  scaler <- fit_group_scaler(X)
  xs <- transform_by_training(X,scaler)
  genes <- select_in_train(parts$selector,xs,y,donor,seed=seed+1L)
  genes <- intersect(colnames(xs),unique(genes))
  if (length(genes)<min_features) stop("SKIP: ",length(genes),
     " genes selected; minimum is ",min_features)
  trained <- fit_classifier(parts$learner,xs[,genes,drop=FALSE],y,donor,seed=seed+2L)
  list(label=parts$label,scaler=scaler,selected_genes=genes,
       learner=trained,train_donors=unique(as.character(donor)))
}

predict_ml_pipeline <- function(model,X) {
  z <- transform_by_training(X,model$scaler)
  predict_class1_prob(model$learner,z[,model$selected_genes,drop=FALSE])
}

auc_fixed_direction <- function(y, p) {
  if (!requireNamespace("pROC",quietly=TRUE)) stop("Install pROC")
  if (length(unique(y))!=2L || any(!is.finite(p))) stop("Invalid ROC inputs")
  as.numeric(pROC::auc(pROC::roc(response=as.numeric(y),predictor=as.numeric(p),
                                 levels=c(0,1),direction="<",quiet=TRUE)))
}

# Bootstrap complete paired donors, not individual biopsies, for 95% AUC CI.
paired_auc_bootstrap <- function(y,p,donor,B=1000L,seed=123456L) {
  d <- as.character(donor)
  tab <- table(d,factor(y,levels=0:1))
  if (any(tab!=1L)) stop("AUC bootstrap requires one CIT and one R2H per donor")
  set.seed(seed); u <- unique(d)
  val <- replicate(B, {
    pick <- sample(u,length(u),replace=TRUE)
    ix <- unlist(lapply(pick,function(x) which(d==x)),use.names=FALSE)
    auc_fixed_direction(y[ix],p[ix])
  })
  as.numeric(stats::quantile(val,c(.025,.975),na.rm=TRUE))
}



 }),envir=.GlobalEnv)
}
load_ml_helpers()

PIPELINE <- new.env(parent=.GlobalEnv)

# =============================================================================
# STAGE figure1_bulk_plots
# =============================================================================
PIPELINE$figure1_bulk_plots <- function() {
# scripts/01_bulk_heatmap_volcano.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

# =============================================================================
# Figure 1A-B | GSE145989 (53 paired CIT/R2H donors)
# Submission-style DEG heatmap + volcano, without ANY gene highlighting/labels.
# Uses existing paired-limma results; DOES NOT re-run differential analysis.
# Matches user script GSE145989_paired_limma_53pairs(2).R.
# =============================================================================

options(stringsAsFactors = FALSE)

# -------- 0. Settings -------------------------------------------------------
# Same directory used in the paired-limma script you uploaded.
ROOT_DIR <- REPO_ROOT
if (!dir.exists(ROOT_DIR)) ROOT_DIR <- getwd()

RESULT_DIR <- repo_output("Figure1/Paired_limma")
DEG_FILE <- repo_input("bulk_paired_results")
META_FILE <- repo_input("bulk_paired_metadata")
EXPR_FILE <- repo_input("bulk_train_expression")
OUT_DIR <- repo_output("Figure1/AB")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- log2(1.5)
TOP_N_UP <- 40L      # number of most significant upregulated genes for heatmap
TOP_N_DOWN <- 40L    # number of most significant downregulated genes
SHOW_GENE_NAMES_HEATMAP <- TRUE  # switch FALSE if 80 gene labels look crowded
DRAW_CUTOFF_LINES <- FALSE       # reference image has no visible cut-off lines
Z_LIMIT <- 2

# -------- 1. Packages -------------------------------------------------------
pkgs <- c("ggplot2", "pheatmap", "gridExtra")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install missing packages before analysis: ", paste(missing_pkgs, collapse=", "))
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(gridExtra)
})

files <- c(DEG_FILE, META_FILE, EXPR_FILE)
if (any(!file.exists(files))) {
  stop("Check source data or configuration. ", paste(files[!file.exists(files)], collapse = "\n"),
       "\n  ROOT_DIR   RESULT_DIR.")
}

# -------- 2. Read the PREVIOUSLY COMPUTED paired-limma results -------------
deg <- read.delim(DEG_FILE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
if (!all(c("Gene", "logFC", "adj.P.Val") %in% names(deg))) {
  stop("  Gene、logFC   adj.P.Val  ")
}
deg$Gene <- trimws(as.character(deg$Gene))
deg$logFC <- as.numeric(deg$logFC)
deg$adj.P.Val <- as.numeric(deg$adj.P.Val)
deg <- deg[nzchar(deg$Gene) & !is.na(deg$Gene) &
             is.finite(deg$logFC) & is.finite(deg$adj.P.Val) &
             deg$adj.P.Val >= 0 & deg$adj.P.Val <= 1, , drop = FALSE]
if (anyDuplicated(deg$Gene)) stop("Check source data or configuration. ")

deg$Sig <- ifelse(deg$adj.P.Val < FDR_CUTOFF & deg$logFC > LOGFC_CUTOFF,
                  "Up", ifelse(deg$adj.P.Val < FDR_CUTOFF &
                                  deg$logFC < -LOGFC_CUTOFF, "Down", "Not"))
deg$Sig <- factor(deg$Sig, levels = c("Down", "Not", "Up"))
deg$negLog10FDR <- -log10(pmax(deg$adj.P.Val, .Machine$double.xmin))
n_up <- sum(deg$Sig == "Up")
n_down <- sum(deg$Sig == "Down")
message("Paired limma DEGs (FDR < 0.05; |log2FC| > log2(1.5)): Up=", n_up,
        "; Down=", n_down, "; Not=", sum(deg$Sig == "Not"))

# -------- 3. Read sample metadata and check pairing ------------------------
meta <- read.delim(META_FILE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
if (!all(c("sample", "donor", "condition") %in% names(meta))) {
  stop("01_sample_metadata_106.tsv  sample、donor、condition  ")
}
meta$sample <- as.character(meta$sample)
meta$donor <- as.character(meta$donor)
meta$condition <- as.character(meta$condition)
if (anyDuplicated(meta$sample)) stop("meta ")
if (nrow(meta) != 106L || length(unique(meta$donor)) != 53L ||
    !identical(sort(unique(meta$condition)), c("CIT", "R2H"))) {
  stop(" 53 CIT/R2H； ")
}
tb <- table(meta$donor, meta$condition)
if (!all(tb[, "CIT"] == 1L & tb[, "R2H"] == 1L)) {
  stop("Check source data or configuration. ")
}

# Deterministic column arrangement: 53 CIT followed by 53 R2H in donor order.
# Matching donor order in both halves makes pairing auditable.
donors <- unique(meta$donor[meta$condition == "CIT"])
meta_key <- paste(meta$donor, meta$condition, sep = "::")
idx <- match(c(paste(donors, "CIT", sep = "::"),
               paste(donors, "R2H", sep = "::")), meta_key)
if (anyNA(idx)) stop("Check source data or configuration. ")
meta <- meta[idx, , drop = FALSE]
stopifnot(all(meta$condition[1:53] == "CIT"),
          all(meta$condition[54:106] == "R2H"),
          identical(meta$donor[1:53], meta$donor[54:106]))

# -------- 4. Select balanced up/down genes objectively ---------------------
raw <- read.delim(EXPR_FILE, header = TRUE, check.names = FALSE,
                  quote = "", comment.char = "", fileEncoding = "UTF-8-BOM")
if (!"id" %in% names(raw)) stop("Check source data or configuration. ")
genes <- as.character(raw$id)
if (anyNA(genes) || anyDuplicated(genes)) stop("Check source data or configuration. ")
if (!all(meta$sample %in% names(raw))) {
  stop("Check source data or configuration. ",
       paste(setdiff(meta$sample, names(raw)), collapse = ", "))
}

# Sort significant genes by BH-FDR, then by the magnitude of logFC.
up <- deg[deg$Sig == "Up" & deg$Gene %in% genes, , drop = FALSE]
down <- deg[deg$Sig == "Down" & deg$Gene %in% genes, , drop = FALSE]
up <- up[order(up$adj.P.Val, -up$logFC), , drop = FALSE]
down <- down[order(down$adj.P.Val, down$logFC), , drop = FALSE]
heat_genes <- unique(c(head(up$Gene, TOP_N_UP), head(down$Gene, TOP_N_DOWN)))
if (length(heat_genes) < 2L) {
  stop(" 2 FDR logFC ")
}

# The input is already log2-normalized: NEVER log2-transform again.
expression <- as.matrix(raw[match(heat_genes, genes), meta$sample, drop = FALSE])
rownames(expression) <- heat_genes
storage.mode(expression) <- "double"
if (any(!is.finite(expression))) stop("Check source data or configuration. ")

# Row Z-score, clipped ONLY FOR DISPLAY (not for any statistical testing).
z <- t(scale(t(expression)))
z[!is.finite(z)] <- 0
# IMPORTANT: do not put a scalar as the first argument to pmin/pmax;
# doing so can drop the matrix dim attribute, causing a pheatmap gaps error.
z[z > Z_LIMIT] <- Z_LIMIT
z[z < -Z_LIMIT] <- -Z_LIMIT
if (!is.matrix(z) || ncol(z) != nrow(meta)) {
  stop("Check source data or configuration. ", paste(dim(z), collapse = "x"))
}
if (ncol(z) != 106L || anyNA(colnames(z))) {
  stop("Check source data or configuration. ")
}
message("Heatmap matrix: ", nrow(z), " genes x ", ncol(z), " samples")
anno <- data.frame(Type = factor(meta$condition, levels = c("CIT", "R2H")))
rownames(anno) <- meta$sample
heat_palette <- colorRampPalette(c("#356EAB", "#F8F8F8", "#C74C4E"))(101)
type_colors <- list(Type = c(CIT = "#59BBC2", R2H = "#E6A2B0"))

# silent=TRUE returns grid gtable for PDF/PNG combining, does not draw instantly.
heat <- pheatmap::pheatmap(
  z,
  color = heat_palette,
  breaks = seq(-Z_LIMIT, Z_LIMIT, length.out = 102),
  annotation_col = anno,
  annotation_colors = type_colors,
  annotation_names_col = FALSE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,            # keep CIT left, R2H right; donor-matched order
  gaps_col = sum(meta$condition == "CIT"),  # separation after CIT columns
  show_colnames = FALSE,
  show_rownames = SHOW_GENE_NAMES_HEATMAP,
  fontsize_row = 5.5,
  fontsize = 9,
  treeheight_row = 32,
  border_color = NA,
  cellwidth = NA,
  cellheight = NA,
  legend_breaks = c(-2, 0, 2),
  legend_labels = c("-2", "0", "+2"),
  silent = TRUE
)

# -------- 5. Volcano (NO gene names, NO highlighted candidates) ------------
volcano_colors <- c(Down = "#3C70AA", Not = "#C2C2C2", Up = "#C64F50")
# Draw grey background first, color-coded significant points second.
deg <- deg[order(deg$Sig != "Not"), , drop = FALSE]

p_volcano <- ggplot2::ggplot(deg, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = Sig), size = 0.85, alpha = 0.83) +
  scale_color_manual(values = volcano_colors, limits = c("Down", "Not", "Up"),
                     drop = FALSE) +
  labs(x = "log2 fold change (R2H - CIT)",
       y = expression(-log[10]("BH-adjusted P")), color = "Sig") +
  theme_bw(base_size = 11.3) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    panel.border = element_rect(color = "#333333", linewidth = 0.55),
    panel.grid.major = element_line(color = "#EEEEEE", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", color = "#202020", size = 11),
    axis.text = element_text(color = "#444444", size = 9.2),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 9),
    legend.key = element_blank(),
    plot.margin = margin(7, 11, 7, 5)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

if (DRAW_CUTOFF_LINES) {
  p_volcano <- p_volcano +
    geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF),
               color = "#777777", linetype = "dashed", linewidth = 0.35) +
    geom_hline(yintercept = -log10(FDR_CUTOFF),
               color = "#777777", linetype = "dashed", linewidth = 0.35)
}

# -------- 6. Export A, B, and side-by-side A-B -------------------------------
# Use vector PDF, plus high-resolution PNG (with native graphics devices).
render_grob <- function(grob, pdf_file, png_file, width, height, dpi = 400) {
  grDevices::pdf(pdf_file, width = width, height = height, family = "sans",
                 useDingbats = FALSE)
  grid::grid.newpage()
  grid::grid.draw(grob)
  grDevices::dev.off()
  grDevices::png(png_file, width = width, height = height, units = "in",
                 res = dpi, bg = "white")
  grid::grid.newpage()
  grid::grid.draw(grob)
  grDevices::dev.off()
}

render_grob(heat$gtable,
            file.path(OUT_DIR, "Figure1A_DEG_heatmap.pdf"),
            file.path(OUT_DIR, "Figure1A_DEG_heatmap.png"),
            width = 10.2, height = 7.5)

ggsave(file.path(OUT_DIR, "Figure1B_volcano.pdf"),
       plot = p_volcano, width = 6.2, height = 5.9, bg = "white",
       device = "pdf")
ggsave(file.path(OUT_DIR, "Figure1B_volcano.png"),
       plot = p_volcano, width = 6.2, height = 5.9, dpi = 400, bg = "white")

# A/B panel labels outside the plot, so no overlap with heatmap legends.
panel_label <- function(letter) {
  grid::textGrob(letter, x = grid::unit(0.02, "npc"),
                 just = "left", gp = grid::gpar(fontface = "bold", fontsize = 24,
                                                col = "#151515"))
}
panel_A <- gridExtra::arrangeGrob(panel_label("A"), heat$gtable,
                                  ncol = 1, heights = c(0.075, 0.925))
panel_B <- gridExtra::arrangeGrob(panel_label("B"), ggplot2::ggplotGrob(p_volcano),
                                  ncol = 1, heights = c(0.075, 0.925))
combined <- gridExtra::arrangeGrob(panel_A, panel_B, ncol = 2,
                                   widths = c(1.28, 1))
render_grob(combined,
            file.path(OUT_DIR, "Figure1_AB_combined.pdf"),
            file.path(OUT_DIR, "Figure1_AB_combined.png"),
            width = 15.5, height = 7.0, dpi = 400)

# Export the exact objectively selected genes and counts for your supplement.
write.csv(data.frame(Gene = heat_genes,
                     Direction = as.character(deg$Sig[match(heat_genes, deg$Gene)])),
          file.path(OUT_DIR, "Figure1A_heatmap_selected_genes.csv"), row.names = FALSE)
write.csv(data.frame(Category = c("Up", "Down", "Not"),
                     Count = c(n_up, n_down, sum(deg$Sig == "Not"))),
          file.path(OUT_DIR, "Figure1B_DEG_counts.csv"), row.names = FALSE)

cat("\nDone! Saved to:\n", normalizePath(OUT_DIR, winslash = "/"), "\n", sep = "")
cat("Heatmap:", length(heat_genes), "genes (", min(TOP_N_UP, nrow(up)),
    " Up, ", min(TOP_N_DOWN, nrow(down)), " Down).\n", sep = "")
cat("Volcano is unlabelled; no gene selected or highlighted.\n")




}

# =============================================================================
# STAGE figure1_wgcna
# =============================================================================
PIPELINE$figure1_wgcna <- function() {
# scripts/02_wgcna.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

# ============================================================================
# GSE145989 - ONE-COHORT WGCNA, 2h CIT/R2H DONOR PAIRED
# Required input: GSE145989.normalize.txt ONLY
# Optional: no other metadata files; the GSM-to-donor map below is GEO-SOFT verified.
# This script replaces 'WGCNA_merge_normalize_corrected(2).R'.
# Results: WGCNA_GSE145989_PAIRED_2H/
# The module-trait significance uses donor-blocked paired limma, not independent
# Pearson correlation P values. Correlation is descriptive only.
# ============================================================================
options(stringsAsFactors = FALSE)

# Set working directory to the folder containing GSE145989.normalize.txt.
# E.g., uncomment and edit:
EXPR_FILE <- repo_input("bulk_train_expression")
OUT_DIR <- repo_output("Figure1/WGCNA")

# Network settings (8,000 genes, split into blocks to reduce Windows RAM use)
TOP_N_MAD_GENES <- 8000L
MIN_MAD <- 0.10
POWERS <- 1:20
NETWORK_TYPE <- "signed"
TOM_TYPE <- "signed"
MAX_BLOCK_SIZE <- 4000L
MIN_MODULE_SIZE <- 30L
DEEP_SPLIT <- 3L
MERGE_CUT_HEIGHT <- 0.20
MM_CUTOFF <- 0.80
GS_COR_CUTOFF <- 0.50   # descriptive threshold for scatterplot only
PAIRED_FDR_CUTOFF <- 0.05
PAIRED_LOGFC_CUTOFF <- log2(1.5)

for (pkg in c("WGCNA", "limma")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Check source data or configuration. ", pkg, ";  Bioconductor limma、CRAN WGCNA")
  }
}
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(limma))

if (!file.exists(EXPR_FILE)) stop("Check source data or configuration. ", EXPR_FILE, "Check source data or configuration. ", getwd())
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (dn in c("01_metadata", "02_sample_QC", "03_soft_power", "04_modules",
             "05_module_trait_paired", "06_gene_GS_MM", "07_module_genes",
             "08_hub_genes", "09_workspace")) {
  dir.create(file.path(OUT_DIR, dn), showWarnings = FALSE)
}
out <- function(...) file.path(OUT_DIR, ...)
write_tsv <- function(x, fname) {
  utils::write.table(x, file = out(fname), sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "")
}

# ==========================================================================
# 1. Original CEL/SOFT-verified donor identity map; 575 & 605 are 2H
# ==========================================================================
PAIR_MAP <- utils::read.table(text = "
DONOR CIT_GSM R2H_GSM
Toronto_571 GSM4343791 GSM4343761
Toronto_575 GSM4343769 GSM4343770
Toronto_578 GSM4343775 GSM4343776
Toronto_589 GSM4343779 GSM4343781
Toronto_602 GSM4343744 GSM4343746
Toronto_605 GSM4343784 GSM4343787
Toronto_606 GSM4343732 GSM4343733
Toronto_607 GSM4343795 GSM4343796
Toronto_614 GSM4343800 GSM4343801
Toronto_618 GSM4343676 GSM4343677
Toronto_626 GSM4343693 GSM4343694
Toronto_627 GSM4343777 GSM4343778
Toronto_635 GSM4343797 GSM4343798
Toronto_639 GSM4343807 GSM4343808
Toronto_645 GSM4343762 GSM4343763
Toronto_652 GSM4343666 GSM4343668
Toronto_654 GSM4343702 GSM4343704
Toronto_659 GSM4343771 GSM4343774
Toronto_667 GSM4343809 GSM4343810
Toronto_679 GSM4343678 GSM4343679
Toronto_684 GSM4343683 GSM4343684
Toronto_686 GSM4343685 GSM4343686
Toronto_690 GSM4343816 GSM4343817
Toronto_692 GSM4343742 GSM4343743
Toronto_699 GSM4343673 GSM4343675
Toronto_704 GSM4343741 GSM4343749
Toronto_705 GSM4343714 GSM4343715
Toronto_717 GSM4343669 GSM4343670
Toronto_719 GSM4343719 GSM4343720
Toronto_728 GSM4343811 GSM4343812
Toronto_740 GSM4343707 GSM4343708
Toronto_744 GSM4343736 GSM4343737
Toronto_747 GSM4343691 GSM4343692
Toronto_757 GSM4343680 GSM4343682
Toronto_758 GSM4343696 GSM4343697
Toronto_763 GSM4343700 GSM4343701
Toronto_772 GSM4343705 GSM4343706
Toronto_780 GSM4343671 GSM4343672
Toronto_785 GSM4343757 GSM4343758
Toronto_786 GSM4343750 GSM4343751
Toronto_793 GSM4343709 GSM4343711
Toronto_796 GSM4343734 GSM4343735
Toronto_824 GSM4343712 GSM4343713
Toronto_828 GSM4343729 GSM4343730
Toronto_832 GSM4343753 GSM4343754
Toronto_833 GSM4343782 GSM4343783
Toronto_834 GSM4343739 GSM4343740
Toronto_853 GSM4343802 GSM4343803
Toronto_855 GSM4343721 GSM4343690
Toronto_857 GSM4343766 GSM4343768
Toronto_858 GSM4343727 GSM4343728
Toronto_864 GSM4343793 GSM4343794
Toronto_898 GSM4343725 GSM4343726
", header = TRUE, stringsAsFactors = FALSE)
stopifnot(nrow(PAIR_MAP) == 53L,
          !anyDuplicated(PAIR_MAP$DONOR),
          !anyDuplicated(c(PAIR_MAP$CIT_GSM, PAIR_MAP$R2H_GSM)))
PAIR_MAP$CIT_column <- paste0(PAIR_MAP$CIT_GSM, "_Control")
PAIR_MAP$R2H_column <- paste0(PAIR_MAP$R2H_GSM, "_Treat")

# ==========================================================================
# 2. Load log2-normalized expression (NO additional log2 transform or ComBat)
# ==========================================================================
message("[1/9] Read expression: ", EXPR_FILE)
rt <- utils::read.delim(EXPR_FILE, header = TRUE, check.names = FALSE,
                        quote = "", comment.char = "", stringsAsFactors = FALSE,
                        fileEncoding = "UTF-8-BOM")
if (ncol(rt) < 3L || tolower(names(rt)[1]) != "id") {
  stop(" id GSMxxx_Control/GSMxxx_Treat")
}
gene_id <- as.character(rt[[1]])
if (anyNA(gene_id) || any(!nzchar(gene_id))) stop("Check source data or configuration. ")
expr <- as.matrix(rt[, -1, drop = FALSE])
suppressWarnings(storage.mode(expr) <- "double")
rownames(expr) <- gene_id
if (any(!is.finite(expr))) stop(" NA、Inf ")
if (anyDuplicated(colnames(expr))) stop("Check source data or configuration. ")
if (!all(grepl("^GSM[0-9]+_(Control|Treat)$", colnames(expr)))) {
  stop(" GSE145989 GSMxxxx_Control/Treat ")
}
if (anyDuplicated(rownames(expr))) {
  message(" gene ID； limma::avereps ")
  expr <- limma::avereps(expr, ID = rownames(expr))
}
message("  Raw matrix: ", nrow(expr), " genes x ", ncol(expr), " columns")
if (min(expr) < -5 || max(expr) > 30) {
  warning("Check source data or configuration. ")
}

# ==========================================================================
# 3. Keep COMPLETE 2h pairs by explicit GSM; exclude 1H, never use GSM parity
# ==========================================================================
has_cit <- PAIR_MAP$CIT_column %in% colnames(expr)
has_r2h <- PAIR_MAP$R2H_column %in% colnames(expr)
if (any(xor(has_cit, has_r2h))) {
  i <- which(xor(has_cit, has_r2h))
  stop("Check source data or configuration. ",
       paste(PAIR_MAP$DONOR[i], collapse = ", "))
}
paired <- PAIR_MAP[has_cit & has_r2h, , drop = FALSE]
if (nrow(paired) < 40L) {
  stop("Check source data or configuration. ", nrow(paired), " 2h GSE145989")
}
if (!nrow(paired) %in% c(52L, 53L)) {
  stop("Check source data or configuration. ", nrow(paired), "Check source data or configuration. ")
}
meta <- rbind(
  data.frame(Sample = paired$CIT_column, GSM = paired$CIT_GSM,
             Donor = paired$DONOR, Condition = "CIT"),
  data.frame(Sample = paired$R2H_column, GSM = paired$R2H_GSM,
             Donor = paired$DONOR, Condition = "R2H")
)
meta$Donor <- factor(meta$Donor, levels = paired$DONOR)
meta$Condition <- factor(meta$Condition, levels = c("CIT", "R2H"))
stopifnot(nrow(meta) == 2L * nrow(paired),
          all(table(meta$Donor, meta$Condition) == 1L))
unused <- setdiff(colnames(expr), meta$Sample)
expr <- expr[, meta$Sample, drop = FALSE]
stopifnot(identical(colnames(expr), meta$Sample))
message("  Original matrix samples: ", length(unused) + ncol(expr),
        "; included CIT/R2H: ", ncol(expr),
        " (", nrow(paired), " matched donors); excluded: ", length(unused))
write_tsv(meta, "01_metadata/GSE145989_CIT_R2H_sample_info.tsv")
write_tsv(paired, "01_metadata/paired_donors_kept.tsv")
write_tsv(
  data.frame(
    Sample = unused,
    Reason = rep(
      "not a complete verified CIT/R2H pair",
      length(unused)
    )
  ),
  "01_metadata/excluded_columns.tsv"
)

# ==========================================================================
# 4. Select variable genes: WGCNA samples x genes; not DEG prefilter
# ==========================================================================
message("[2/9] Select top variable genes by MAD ...")
mad_gene <- apply(expr, 1L, stats::mad)
good <- is.finite(mad_gene) & mad_gene > MIN_MAD
if (sum(good) < 200L) stop("Check source data or configuration. ")
mad_ranked <- sort(mad_gene[good], decreasing = TRUE)
n_top <- min(TOP_N_MAD_GENES, length(mad_ranked))
selected_genes <- names(mad_ranked)[seq_len(n_top)]
write_tsv(data.frame(Gene = names(mad_ranked), MAD = unname(mad_ranked)),
          "01_metadata/variable_gene_MAD_ranking.tsv")
datExpr <- t(expr[selected_genes, , drop = FALSE])
stopifnot(nrow(datExpr) == nrow(meta),
          identical(rownames(datExpr), meta$Sample))
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 2)
if (!isTRUE(gsg$allOK)) {
  if (any(!gsg$goodSamples)) {
    stop("WGCNA ",
         paste(rownames(datExpr)[!gsg$goodSamples], collapse = ","))
  }
  datExpr <- datExpr[, gsg$goodGenes, drop = FALSE]
  message("Check source data or configuration. ", sum(!gsg$goodGenes), "Check source data or configuration. ")
}
if (ncol(datExpr) < 200L) stop(" WGCNA")
trait_binary <- as.numeric(meta$Condition == "R2H")
write_tsv(data.frame(Sample = meta$Sample, Donor = meta$Donor,
                     Condition = meta$Condition, R2H = trait_binary),
          "01_metadata/WGCNA_traits.tsv")

# ==========================================================================
# 5. Sample QC and soft power
# ==========================================================================
message("[3/9] Sample clustering/QC ...")
sampleTree <- stats::hclust(stats::dist(datExpr), method = "average")
grDevices::pdf(out("02_sample_QC/sample_dendrogram.pdf"), width = 13, height = 7)
graphics::par(mar = c(0, 4, 3, 0))
graphics::plot(sampleTree, cex = 0.65, main = "GSE145989: CIT/R2H paired biopsies",
               xlab = "", sub = "")
grDevices::dev.off()

trait_colors <- WGCNA::numbers2colors(trait_binary, signed = FALSE)
grDevices::pdf(out("02_sample_QC/sample_dendrogram_with_trait.pdf"),
               width = 13, height = 7)
WGCNA::plotDendroAndColors(sampleTree, colors = trait_colors,
                           groupLabels = "R2H (1) vs CIT (0)",
                           main = "CIT/R2H sample clustering")
grDevices::dev.off()

message("[4/9] Soft threshold selection ...")
sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = POWERS,
                                networkType = NETWORK_TYPE, verbose = 5)
fitIndex <- sft$fitIndices
write_tsv(fitIndex, "03_soft_power/power_fit_indices.tsv")
fit_r2 <- -sign(fitIndex[, 3]) * fitIndex[, 2]
pass <- which(is.finite(fit_r2) & fit_r2 >= 0.80)
if (length(pass)) {
  power <- fitIndex[pass[1], 1]
} else {
  valid <- which(is.finite(fit_r2))
  if (!length(valid)) stop(" soft power ")
  power <- fitIndex[valid[which.max(fit_r2[valid])], 1]
  warning(" power R2 >= 0.8, power=", power,
          "Check source data or configuration. ")
}
writeLines(as.character(power), out("03_soft_power/chosen_power.txt"))
grDevices::pdf(out("03_soft_power/soft_threshold.pdf"), width = 11, height = 5)
graphics::par(mfrow = c(1, 2))
graphics::plot(fitIndex[, 1], fit_r2, pch = 19, col = "#2E6DA4",
               xlab = "Soft threshold power", ylab = "Signed scale-free fit R2")
graphics::text(fitIndex[, 1], fit_r2, labels = fitIndex[, 1], pos = 3, cex = 0.7)
graphics::abline(h = 0.8, lty = 2, col = "red")
graphics::plot(fitIndex[, 1], fitIndex[, 5], pch = 19, col = "#2E6DA4",
               xlab = "Soft threshold power", ylab = "Mean connectivity")
graphics::par(mfrow = c(1, 1))
grDevices::dev.off()
message("  Selected power: ", power)

# ==========================================================================
# 6. Detect co-expression modules (signed, blockwise for lower RAM use)
# ==========================================================================
message("[5/9] Detect WGCNA modules ...")
net <- WGCNA::blockwiseModules(
  datExpr, power = power, networkType = NETWORK_TYPE, TOMType = TOM_TYPE,
  maxBlockSize = min(MAX_BLOCK_SIZE, ncol(datExpr)),
  minModuleSize = MIN_MODULE_SIZE, deepSplit = DEEP_SPLIT,
  mergeCutHeight = MERGE_CUT_HEIGHT, reassignThreshold = 0,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = FALSE, verbose = 3
)
moduleColors <- WGCNA::labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)
MEs <- WGCNA::orderMEs(WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes)
write_tsv(data.frame(Gene = names(moduleColors), Module = unname(moduleColors)),
          "04_modules/gene_module_assignments.tsv")
write_tsv(as.data.frame(table(moduleColors)), "04_modules/module_sizes.tsv")
for (b in seq_along(net$dendrograms)) {
  grDevices::pdf(out(sprintf("04_modules/gene_dendrogram_block_%02d.pdf", b)),
                 width = 13, height = 6.5)
  WGCNA::plotDendroAndColors(
    net$dendrograms[[b]],
    moduleColors[net$blockGenes[[b]]],
    "Modules", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
    guideHang = 0.05, main = sprintf("Gene dendrogram - block %d", b)
  )
  grDevices::dev.off()
}

# ==========================================================================
# 7. MODULE association: limma ~ Donor + Condition (paired inference)
# ==========================================================================
message("[6/9] Donor-paired module-trait analysis ...")
if (!identical(rownames(MEs), rownames(datExpr))) {
  MEs <- MEs[rownames(datExpr), , drop = FALSE]
}
ME_names <- setdiff(colnames(MEs), "MEgrey")
if (!length(ME_names)) stop(" grey ")
MEs_for_test <- MEs[, ME_names, drop = FALSE]
design <- stats::model.matrix(~ Donor + Condition, data = meta)
if (qr(design)$rank != ncol(design) ||
    !"ConditionR2H" %in% colnames(design)) stop("Check source data or configuration. ")
fit_ME <- limma::eBayes(limma::lmFit(t(MEs_for_test), design))
tt_ME <- limma::topTable(fit_ME, coef = "ConditionR2H", number = Inf,
                          adjust.method = "BH", sort.by = "none")
ME_cor <- vapply(ME_names, function(nm) {
  stats::cor(as.numeric(MEs_for_test[, nm]), trait_binary)
}, numeric(1))
moduleStats <- data.frame(
  Module = sub("^ME", "", ME_names),
  ME_name = ME_names,
  Mean_ME_R2H_minus_CIT = unname(tt_ME[ME_names, "logFC"]),
  Descriptive_correlation_R2H = as.numeric(ME_cor[ME_names]),
  Paired_t = unname(tt_ME[ME_names, "t"]),
  Paired_P = unname(tt_ME[ME_names, "P.Value"]),
  Paired_FDR = unname(tt_ME[ME_names, "adj.P.Val"]),
  N_donors = nlevels(meta$Donor), stringsAsFactors = FALSE
)
moduleStats <- moduleStats[order(moduleStats$Paired_FDR), , drop = FALSE]
write_tsv(moduleStats, "05_module_trait_paired/module_trait_paired_limma.tsv")

# Descriptive module-trait correlation colors, paired limma FDR labels.
mat <- matrix(moduleStats$Descriptive_correlation_R2H, ncol = 1,
              dimnames = list(paste0("ME", moduleStats$Module), "R2H vs CIT"))
txt <- matrix(paste0(sprintf("%.2f", moduleStats$Descriptive_correlation_R2H),
                     "\nFDR=", formatC(moduleStats$Paired_FDR, digits = 2, format = "g")),
              ncol = 1)
grDevices::pdf(out("05_module_trait_paired/module_trait_paired_FDR_heatmap.pdf"),
               width = 6.5, height = max(5.5, 0.40 * nrow(mat) + 1.5))
graphics::par(mar = c(5, 10, 3, 3))
WGCNA::labeledHeatmap(Matrix = mat, xLabels = colnames(mat),
                       yLabels = rownames(mat), ySymbols = rownames(mat),
                       textMatrix = txt, colors = WGCNA::blueWhiteRed(50),
                       colorLabels = FALSE, setStdMargins = FALSE,
                       cex.text = 0.8, zlim = c(-1, 1),
                       main = "Module-trait association\npaired donor-adjusted FDR")
grDevices::dev.off()

# ==========================================================================
# 8. GENE significance: paired limma P values, not unpaired correlation P
# ==========================================================================
message("[7/9] Gene significance, module membership ...")
fit_G <- limma::eBayes(limma::lmFit(t(datExpr), design))
tt_G <- limma::topTable(fit_G, coef = "ConditionR2H", number = Inf,
                         adjust.method = "BH", sort.by = "none")
stopifnot(all(rownames(tt_G) %in% colnames(datExpr)))
gene_names <- colnames(datExpr)
mm <- stats::cor(datExpr, MEs, use = "pairwise.complete.obs")
colnames(mm) <- paste0("MM_", sub("^ME", "", colnames(mm)))
gs_r <- as.numeric(stats::cor(datExpr, trait_binary, use = "pairwise.complete.obs"))
mm_own <- vapply(seq_along(gene_names), function(i) {
  nm <- paste0("MM_", moduleColors[gene_names[i]])
  if (nm %in% colnames(mm)) mm[i, nm] else NA_real_
}, numeric(1))
geneStats <- data.frame(
  Gene = gene_names,
  Module = unname(moduleColors[gene_names]),
  MM_own_module = mm_own,
  GS_R2H_descriptive_correlation = gs_r,
  Paired_log2FC_R2H_minus_CIT = unname(tt_G[gene_names, "logFC"]),
  Paired_t = unname(tt_G[gene_names, "t"]),
  Paired_P = unname(tt_G[gene_names, "P.Value"]),
  Paired_FDR = unname(tt_G[gene_names, "adj.P.Val"]),
  stringsAsFactors = FALSE
)
# MM matrix (all module memberships), ordered in the same row order.
geneStats <- cbind(geneStats, as.data.frame(mm, check.names = FALSE))
write_tsv(geneStats, "06_gene_GS_MM/all_gene_GS_MM_paired.tsv")

for (mod in setdiff(unique(moduleColors), "grey")) {
  z <- geneStats[geneStats$Module == mod, , drop = FALSE]
  if (nrow(z) < 3L) next
  grDevices::pdf(out(paste0("06_gene_GS_MM/MM_vs_GS_", mod, ".pdf")),
                 width = 6, height = 5.5)
  graphics::plot(abs(z$MM_own_module), abs(z$GS_R2H_descriptive_correlation),
                 pch = 16, cex = 0.6, col = mod, xlim = c(0, 1), ylim = c(0, 1),
                 xlab = paste("|Module membership|:", mod),
                 ylab = "|Descriptive gene-condition correlation|",
                 main = paste("MM vs GS:", mod))
  graphics::abline(v = MM_CUTOFF, h = GS_COR_CUTOFF,
                   col = "red", lty = 2)
  grDevices::dev.off()
}

# ==========================================================================
# 9. Module and hub gene exports: no predetermined target genes
# ==========================================================================
message("[8/9] Export module lists / donor-aware hub candidate sets ...")
hubCounts <- list()
for (mod in unique(moduleColors)) {
  z <- geneStats[geneStats$Module == mod, , drop = FALSE]
  z <- z[order(-abs(z$MM_own_module)), , drop = FALSE]
  write_tsv(z, paste0("07_module_genes/module_", mod, "_all_genes.tsv"))
  if (mod == "grey") next
  # Candidate hub definition: strong own-module membership AND paired-DE signal.
  # This is exploratory, not a validated causal driver or independent replication.
  hub <- z[is.finite(z$MM_own_module) & abs(z$MM_own_module) >= MM_CUTOFF &
           z$Paired_FDR < PAIRED_FDR_CUTOFF &
           abs(z$Paired_log2FC_R2H_minus_CIT) >= PAIRED_LOGFC_CUTOFF,
           , drop = FALSE]
  hub <- hub[order(hub$Paired_FDR, -abs(hub$MM_own_module)), , drop = FALSE]
  write_tsv(hub, paste0("08_hub_genes/module_", mod, "_paired_hub_candidates.tsv"))
  hubCounts[[mod]] <- data.frame(Module = mod, Candidate_hub_count = nrow(hub),
                                 Top_genes = paste(head(hub$Gene, 10), collapse = ","))
}
if (length(hubCounts)) {
  write_tsv(do.call(rbind, hubCounts), "08_hub_genes/hub_candidate_summary.tsv")
}

saveRDS(list(
  meta = meta, paired = paired, datExpr = datExpr, softPower = power,
  moduleColors = moduleColors, MEs = MEs, moduleStats = moduleStats,
  geneStats = geneStats, settings = list(N_TOP = TOP_N_MAD_GENES,
  networkType = NETWORK_TYPE, TOMType = TOM_TYPE)
), out("09_workspace/GSE145989_WGCNA_paired_results.rds"))
writeLines(capture.output(sessionInfo()), out("09_workspace/sessionInfo.txt"))
message("[9/9] Completed.  Pairs: ", nlevels(meta$Donor),
        "; samples: ", nrow(meta), "; WGCNA genes: ", ncol(datExpr),
        "; modules (excluding grey): ", length(ME_names))
message("Results: ", normalizePath(OUT_DIR, winslash = "/"))
message("MOST IMPORTANT: 05_module_trait_paired/module_trait_paired_limma.tsv")




# ============================================================================
# GSE145989: RE-DRAW FIGURE C + D ONLY. Run in the CURRENT R session.
# No data reads, no network inference, no re-normalization, no limma fitting.
# Required objects from existing WGCNA run: net, moduleColors, moduleStats.
# ============================================================================

stopifnot(exists("net"), exists("moduleColors"), exists("moduleStats"))
if (!requireNamespace("WGCNA", quietly = TRUE)) stop(" WGCNA ")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("  install.packages('ggplot2')")

fig_dir <- if (exists("OUT_DIR") && dir.exists(OUT_DIR)) {
  file.path(OUT_DIR, "10_pretty_C_D")
} else {
  "WGCNA_pretty_C_D"
}
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------- C: gene dendrogram + module colors -------------------
# The original network used blockwiseModules(maxBlockSize = 4000):
# each block has a DIFFERENT dendrogram. Never show block 1 as all genes.
if (length(net$dendrograms) != length(net$blockGenes)) {
  stop("net$dendrograms net$blockGenes ")
}

plot_C <- function(block_idx) {
  gene_idx <- net$blockGenes[[block_idx]]
  cols <- unname(moduleColors[gene_idx])
  if (length(cols) != length(gene_idx) || anyNA(cols)) {
    stop("C block block ", block_idx)
  }
  WGCNA::plotDendroAndColors(
    dendro = net$dendrograms[[block_idx]],
    colors = cols,
    groupLabels = "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    cex.colorLabels = 0.85,
    main = sprintf("Gene dendrogram and module colors  |  Block %d/%d",
                   block_idx, length(net$dendrograms))
  )
}

for (b in seq_along(net$dendrograms)) {
  fname <- sprintf("Figure_C_gene_dendrogram_block_%02d", b)
  grDevices::pdf(file.path(fig_dir, paste0(fname, ".pdf")),
                 width = 12, height = 5.8, useDingbats = FALSE)
  plot_C(b)
  grDevices::dev.off()
  
  grDevices::png(file.path(fig_dir, paste0(fname, ".png")),
                 width = 3600, height = 1740, res = 300, bg = "white")
  plot_C(b)
  grDevices::dev.off()
}

# ---------------------- D: module-trait heatmap ------------------------------
# Correlation is descriptive; q-values are from the DONOR-PAIRED limma.
required <- c("Module", "Descriptive_correlation_R2H", "Paired_FDR")
if (!all(required %in% colnames(moduleStats))) {
  stop("moduleStats : ", paste(setdiff(required, names(moduleStats)), collapse = ", "))
}
ms <- as.data.frame(moduleStats[, required, drop = FALSE])
ms$Module <- as.character(ms$Module)
ms$Descriptive_correlation_R2H <- as.numeric(ms$Descriptive_correlation_R2H)
ms$Paired_FDR <- as.numeric(ms$Paired_FDR)
ms <- ms[is.finite(ms$Descriptive_correlation_R2H) &
           is.finite(ms$Paired_FDR), , drop = FALSE]
if (!nrow(ms)) stop("Check source data or configuration. ")
# Keep the original paired-FDR ranking from moduleStats.
n <- nrow(ms)
ms$y <- n:1
ms$module_label <- paste0("ME", ms$Module)

# CIT and R2H are complementary binary labels: correlations have opposite sign.
# Their q-values are the SAME paired comparison, not separate tests.
plot_data <- rbind(
  data.frame(x = 1, y = ms$y, Trait = "CIT",
             Correlation = -ms$Descriptive_correlation_R2H,
             q = ms$Paired_FDR),
  data.frame(x = 2, y = ms$y, Trait = "R2H",
             Correlation = ms$Descriptive_correlation_R2H,
             q = ms$Paired_FDR)
)
format_q <- function(q) {
  ifelse(q < 0.001,
         formatC(q, format = "e", digits = 1),
         sprintf("%.3f", q))
}
plot_data$label <- sprintf("%.2f\n(q=%s)",
                           plot_data$Correlation, format_q(plot_data$q))

# R/WGCNA color names are used without remapping module membership.
color_strip <- data.frame(y = ms$y, color = ms$Module)
tryCatch(grDevices::col2rgb(color_strip$color),
         error = function(e) stop("Check source data or configuration. ", conditionMessage(e)))

pD <- ggplot2::ggplot() +
  ggplot2::geom_tile(
    data = plot_data,
    ggplot2::aes(x = x, y = y, fill = Correlation),
    width = 0.975, height = 0.965, color = "white", linewidth = 0.45
  ) +
  ggplot2::geom_tile(
    data = color_strip,
    ggplot2::aes(x = 0.38, y = y),
    fill = color_strip$color,
    width = 0.15, height = 0.965,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_text(
    data = plot_data,
    ggplot2::aes(x = x, y = y, label = label),
    size = 2.85, color = "#222222", lineheight = 0.92
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#406F9D", mid = "#F7F7F7", high = "#BD5960",
    midpoint = 0, limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Correlation\n(descriptive)"
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(1, 2), labels = c("CIT", "R2H"),
    limits = c(0.24, 2.52), expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq_len(n), labels = rev(ms$module_label),
    expand = ggplot2::expansion(add = 0.5)
  ) +
  ggplot2::labs(
    title = "Module–trait relationships",
    subtitle = "53 donor-matched CIT / 2h reperfusion biopsies",
    x = NULL, y = NULL,
    caption = "Cell text: descriptive r; q = donor-paired limma BH-FDR.\nCIT/R2H columns encode the same paired contrast."
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(face = "bold", size = 9.3,
                                        color = "#202020",
                                        margin = ggplot2::margin(r = 7)),
    axis.text.x = ggplot2::element_text(face = "bold", size = 11,
                                        color = "#202020"),
    plot.title = ggplot2::element_text(face = "bold", size = 14,
                                       hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = 9, color = "#555555",
                                          hjust = 0.5),
    plot.caption = ggplot2::element_text(size = 8, color = "#666666",
                                         hjust = 0, margin = ggplot2::margin(t = 8)),
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 9),
    plot.margin = ggplot2::margin(10, 15, 10, 10)
  )

figure_D_height <- max(6.5, 0.39 * n + 2.4)
ggplot2::ggsave(file.path(fig_dir, "Figure_D_module_trait_paired_pretty.pdf"),
                plot = pD, width = 6.6, height = figure_D_height,
                device = "pdf", bg = "white")
ggplot2::ggsave(file.path(fig_dir, "Figure_D_module_trait_paired_pretty.png"),
                plot = pD, width = 6.6, height = figure_D_height,
                dpi = 400, bg = "white")

print(pD)
message(" WGCNA limma.\n : ",
        normalizePath(fig_dir, winslash = "/"))

saveRDS(list(datExpr=datExpr,moduleColors=moduleColors,meta=meta,MEs=MEs,moduleStats=moduleStats,geneStats=geneStats), file.path(repo_output("Figure1/WGCNA"),"wgcna_checkpoint.rds"))

}

# =============================================================================
# STAGE figure1_blue_enrichment
# =============================================================================
PIPELINE$figure1_blue_enrichment <- function() {
# scripts/02_blue_module_enrichment.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.
restore_named_list(repo_input("wgcna_checkpoint"))
# =============================================================================
# GSE145989: BLUE module | GO (BP/CC/MF) + KEGG + 3 polished WGCNA figures
# =============================================================================
# Run AFTER WGCNA_GSE145989_only_strict_2h_PAIRED_FIXED.R in the SAME R session.
# This script DOES NOT: setwd(), load old RData, read the expression matrix,
# rebuild the WGCNA network, or re-run limma.
# It uses current objects: datExpr, moduleColors, meta, MEs, moduleStats, geneStats.
# The only optional input file is module_blue.txt (or the new blue-gene TSV).
# Enrichment background = ALL genes actually used in the WGCNA (datExpr columns).
# Statistical labels in the module-trait heatmap are EXISTING paired limma FDR.
# KEGG uses current online KEGG annotations and may require internet access.
# =============================================================================

options(stringsAsFactors = FALSE)

# ----------------------------- USER OPTIONS -----------------------------------
BLUE_FILE <- NA_character_  # e.g. "07_module_genes/module_blue.txt"; NA = auto-locate
OUT_DIR <- repo_output("Figure1/Blue_module_enrichment")
TOP_TERMS <- 15L
FDR_CUTOFF <- 0.05
MM_CUTOFF <- 0.80
GS_COR_CUTOFF <- 0.50  # descriptive correlation; NOT a paired P value
# -----------------------------------------------------------------------------

required_packages <- c("ggplot2", "pheatmap", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Check source data or configuration. ", paste(missing_packages, collapse = ", "),
       "Check source data or configuration. ",
       "install.packages(c('ggplot2','pheatmap'))\n",
       "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
       "BiocManager::install(c('clusterProfiler','org.Hs.eg.db','AnnotationDbi'), ask=FALSE, update=FALSE)")
}

need <- c("datExpr", "moduleColors", "meta", "MEs", "moduleStats", "geneStats")
absent <- need[!vapply(need, exists, logical(1), envir = .GlobalEnv,
                       inherits = FALSE)]
if (length(absent)) {
  stop("Check source data or configuration. ", paste(absent, collapse = ", "),
       ". RStudio  GSE145989  WGCNA \n",
       "  WGCNA_merge_normalize_corrected_workspace.RData.")
}

# ------------------------- 1. Verify object provenance ------------------------
stopifnot(is.matrix(datExpr) || is.data.frame(datExpr))
stopifnot(is.data.frame(meta), is.data.frame(moduleStats), is.data.frame(geneStats))
stopifnot(nrow(datExpr) == nrow(meta), length(moduleColors) == ncol(datExpr))
stopifnot(all(c("Sample", "Donor", "Condition") %in% colnames(meta)))
stopifnot(all(c("Module", "Descriptive_correlation_R2H", "Paired_FDR") %in%
              colnames(moduleStats)))
stopifnot(all(c("Gene", "Module", "MM_own_module",
                "GS_R2H_descriptive_correlation") %in% colnames(geneStats)))
stopifnot(!anyDuplicated(colnames(datExpr)), !anyDuplicated(geneStats$Gene))
if (!identical(rownames(datExpr), as.character(meta$Sample))) {
  stop("datExpr meta$Sample / ")
}
if (!all(grepl("^GSM434[0-9]+_(Control|Treat)$", meta$Sample))) {
  stop(" GSE145989 GSM434..._Control/Treat ",
       " WGCNA ")
}
if (!setequal(as.character(meta$Condition), c("CIT", "R2H"))) {
  stop(" meta CIT/R2H, ")
}
pair_tab <- table(meta$Donor, meta$Condition)
if (nrow(pair_tab) < 50L || any(pair_tab != 1L)) {
  stop(" meta GSE145989 CIT/R2H ")
}
if (!"blue" %in% as.character(moduleStats$Module) ||
    !"blue" %in% as.character(geneStats$Module)) {
  stop(" WGCNA blue ")
}
if (anyDuplicated(moduleStats$Module)) stop("moduleStats ")
if (anyNA(moduleStats$Paired_FDR)) {
  stop("moduleStats$Paired_FDR NA； limma FDR.")
}

wgcna_genes <- colnames(datExpr)
if (is.null(names(moduleColors))) {
  if (length(moduleColors) != length(wgcna_genes)) stop("Check source data or configuration. ")
  names(moduleColors) <- wgcna_genes
}
moduleColors <- moduleColors[wgcna_genes]
if (anyNA(moduleColors)) stop("moduleColors datExpr ")
true_blue <- wgcna_genes[as.character(moduleColors) == "blue"]
if (length(true_blue) < 10L) stop(" blue ")
if (!setequal(true_blue, geneStats$Gene[geneStats$Module == "blue"])) {
  stop("geneStats moduleColors blue ")
}

cat("[QC] Paired donors:", nrow(pair_tab), "/ samples:", nrow(meta), "\n")
cat("[QC] WGCNA genes:", length(wgcna_genes),
    "/ current blue module:", length(true_blue), "\n")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(OUT_DIR, "01_beautified_WGCNA_plots")
enrich_dir <- file.path(OUT_DIR, "02_blue_GO_KEGG")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(obj, path) {
  utils::write.table(obj, path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "")
}

BLUE <- "#3F6FA6"
LIGHT_BLUE <- "#A3C2D8"
RED <- "#B64E54"
LIGHT_RED <- "#E3ADB0"
WHITE <- "#FAFAFA"
GREY <- "#BFC5CC"

pretty_theme <- ggplot2::theme_classic(base_size = 12, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#606B78", size = 10),
    axis.text = ggplot2::element_text(color = "#333A44"),
    axis.title = ggplot2::element_text(face = "bold", color = "#23272C"),
    panel.grid.major.x = ggplot2::element_line(color = "#EDF1F5", linewidth = 0.3),
    legend.title = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(12, 15, 12, 12)
  )

save_gg <- function(p, stem, width = 8, height = 6) {
  ggplot2::ggsave(file.path(plot_dir, paste0(stem, ".pdf")), p,
                  width = width, height = height, bg = "white", device = "pdf",
                  useDingbats = FALSE)
  ggplot2::ggsave(file.path(plot_dir, paste0(stem, ".png")), p,
                  width = width, height = height, dpi = 400, bg = "white")
}

# ---------------- 2. Polished plot 1: paired module–trait heatmap -------------
# Do not use corPvalueStudent() on 106 dependent biopsies.
# CIT and R2H correlation signs are opposite, but both are ONE paired contrast;
# its q-value is the precomputed paired-limma BH FDR from moduleStats.
ms <- moduleStats
ms$Module <- as.character(ms$Module)
ms <- ms[order(ms$Module != "blue", ms$Paired_FDR, -abs(ms$Descriptive_correlation_R2H)),
         , drop = FALSE]
mod_name <- paste0("ME", ms$Module)
cor_R2H <- as.numeric(ms$Descriptive_correlation_R2H)
cor_mat <- cbind(CIT = -cor_R2H, R2H = cor_R2H)
rownames(cor_mat) <- mod_name

format_q <- function(q) {
  ifelse(q < 0.001, formatC(q, format = "e", digits = 1),
         formatC(q, format = "f", digits = 3))
}
q_lab <- format_q(ms$Paired_FDR)
lab <- matrix("", nrow = nrow(cor_mat), ncol = 2L,
              dimnames = dimnames(cor_mat))
for (i in seq_len(nrow(cor_mat))) {
  lab[i, ] <- paste0(sprintf("%.2f", cor_mat[i, ]), "\nq=", q_lab[i])
}

row_anno <- data.frame(Module = factor(ms$Module, levels = unique(ms$Module)))
rownames(row_anno) <- mod_name
known_colors <- unique(ms$Module)
# WGCNA's named modules are R color strings.
for (cc in known_colors) {
  tryCatch(grDevices::col2rgb(cc), error = function(e) {
    stop("Check source data or configuration. ", cc)
  })
}
anno_colors <- list(Module = stats::setNames(known_colors, known_colors))
heat_cols <- grDevices::colorRampPalette(c("#3E6D9C", "#C9DCEB", "#FAFAFA",
                                           "#E7BFC1", "#B95159"))(100)
plot_heat <- pheatmap::pheatmap(
  mat = cor_mat, color = heat_cols,
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_row = row_anno, annotation_colors = anno_colors,
  annotation_names_row = FALSE, annotation_legend = FALSE,
  display_numbers = lab, number_color = "#20262B",
  fontsize = 10, fontsize_row = 10, fontsize_col = 11,
  fontsize_number = 7.4,
  cellwidth = 83, cellheight = 23,
  border_color = "#F0F0F0", angle_col = 0,
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  legend_labels = c("-1", "-0.5", "0", "0.5", "1"),
  main = "Module–trait relationships (paired donors)",
  silent = TRUE
)

heat_height <- max(6.0, 0.33 * nrow(cor_mat) + 2.6)
for (ext in c("pdf", "png")) {
  fname <- file.path(plot_dir, paste0("Figure_1_module_trait_paired_BLUE_focus.", ext))
  if (ext == "pdf") grDevices::pdf(fname, width = 7.2, height = heat_height,
                                    useDingbats = FALSE)
  else grDevices::png(fname, width = 7.2, height = heat_height,
                      units = "in", res = 400, bg = "white")
  grid::grid.newpage()
  grid::grid.draw(plot_heat$gtable)
  grDevices::dev.off()
}
write_tsv(data.frame(Module = ms$Module, Cor_CIT = -cor_R2H,
                     Cor_R2H = cor_R2H, Paired_FDR = ms$Paired_FDR),
          file.path(plot_dir, "Figure_1_module_trait_source.tsv"))

# ----------- 3. Polished plot 2: module significance bar plot -----------------
# Descriptive mean |gene–R2H correlation|, NOT inferential paired p-value.
gs <- geneStats[, c("Gene", "Module", "GS_R2H_descriptive_correlation")]
names(gs)[3] <- "GS"
gs <- gs[is.finite(gs$GS) & !is.na(gs$Module), , drop = FALSE]
module_sum <- stats::aggregate(abs(gs$GS), by = list(Module = gs$Module), FUN = mean)
names(module_sum)[2] <- "Mean_abs_GS"
module_count <- table(gs$Module)
module_sum$N_genes <- as.integer(module_count[module_sum$Module])
module_sum <- module_sum[order(module_sum$Mean_abs_GS), , drop = FALSE]
module_sum$Module <- factor(module_sum$Module, levels = module_sum$Module)
module_sum$Focus <- ifelse(as.character(module_sum$Module) == "blue", "Blue", "Other")

p_ms <- ggplot2::ggplot(module_sum,
                        ggplot2::aes(x = Module, y = Mean_abs_GS, fill = Focus)) +
  ggplot2::geom_col(width = 0.68, color = "white", linewidth = 0.28) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = c(Blue = BLUE, Other = "#BCCDDC"), guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
  ggplot2::labs(title = "Module significance", subtitle = "Blue module highlighted",
                x = NULL, y = "Mean |gene–R2H correlation| (descriptive)") +
  pretty_theme +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

save_gg(p_ms, "Figure_2_module_significance_BLUE_focus", 8.4,
        max(5.3, 0.29 * nrow(module_sum) + 1.3))
write_tsv(module_sum, file.path(plot_dir, "Figure_2_module_significance_source.tsv"))

# -------------- 4. Polished plot 3: BLUE module MM–GS -------------------------
blue_gs <- geneStats[as.character(geneStats$Module) == "blue", , drop = FALSE]
blue_gs <- blue_gs[is.finite(blue_gs$MM_own_module) &
                     is.finite(blue_gs$GS_R2H_descriptive_correlation), , drop = FALSE]
if (nrow(blue_gs) < 3L) stop("blue MM–GS ")
blue_scatter <- data.frame(Gene = blue_gs$Gene,
                           MM = abs(blue_gs$MM_own_module),
                           GS = abs(blue_gs$GS_R2H_descriptive_correlation))

p_mmgs <- ggplot2::ggplot(blue_scatter, ggplot2::aes(x = MM, y = GS)) +
  ggplot2::annotate("rect", xmin = MM_CUTOFF, xmax = 1,
                     ymin = GS_COR_CUTOFF, ymax = 1,
                     fill = "#D8E8F4", alpha = 0.65) +
  ggplot2::geom_hline(yintercept = GS_COR_CUTOFF, linetype = "dashed",
                      color = "#A96060", linewidth = 0.55) +
  ggplot2::geom_vline(xintercept = MM_CUTOFF, linetype = "dashed",
                      color = "#A96060", linewidth = 0.55) +
  ggplot2::geom_point(color = BLUE, alpha = 0.70, size = 1.7) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  ggplot2::labs(title = "Blue module: membership vs. gene significance",
                subtitle = paste0(nrow(blue_scatter), " genes; donor-paired inference not implied"),
                x = "|Module membership (MM)|",
                y = "|Gene–R2H correlation (GS)|") +
  pretty_theme

save_gg(p_mmgs, "Figure_3_BLUE_module_MM_vs_GS", 7.3, 6.3)
write_tsv(blue_scatter, file.path(plot_dir, "Figure_3_blue_MM_GS_source.tsv"))
cat("[Plots] 3 beautified WGCNA figures saved to:", plot_dir, "\n")

# ---------------------- 5. BLUE module gene file -----------------------------
# Prefer the user's actual module_blue.txt. If absent, the current network
# moduleColors can regenerate exactly the same full blue-gene set.
search_paths <- c(
  "module_blue.txt",
  file.path("07_module_genes", "module_blue.txt"),
  file.path("WGCNA_GSE145989_PAIRED_2H", "07_module_genes", "module_blue.txt"),
  file.path("WGCNA_GSE145989_PAIRED_2H", "07_module_genes", "module_blue_all_genes.tsv")
)
if (is.character(BLUE_FILE) && length(BLUE_FILE) == 1L &&
    !is.na(BLUE_FILE) && nzchar(BLUE_FILE)) {
  if (!file.exists(BLUE_FILE)) stop(" BLUE_FILE ", BLUE_FILE)
  blue_source <- BLUE_FILE
} else {
  found <- search_paths[file.exists(search_paths)]
  blue_source <- if (length(found)) found[1] else NA_character_
}

read_blue_list <- function(path) {
  if (grepl("\\.tsv$", path, ignore.case = TRUE)) {
    tmp <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!"Gene" %in% colnames(tmp)) stop(" TSV Gene ", path)
    z <- as.character(tmp$Gene)
  } else {
    z <- readLines(path, warn = FALSE, encoding = "UTF-8")
  }
  z <- trimws(gsub('^"|"$', "", z))
  unique(z[!is.na(z) & nzchar(z) & z != "Gene"])
}

if (!is.na(blue_source)) {
  blue_genes <- read_blue_list(blue_source)
  missing_from_file <- setdiff(true_blue, blue_genes)
  extra_in_file <- setdiff(blue_genes, true_blue)
  if (length(missing_from_file) || length(extra_in_file)) {
    stop("module_blue GSE145989 blue \n ", blue_source,
         "\n blue  ( 10 ): ",
         paste(head(missing_from_file, 10), collapse = ", "),
         "Check source data or configuration. ",
         paste(head(extra_in_file, 10), collapse = ", "),
         "\n module_blue.txt, ")
  }
  message("[BLUE]  ", blue_source)
} else {
  blue_genes <- true_blue
  message("[BLUE]  module_blue.txt； moduleColors ")
}
blue_genes <- intersect(wgcna_genes, blue_genes)
writeLines(blue_genes, file.path(enrich_dir, "module_blue_verified_genes.txt"))

# --------------------- 6. SYMBOL -> ENTREZ mapping ---------------------------
# Background MUST be the gene universe entering the WGCNA network,
# not all annotated genes and not just the DEGs.
message("[GO/KEGG] Mapping SYMBOL -> ENTREZID for human genes ...")
map_table <- function(symbols) {
  symbols <- unique(as.character(symbols))
  z <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                               keys = symbols,
                                               keytype = "SYMBOL",
                                               columns = c("SYMBOL", "ENTREZID")))
  z <- z[!is.na(z$SYMBOL) & !is.na(z$ENTREZID) & nzchar(z$ENTREZID),
         c("SYMBOL", "ENTREZID"), drop = FALSE]
  unique(z)
}
background_map <- map_table(wgcna_genes)
blue_map <- background_map[background_map$SYMBOL %in% blue_genes, , drop = FALSE]
background_entrez <- unique(as.character(background_map$ENTREZID))
blue_entrez <- intersect(unique(as.character(blue_map$ENTREZID)), background_entrez)

write_tsv(background_map, file.path(enrich_dir, "WGCNA_background_SYMBOL_ENTREZ.tsv"))
write_tsv(blue_map, file.path(enrich_dir, "BLUE_SYMBOL_ENTREZ.tsv"))
coverage <- data.frame(
  Item = c("Paired_donors", "WGCNA_genes_SYMBOL", "Blue_genes_SYMBOL",
           "Mapped_background_SYMBOL", "Mapped_blue_SYMBOL",
           "Mapped_background_ENTREZ", "Mapped_blue_ENTREZ"),
  Count = c(nrow(pair_tab), length(wgcna_genes), length(blue_genes),
            length(unique(background_map$SYMBOL)), length(unique(blue_map$SYMBOL)),
            length(background_entrez), length(blue_entrez))
)
write_tsv(coverage, file.path(enrich_dir, "enrichment_gene_counts.tsv"))
print(coverage, row.names = FALSE)
if (length(blue_entrez) < 10L) {
  stop("blue Entrez 10 SYMBOL（ ID）.")
}

# --------------------------- 7. Enrichment plots -----------------------------
parse_ratio <- function(s) {
  vapply(strsplit(as.character(s), "/", fixed = TRUE), function(z) {
    if (length(z) != 2L) return(NA_real_)
    as.numeric(z[1]) / as.numeric(z[2])
  }, numeric(1))
}
wrap_one <- function(x, width = 49L) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

plot_enrich <- function(tab, title, filepath_stem, n_show = TOP_TERMS) {
  tab <- tab[is.finite(tab$p.adjust) & tab$p.adjust < FDR_CUTOFF, , drop = FALSE]
  if (!nrow(tab)) {
    message("[No significant pathways] ", title, " (BH FDR < ", FDR_CUTOFF, ")")
    return(invisible(NULL))
  }
  tab <- tab[order(tab$p.adjust, -tab$Count), , drop = FALSE]
  tab <- head(tab, n_show)
  tab$GeneRatioValue <- parse_ratio(tab$GeneRatio)
  tab$negLog10FDR <- -log10(pmax(tab$p.adjust, .Machine$double.xmin))
  tab$term <- vapply(tab$Description, wrap_one, character(1))
  tab$term <- make.unique(tab$term)
  tab$term <- factor(tab$term, levels = rev(tab$term))
  tab <- tab[is.finite(tab$GeneRatioValue), , drop = FALSE]
  if (!nrow(tab)) return(invisible(NULL))

  p <- ggplot2::ggplot(tab, ggplot2::aes(x = GeneRatioValue, y = term)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = negLog10FDR), alpha = 0.94) +
    ggplot2::scale_color_gradient(low = "#9CC1D8", high = "#2F5889",
                                   name = expression(-log[10]("BH FDR"))) +
    ggplot2::scale_size_continuous(range = c(3.4, 8.5), name = "Gene count") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.10))) +
    ggplot2::labs(title = title, subtitle = "GSE145989 | WGCNA blue module",
                  x = "Gene ratio", y = NULL) +
    pretty_theme +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 10))
  height <- max(5.3, 0.38 * nrow(tab) + 1.6)
  ggplot2::ggsave(paste0(filepath_stem, ".pdf"), p,
                  width = 10.2, height = height, bg = "white",
                  device = "pdf", useDingbats = FALSE)
  ggplot2::ggsave(paste0(filepath_stem, ".png"), p,
                  width = 10.2, height = height, dpi = 400, bg = "white")
  invisible(p)
}

# ------------------------- 8. GO: BP, CC, MF ----------------------------------
all_GO <- list()
for (ont in c("BP", "CC", "MF")) {
  message("[GO] Analyzing ontology: ", ont)
  fit <- clusterProfiler::enrichGO(
    gene = blue_entrez, universe = background_entrez,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
    ont = ont, pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 500,
    readable = TRUE
  )
  tab <- as.data.frame(fit)
  if (!nrow(tab)) {
    message("[GO] ", ont, "Check source data or configuration. ")
    next
  }
  tab$Ontology <- ont
  all_GO[[ont]] <- tab
  write_tsv(tab, file.path(enrich_dir, paste0("GO_", ont, "_BLUE_all.tsv")))
  significant <- tab[!is.na(tab$p.adjust) & tab$p.adjust < FDR_CUTOFF, , drop = FALSE]
  write_tsv(significant,
            file.path(enrich_dir, paste0("GO_", ont, "_BLUE_FDR05.tsv")))
  plot_enrich(tab, paste0("GO ", ont, " enrichment"),
              file.path(enrich_dir, paste0("GO_", ont, "_BLUE_dotplot")))
}
if (length(all_GO)) {
  write_tsv(do.call(rbind, all_GO), file.path(enrich_dir, "GO_ALL_BLUE_full.tsv"))
}

# ----------------------------- 9. KEGG ---------------------------------------
message("[KEGG] Querying KEGG pathways (internet may be required) ...")
kegg_error <- NULL
fit_kegg <- tryCatch(
  clusterProfiler::enrichKEGG(
    gene = blue_entrez, universe = background_entrez,
    organism = "hsa", keyType = "ncbi-geneid", pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 500
  ),
  error = function(e) {
    kegg_error <<- conditionMessage(e)
    NULL
  }
)
if (!is.null(kegg_error)) {
  writeLines(c("KEGG query failed:", kegg_error,
               "Please verify access to KEGG and retry this section."),
             file.path(enrich_dir, "KEGG_query_error.txt"))
  warning("KEGG GO ", kegg_error)
} else {
  kegg_tab <- as.data.frame(fit_kegg)
  if (nrow(kegg_tab)) {
    write_tsv(kegg_tab, file.path(enrich_dir, "KEGG_BLUE_all.tsv"))
    write_tsv(kegg_tab[!is.na(kegg_tab$p.adjust) &
                             kegg_tab$p.adjust < FDR_CUTOFF, , drop = FALSE],
              file.path(enrich_dir, "KEGG_BLUE_FDR05.tsv"))
    plot_enrich(kegg_tab, "KEGG pathway enrichment",
                file.path(enrich_dir, "KEGG_BLUE_dotplot"))
  } else {
    message("[KEGG] No pathways returned for mapped blue genes")
    writeLines("No KEGG pathways were returned; no plot drawn.",
               file.path(enrich_dir, "KEGG_no_results.txt"))
  }
}

# --------------------------- 10. Reproducibility ------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("\n========= FINISHED =========\n")
cat("Module: BLUE | ", length(blue_genes), " symbols | ",
    length(blue_entrez), " Entrez IDs\n", sep = "")
cat("Enrichment background: ", length(wgcna_genes), " WGCNA symbols\n", sep = "")
cat("Output: ", normalizePath(OUT_DIR, winslash = "/"), "\n", sep = "")
cat("Paired module FDR was reused; it was NOT recomputed with unpaired tests.\n")


## =========================================================
## FINAL VERSION
## Blue module GO BP plot
## =========================================================

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

## ---------------------------
## ---------------------------
go_bp <- read_tsv("GO_BP_BLUE_all.tsv", show_col_types = FALSE)

go_bp <- go_bp %>%
  filter(Ontology == "BP")

ratio_to_num <- function(x) {
  sapply(strsplit(x, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]))
}

go_bp <- go_bp %>%
  mutate(
    GeneRatio_num = ratio_to_num(GeneRatio),
    neglog10FDR = -log10(p.adjust)
  )

## ---------------------------
## ---------------------------
target_order <- c(
  "cellular response to cytokine stimulus",
  "apoptotic signaling pathway",
  "cytokine-mediated signaling pathway",
  "inflammatory response",
  "regulation of cell-cell adhesion",
  "stress-activated MAPK cascade",
  "canonical NF-kappaB signal transduction",
  "response to tumor necrosis factor",
  "positive regulation of cell migration",
  "response to hypoxia",
  "regulation of blood vessel endothelial cell migration",
  "response to reactive oxygen species",
  "protein modification by small protein conjugation or removal",
  "post-translational protein modification",
  "protein ubiquitination"
)

## ---------------------------
## ---------------------------
plot_df <- go_bp %>%
  filter(Description %in% target_order) %>%
  distinct(Description, .keep_all = TRUE)

missing_terms <- setdiff(target_order, plot_df$Description)
if (length(missing_terms) > 0) {
  message("Check source data or configuration. ")
  print(missing_terms)
}

present_order <- target_order[target_order %in% plot_df$Description]

plot_df <- plot_df %>%
  mutate(
    Description = factor(Description, levels = rev(present_order)),
    Description_wrap = factor(
      stringr::str_wrap(as.character(Description), width = 34),
      levels = stringr::str_wrap(rev(present_order), width = 34)
    )
  )

## ---------------------------
## ---------------------------
p_go_final <- ggplot(
  plot_df,
  aes(
    x = GeneRatio_num,
    y = Description_wrap,
    size = Count,
    color = neglog10FDR
  )
) +
  geom_point(alpha = 0.95) +
  scale_size_continuous(
    name = "Gene count",
    range = c(3.5, 8.5),
    breaks = pretty(plot_df$Count, n = 4)
  ) +
  scale_color_gradient2(
    name = expression(-log[10](FDR)),
    low = "#8FB6D9",
    mid = "#F3F3F3",
    high = "#C85D57",
    midpoint = median(plot_df$neglog10FDR, na.rm = TRUE)
  ) +
  scale_x_continuous(
    name = "Gene ratio",
    labels = number_format(accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  ylab(NULL) +
  ggtitle("GO biological process") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(color = "#D9D9D9", linewidth = 0.6),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 0.8),
    
    axis.text.x = element_text(color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 10, face = "bold"),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_blank(),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14,
      color = "black"
    ),
    
    legend.title = element_text(face = "bold", size = 11, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    
    plot.margin = margin(8, 14, 8, 8)
  ) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  )

print(p_go_final)

## ---------------------------
## ---------------------------
ggsave("GO_BP_BLUE_final_15terms_reference_style.pdf",
       p_go_final, width = 8.0, height = 6.3)

ggsave("GO_BP_BLUE_final_15terms_reference_style.png",
       p_go_final, width = 8.0, height = 6.3, dpi = 600)

write_tsv(
  plot_df %>%
    mutate(Description = as.character(Description_wrap)) %>%
    select(ID, Description, GeneRatio, GeneRatio_num, Count, pvalue, p.adjust, neglog10FDR),
  "GO_BP_BLUE_final_15terms_reference_style.tsv"
)




}

# =============================================================================
# STAGE figure1_ml
# =============================================================================
PIPELINE$figure1_ml <- function() {
# scripts/03_machine_learning.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.
validate_ml_inputs(repo_input("ml_train"), repo_input("ml_test"), repo_input("ml_methods"))
# =============================================================================
# ONE-ENTRY-POINT COMPLETE RE-ANALYSIS: donor-matched GSE145989 -> GSE18995
# Revised 2026-10-21. THIS SCRIPT DOES NOT READ ANY PREVIOUS RANKING OR PREDICTIONS.
# 53 donor pairs (106 biopsies) for CV/selection; 16 independently mapped donor
# pairs (32 biopsies) only for external evaluation AFTER model selection.
# In every outer CV fold, preprocessing, supervised feature selection and
# supported inner tuning are performed on TRAINING DONORS ONLY.
# Every selected algorithm is refitted on the 53 donor pairs to produce an
# exploratory GSE18995 AUC column, but GSE18995 does NOT select the winner.
# Results are conditional on the user-provided candidate gene panel: if it was
# selected using all 53 donors, OOF CV is not a fully nested upstream-DEG CV.
# =============================================================================
options(stringsAsFactors = FALSE)

# ===================== CONFIGURE ONLY THIS SECTION ===========================
# In RStudio, use Session > Set Working Directory > To Source File Location,
# or setwd("YOUR/FOLDER") first. Place these scripts + matrices + pair maps here.
TRAIN_FILE <- repo_input("ml_train")
TEST_FILE <- repo_input("ml_test")
# Supply the actual 53-donor crosswalk from the final sample audit.
# It covers all 106 samples. Read the original-stage review notices below.
TRAIN_PAIR_FILE <- file.path(REPO_ROOT,"metadata/GSE145989_pairs.csv")
TEST_PAIR_FILE <- file.path(REPO_ROOT,"metadata/GSE18995_pairs.csv")
HELPER_FILE <- SCRIPT_FILE
METHOD_FILE <- repo_input("ml_methods")
OUTPUT_ROOT <- repo_output("Figure1/ML")
RANDOM_SEED <- 123456L
OUTER_FOLDS <- 5L
MIN_FEATURES <- 2L
MAX_METHODS <- Inf  # QUICK CHECK: 10L; FINAL FULL RUN: Inf
PAIRED_BOOTSTRAP_B <- 1000L
CHECKPOINT <- TRUE
RUN_EXTERNAL_FOR_ALL_SUCCESSFUL <- TRUE
MAKE_FIGURES <- TRUE
# The historical audit flags Toronto_575 and Toronto_605 with stage 'R' rather
# than explicitly '2H' in CEL names. The user's prior screening retained them.
# We do NOT silently claim that the original CEL name proves 2h for these two.
# Their original stage flags remain visible in the bundled CSV and audit note.
# ============================================================================

if (!file.exists(HELPER_FILE)) stop("Check source data or configuration. ", HELPER_FILE,
   "Check source data or configuration. ",getwd(),"\n ZIP ")
load_ml_helpers()
required <- c(TRAIN_FILE,TEST_FILE,TEST_PAIR_FILE,METHOD_FILE)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Check source data or configuration. ",paste(missing,collapse="\n"),
                         "Check source data or configuration. ",getwd())
if (!file.exists(TRAIN_PAIR_FILE)) stop("Check source data or configuration. ",TRAIN_PAIR_FILE)
# Report—not conceal—ambiguities inherited from the original CEL audit.
stage_audit <- utils::read.csv(TRAIN_PAIR_FILE,check.names=FALSE,stringsAsFactors=FALSE,
                              fileEncoding="UTF-8-BOM")
if ("stage_review_required" %in% names(stage_audit)) {
  flagged <- stage_audit$donor_id[as.character(stage_audit$stage_review_required) %in% c("True","TRUE","1")]
  if (length(flagged)) warning("GEO CEL : ",
    paste(flagged,collapse=", "),
    ".  GSM R 2H. ")
}
if (!requireNamespace("pROC",quietly=TRUE)) stop(" pROC, install.packages('pROC')")

read_input <- function(path,cohort,n_pairs) {
  dat <- utils::read.delim(path,header=TRUE,row.names=1,sep="\t",
                           check.names=FALSE,stringsAsFactors=FALSE,
                           quote="",comment.char="")
  if (anyDuplicated(rownames(dat)) || !"Type" %in% names(dat))
    stop(path,":   Type  ")
  lab <- trimws(as.character(dat$Type))
  y <- ifelse(tolower(lab) %in% c("0","cit","control"),0L,
         ifelse(tolower(lab) %in% c("1","r2h","treat","hr"),1L,NA_integer_))
  if (anyNA(y)) stop(path,":   Type: ",paste(unique(lab[is.na(y)]),collapse=","))
  id <- rownames(dat)
  if (!all(startsWith(id,paste0(cohort,"_"))) ||
      length(id)!=2*n_pairs ||
      !identical(as.integer(table(factor(y,levels=0:1))),rep(as.integer(n_pairs),2L)))
    stop(path,"Check source data or configuration. ",cohort,"Check source data or configuration. ",n_pairs,"Check source data or configuration. ",2*n_pairs,"Check source data or configuration. ")
  feat <- setdiff(names(dat),"Type")
  X <- as.data.frame(dat[,feat,drop=FALSE],check.names=FALSE)
  X[] <- lapply(X,function(v) suppressWarnings(as.numeric(as.character(v))))
  if (anyNA(X) || any(!is.finite(as.matrix(X)))) stop(path,":  /NA/Inf ")
  list(X=X,y=as.integer(y),id=id)
}
train <- read_input(TRAIN_FILE,"GSE145989",53L)
test  <- read_input(TEST_FILE,"GSE18995",16L)
train_donor <- read_verified_pairs(TRAIN_PAIR_FILE,train$id,train$y,53L)
test_donor  <- read_verified_pairs(TEST_PAIR_FILE,test$id,test$y,16L)
if (length(intersect(train$id,test$id))) stop("Check source data or configuration. ")
if (length(intersect(train_donor,test_donor))) stop("Check source data or configuration. ")
allgenes <- intersect(colnames(train$X),colnames(test$X))
if (length(allgenes)<2L) stop("Check source data or configuration. ")
train$X <- train$X[,allgenes,drop=FALSE]
test$X  <- test$X[,allgenes,drop=FALSE]
cat("\n Train ",nrow(train$X)," samples / ",length(unique(train_donor)),
    " donors; Test ",nrow(test$X)," samples / ",length(unique(test_donor)),
    " donors; ",length(allgenes)," overlapping candidate genes\n",sep="")

mt <- utils::read.delim(METHOD_FILE,sep="\t",header=TRUE,
                        check.names=FALSE,stringsAsFactors=FALSE)
if (!"Model" %in% names(mt)) stop(METHOD_FILE,"   Model  ")
methods <- unique(gsub("[[:space:]]|-","",as.character(mt$Model)))
methods <- methods[!is.na(methods)&nzchar(methods)]
if (is.finite(MAX_METHODS)) methods <- head(methods,as.integer(MAX_METHODS))
if (length(methods)<2L) stop("Check source data or configuration. ")
message("Check source data or configuration. ",length(methods),"Check source data or configuration. ")

# Version-locked resume folders: no mixing different inputs, helper, settings.
inputs <- c(TRAIN_FILE,TEST_FILE,TRAIN_PAIR_FILE,TEST_PAIR_FILE,METHOD_FILE,
            HELPER_FILE,"01_RUN_FROM_SCRATCH_DONOR_GROUPED_ML.R")
inputs <- inputs[file.exists(inputs)]
hash <- paste(substr(unname(tools::md5sum(inputs)),1L,8L),collapse="-")
key <- paste0(hash,"_seed",RANDOM_SEED,"_k",OUTER_FOLDS,"_min",MIN_FEATURES,
              "_n",length(methods))
RUN_DIR <- file.path(OUTPUT_ROOT,paste0("run_",key))
dir.create(RUN_DIR,showWarnings=FALSE,recursive=TRUE)
CV_CACHE <- file.path(RUN_DIR,"checkpoints_CV")
EXT_CACHE <- file.path(RUN_DIR,"checkpoints_external")
dir.create(CV_CACHE,showWarnings=FALSE,recursive=TRUE)
dir.create(EXT_CACHE,showWarnings=FALSE,recursive=TRUE)
cat("Check source data or configuration. ",normalizePath(RUN_DIR),"\n",sep="")
writeLines(c("Full from-scratch donor-grouped ML screening", paste("Date",Sys.Date()),
 "Train GSE145989 53 pairs, external GSE18995 16 pairs",
 "CIT=0, R2H=1. Models selected only by pooled donor-grouped OOF AUC.",
 "Candidate panel is predeclared; upstream supervised gene-panel selection is not nested.",
 "RF-selector uses fixed in-fold positive permutation importance (not original var.select).",
 "plsRglm fixed nt<=2, no row-CV tuning (version-dependent predictions).",
 "SVM uses oriented sigmoid of decision margin as an AUC score, NOT calibrated probability.",
 "External AUC for non-winners is exploratory. Average train/external AUC is display ONLY.",
 "Checkpoints tied to MD5 fingerprints of input matrices, maps, code, method list."),
 file.path(RUN_DIR,"00_analysis_manifest.txt"))
utils::write.table(data.frame(Gene=allgenes),file.path(RUN_DIR,"00_model_input_genes.tsv"),
                   sep="\t",quote=FALSE,row.names=FALSE)

# One donor is indivisible; both biopsies always in the same outer fold.
foldid <- make_donor_folds(train_donor,k=OUTER_FOLDS,seed=RANDOM_SEED)
fa <- data.frame(Sample=train$id,Donor=as.character(train_donor),
                 Condition=ifelse(train$y==0,"CIT","R2H"),Fold=foldid)
if (!all(vapply(split(fa$Fold,fa$Donor),function(x)length(unique(x))==1L,logical(1))))
  stop("FATAL: donor split across CV folds")
utils::write.table(fa,file.path(RUN_DIR,"01_verified_outer_donor_folds.tsv"),
                   sep="\t",quote=FALSE,row.names=FALSE)
cat("Check source data or configuration. ",paste(as.integer(table(factor(fa$Fold,levels=1:OUTER_FOLDS))/2),collapse="/"),"Check source data or configuration. ",sep="")

run_cv <- function(method,index) {
  prob <- rep(NA_real_,length(train$y))
  fold_auc <- rep(NA_real_,OUTER_FOLDS)
  genes_by_fold <- character(OUTER_FOLDS)
  for (v in seq_len(OUTER_FOLDS)) {
    ii <- which(foldid!=v); jj <- which(foldid==v)
    # selector and classifier rebuilt within EACH training fold; all tuning uses donor IDs.
    mod <- fit_ml_pipeline(method,train$X[ii,,drop=FALSE],train$y[ii],
             train_donor[ii],seed=RANDOM_SEED+50L*index+v,
             min_features=MIN_FEATURES)
    prob[jj] <- predict_ml_pipeline(mod,train$X[jj,,drop=FALSE])
    fold_auc[v] <- auc_fixed_direction(train$y[jj],prob[jj])
    genes_by_fold[v] <- paste(mod$selected_genes,collapse=";")
  }
  if (any(!is.finite(prob))) stop("CV contains missing/non-finite OOF scores")
  list(Status="OK", Method=method, CV_AUC=auc_fixed_direction(train$y,prob),
       Mean_Fold_AUC=mean(fold_auc),SD_Fold_AUC=stats::sd(fold_auc),
       Fold_AUC=fold_auc,GenesByFold=genes_by_fold,OOF=prob,Reason="")
}
status_record <- function(method,why) list(
  Status=if (grepl("^SKIP:|unavailable|Unsupported",why)) "SKIPPED" else "FAILED",
  Method=method,CV_AUC=NA_real_,Mean_Fold_AUC=NA_real_,SD_Fold_AUC=NA_real_,
  Fold_AUC=rep(NA_real_,OUTER_FOLDS),GenesByFold=rep("",OUTER_FOLDS),
  OOF=rep(NA_real_,length(train$y)),Reason=why)

cv_results <- vector("list",length(methods))
for (i in seq_along(methods)) {
  ck <- file.path(CV_CACHE,sprintf("model_%03d.rds",i))
  if (CHECKPOINT && file.exists(ck)) {
    x <- tryCatch(readRDS(ck),error=function(e) NULL)
    if (!is.null(x) && identical(x$Method,methods[i]) &&
        length(x$OOF)==length(train$y)) {
      cv_results[[i]] <- x
      cat(sprintf("[CV %d/%d] resumed %s: %s\n",i,length(methods),methods[i],x$Status))
      next
    }
  }
  cat(sprintf("[CV %d/%d] %s ... ",i,length(methods),methods[i]))
  t0 <- proc.time()[["elapsed"]]
  x <- tryCatch(run_cv(methods[i],i),error=function(e)
    status_record(methods[i],conditionMessage(e)))
  x$Elapsed_s <- proc.time()[["elapsed"]]-t0
  saveRDS(x,ck)
  cv_results[[i]] <- x
  cat(x$Status,if(is.finite(x$CV_AUC))sprintf(" AUC=%.4f",x$CV_AUC) else "",
      " (",round(x$Elapsed_s,1),"s)","\n",sep="")
  if (x$Status!="OK") message("Check source data or configuration. ",x$Reason)
}

cv_table <- do.call(rbind,lapply(cv_results,function(v) data.frame(
  Method=v$Method,Status=v$Status,CV_AUC=v$CV_AUC,
  Mean_Fold_AUC=v$Mean_Fold_AUC,SD_Fold_AUC=v$SD_Fold_AUC,
  Fold1=v$Fold_AUC[1],Fold2=v$Fold_AUC[2],Fold3=v$Fold_AUC[3],
  Fold4=v$Fold_AUC[4],Fold5=v$Fold_AUC[5],
  Runtime_seconds=if(is.null(v$Elapsed_s))NA_real_ else v$Elapsed_s,
  Reason=v$Reason,check.names=FALSE)))
cv_table <- cv_table[order(!is.finite(cv_table$CV_AUC),-cv_table$CV_AUC,
                            -cv_table$Mean_Fold_AUC,cv_table$Method),,drop=FALSE]
rownames(cv_table) <- NULL
CV_FILE <- file.path(RUN_DIR,"02_grouped_CV_model_ranking.tsv")
utils::write.table(cv_table,CV_FILE,sep="\t",quote=FALSE,row.names=FALSE,na="")
utils::write.table(cv_table[cv_table$Status!="OK",,drop=FALSE],
   file.path(RUN_DIR,"03_failed_or_skipped_models.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="")

valid <- cv_table[cv_table$Status=="OK" & is.finite(cv_table$CV_AUC),,drop=FALSE]
if (!nrow(valid)) stop("Check source data or configuration. ")
WINNER <- valid$Method[1L]  # LOCKED before any GSE18995 model evaluation.
win <- cv_results[[match(WINNER,methods)]]
cat("\n*****  TRAIN CV : ",WINNER,
    "; OOF AUC=",sprintf("%.4f",win$CV_AUC)," *****\n",sep="")
writeLines(c(paste("Selected before external validation:",WINNER),
             paste("Training grouped OOF AUC:",sprintf("%.6f",win$CV_AUC)),
             paste("Successfully audited:",nrow(valid),"of",length(methods)),
             "External cohort was NOT used for model selection."),
           file.path(RUN_DIR,"04_SELECTED_MODEL_BEFORE_EXTERNAL.txt"))
train_pred <- data.frame(Sample=train$id,Donor=train_donor,Type=train$y,Fold=foldid,
                         Score_R2H=win$OOF,Predicted=as.integer(win$OOF>=0.5))
utils::write.table(train_pred,file.path(RUN_DIR,"04_selected_model_train_OOF_predictions.tsv"),
                   sep="\t",quote=FALSE,row.names=FALSE)

# Refit EVERY CV-successful model on all 53 training donors, without reading
# external labels until after training is frozen. Do not promote by GSE AUC.
ext_results <- vector("list",nrow(cv_table))
all_external <- list()
selected_frozen <- NULL
for (i in seq_len(nrow(cv_table))) {
  nam <- cv_table$Method[i]
  if (cv_table$Status[i]!="OK") {
    ext_results[[i]] <- list(Method=nam,Status="NOT_CV_SUCCESS",AUC=NA_real_,
                             Genes="",Reason="Failed/skipped training CV")
    next
  }
  if (!RUN_EXTERNAL_FOR_ALL_SUCCESSFUL && nam!=WINNER) {
    ext_results[[i]] <- list(Method=nam,Status="NOT_REQUESTED",AUC=NA_real_,
                             Genes="",Reason="External evaluation limited to winner")
    next
  }
  ck <- file.path(EXT_CACHE,sprintf("method_%03d.rds",i))
  e <- NULL
  if (CHECKPOINT && file.exists(ck)) {
    c <- tryCatch(readRDS(ck),error=function(er) NULL)
    if (!is.null(c) && identical(c$Method,nam) && length(c$Scores)==length(test$y)) e<-c
  }
  if (is.null(e)) {
    cat(sprintf("[EXTERNAL %d/%d] %s ... ",i,nrow(cv_table),nam))
    t0<-proc.time()[["elapsed"]]
    e <- tryCatch({
      # Frozen full 53-pair model; test labels NOT used in fit_ml_pipeline.
      trained <- fit_ml_pipeline(nam,train$X,train$y,train_donor,
          seed=RANDOM_SEED+1234L,min_features=MIN_FEATURES)
      scores <- predict_ml_pipeline(trained,test$X)
      if (nam==WINNER) {
        selected_frozen <- trained
        saveRDS(trained,file.path(RUN_DIR,"05_selected_frozen_model.rds"))
      }
      list(Method=nam,Status="OK",AUC=auc_fixed_direction(test$y,scores),
           Scores=scores,Genes=paste(trained$selected_genes,collapse=";"),Reason="")
    },error=function(er) list(Method=nam,Status="FAILED",AUC=NA_real_,
      Scores=rep(NA_real_,length(test$y)),Genes="",Reason=conditionMessage(er)))
    e$Elapsed_s <- proc.time()[["elapsed"]]-t0
    saveRDS(e,ck)
    cat(e$Status,if(is.finite(e$AUC))sprintf(" AUC=%.4f",e$AUC) else "",
        "\n",sep="")
    if (e$Status!="OK") message("Check source data or configuration. ",e$Reason)
  }
  if (e$Status=="OK") {
    all_external[[length(all_external)+1L]] <- data.frame(
       Method=nam,Sample=test$id,Donor=test_donor,Type=test$y,
       Score_R2H=e$Scores,Predicted=as.integer(e$Scores>=0.5),stringsAsFactors=FALSE)
  }
  ext_results[[i]]<-e
}
# Resume case: the winning model's RDS might already exist from prior run.
if (!file.exists(file.path(RUN_DIR,"05_selected_frozen_model.rds"))) {
  message("Check source data or configuration. ")
  selected_frozen <- fit_ml_pipeline(WINNER,train$X,train$y,train_donor,
                                    seed=RANDOM_SEED+1234L,min_features=MIN_FEATURES)
  saveRDS(selected_frozen,file.path(RUN_DIR,"05_selected_frozen_model.rds"))
}

cv_table$GSE18995_AUC <- vapply(ext_results,`[[`,numeric(1),"AUC")
cv_table$External_Status <- vapply(ext_results,`[[`,character(1),"Status")
cv_table$External_Reason <- vapply(ext_results,`[[`,character(1),"Reason")
cv_table$Final_Model_Genes <- vapply(ext_results,`[[`,character(1),"Genes")
cv_table$Display_Mean_AUC <- ifelse(is.finite(cv_table$CV_AUC) &
  is.finite(cv_table$GSE18995_AUC),
  (cv_table$CV_AUC+cv_table$GSE18995_AUC)/2,NA_real_)
cv_table$CV_Rank <- NA_integer_
cv_table$CV_Rank[which(cv_table$Status=="OK" & is.finite(cv_table$CV_AUC))] <-
  seq_len(sum(cv_table$Status=="OK" & is.finite(cv_table$CV_AUC)))
FULL_RANK_FILE <- file.path(RUN_DIR,"06_CV_ranked_with_ALL_GSE18995_AUC.tsv")
utils::write.table(cv_table,FULL_RANK_FILE,sep="\t",quote=FALSE,row.names=FALSE,na="")
if (length(all_external)) {
  expred <- do.call(rbind,all_external)
  utils::write.table(expred,file.path(RUN_DIR,"07_ALL_models_GSE18995_predictions.tsv"),
                     sep="\t",quote=FALSE,row.names=FALSE)
}

winner_external <- ext_results[[match(WINNER,cv_table$Method)]]
if (is.null(winner_external) || !identical(winner_external$Status,"OK")) {
  warning("Check source data or configuration. ")
} else {
  etab <- data.frame(Sample=test$id,Donor=test_donor,Type=test$y,
                     Score_R2H=winner_external$Scores,
                     Predicted=as.integer(winner_external$Scores>=0.5))
  utils::write.table(etab,file.path(RUN_DIR,"08_SELECTED_model_GSE18995_predictions.tsv"),
                     sep="\t",quote=FALSE,row.names=FALSE)
  external_auc <- winner_external$AUC
  cv_ci <- paired_auc_bootstrap(train$y,win$OOF,train_donor,
                               B=PAIRED_BOOTSTRAP_B,seed=RANDOM_SEED+321L)
  ext_ci <- paired_auc_bootstrap(test$y,winner_external$Scores,test_donor,
                                B=PAIRED_BOOTSTRAP_B,seed=RANDOM_SEED+654L)
  perf <- data.frame(Dataset=c("GSE145989 grouped 5-fold CV","GSE18995 independent test"),
      Donors=c(53L,16L),Samples=c(106L,32L),
      AUC=c(win$CV_AUC,external_auc),
      ClusterBootstrap_CI_Low=c(cv_ci[1],ext_ci[1]),
      ClusterBootstrap_CI_High=c(cv_ci[2],ext_ci[2]))
  utils::write.table(perf,file.path(RUN_DIR,"09_SELECTED_model_ROC_summary.tsv"),
                     sep="\t",quote=FALSE,row.names=FALSE)
  cat("\n===== FINAL WINNER (train CV selection only) =====\n",
      "Model = ",WINNER,"\nGrouped OOF AUC = ",sprintf("%.4f",win$CV_AUC),
      "\nGSE18995 external AUC = ",sprintf("%.4f",external_auc),
      "\nNote: non-winner external AUCs are descriptive, not for re-ranking.\n",sep="")
}

sink(file.path(RUN_DIR,"99_sessionInfo.txt"));print(utils::sessionInfo());sink()
# Run the embedded figure1_ml_plots stage after this stage.

cat("Check source data or configuration. ",normalizePath(RUN_DIR),"\n",sep="")

saveRDS(list(RUN_DIR=RUN_DIR),file.path(repo_output("Figure1/ML"),"selected_run_pointer.rds"))

}

# =============================================================================
# STAGE figure1_ml_plots
# =============================================================================
PIPELINE$figure1_ml_plots <- function() {
# scripts/03_machine_learning_plots.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.
RUN_DIR <- repo_ml_run()
# =============================================================================
# FIGURE 1I-J (+ K) | ORIGINAL REFERENCE STYLE, donor-grouped analysis
# Figure I: ComplexHeatmap SimpleHeatmap style of refer.ML.R.
# Figure J: ggplot confusion-matrix tile style of original publication script.
# Figure K: ROC in the same publication-figure palette.
# Inputs: ONLY the outputs of 01_RUN_FROM_SCRATCH_DONOR_GROUPED_ML.R
# =============================================================================
options(stringsAsFactors = FALSE)
# Cosmetic adjustments only: model predictions/AUC/ranking remain unchanged.
# All labels (I bar/cells, J metrics, K AUC) are displayed with three decimals.
LABEL_DIGITS <- 3L
HEAT_CELL_LINEWIDTH <- 1.30  # reference: 1.0
CONFUSION_TILE_LINEWIDTH <- 1.65  # reference: 1.3
ROC_LINEWIDTH <- 1.65  # a little stronger than previous 1.55
BAR_WIDTH_CM <- 2.30   # shorter right-hand AUC bars (previous 3.20)
BAR_LABEL_FONT_PT <- 8.2
HEAT_LABEL_FONT_PT <- 8.0
# All text in panels I/J/K will be black and bold. Curve/tile colors remain unchanged.

if (!exists("RUN_DIR", inherits=TRUE)) {
  root <- "ML_PAIR_FROM_SCRATCH_results"
  runs <- if (dir.exists(root)) list.dirs(root,recursive=FALSE,full.names=TRUE) else character()
  good <- runs[file.exists(file.path(runs,"06_CV_ranked_with_ALL_GSE18995_AUC.tsv"))]
  if (!length(good)) stop(" 01_RUN_FROM_SCRATCH_DONOR_GROUPED_ML.R")
  if (length(good) > 1L) {
    # Use the newest completed ranking; user can override RUN_DIR explicitly.
    mt <- file.info(file.path(good,"06_CV_ranked_with_ALL_GSE18995_AUC.tsv"))$mtime
    good <- good[order(mt,decreasing=TRUE)]
    message("Check source data or configuration. ", good[1],
            ".  RUN_DIR.")
  }
  RUN_DIR <- good[1]
}
if (!dir.exists(RUN_DIR)) stop("Check source data or configuration. ",RUN_DIR)
OUT <- file.path(RUN_DIR,"FIGURE1_IJK")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

rank_file <- file.path(RUN_DIR,"06_CV_ranked_with_ALL_GSE18995_AUC.tsv")
if (!file.exists(rank_file)) stop("Check source data or configuration. ",rank_file)
r <- utils::read.delim(rank_file,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
expected <- c("Method","Status","CV_AUC","Mean_Fold_AUC","GSE18995_AUC","Display_Mean_AUC")
if (!all(expected %in% names(r))) stop("Check source data or configuration. ",paste(setdiff(expected,names(r)),collapse=","))
for (nm in c("CV_AUC","Mean_Fold_AUC","GSE18995_AUC","Display_Mean_AUC"))
  r[[nm]] <- suppressWarnings(as.numeric(r[[nm]]))
vis <- r[r$Status=="OK" & is.finite(r$CV_AUC) & is.finite(r$GSE18995_AUC) &
         is.finite(r$Display_Mean_AUC),,drop=FALSE]
vis <- vis[order(-vis$CV_AUC,-vis$Mean_Fold_AUC,vis$Method),,drop=FALSE]
if (!nrow(vis)) stop(" CV GSE18995 AUC ")

# Set DRAW_FIGURE_I <- FALSE before source(...) to regenerate only J/K.
if (!exists("DRAW_FIGURE_I", inherits=TRUE)) DRAW_FIGURE_I <- TRUE
if (isTRUE(DRAW_FIGURE_I)) {
# -------- FIGURE I: preserve original ComplexHeatmap layout -----------------
if (!requireNamespace("ComplexHeatmap",quietly=TRUE)) stop(" ComplexHeatmap: BiocManager::install('ComplexHeatmap')")
if (!requireNamespace("circlize",quietly=TRUE)) stop(" circlize")
M <- cbind("Train"=vis$CV_AUC, "GSE18995"=vis$GSE18995_AUC)
rownames(M) <- vis$Method
# Average is ONLY for reference's steelblue display bar, not selecting winner.
avg_AUC <- vis$Display_Mean_AUC
CohortCol <- c("Train"="#CB5746", "GSE18995"="#4195C1")
col_ha <- ComplexHeatmap::columnAnnotation(
   Cohort=colnames(M), col=list(Cohort=CohortCol),
   show_annotation_name=FALSE,
   annotation_legend_param=list(
     Cohort=list(labels_gp=grid::gpar(col="black",fontface="bold",fontsize=9),
                 title_gp=grid::gpar(col="black",fontface="bold",fontsize=10))
   ))
# IMPORTANT: anno_barplot(add_numbers=TRUE) prints the full floating-point
# AUC, not exactly 3 decimals. Hide automatic values and overlay exact 3-decimal labels CENTERED
# INSIDE shorter steelblue bars. Only width/label position changes; AUC unchanged.
row_ha <- ComplexHeatmap::rowAnnotation(
  bar=ComplexHeatmap::anno_barplot(
     avg_AUC,
     bar_width=0.80,
     border=FALSE,
     gp=grid::gpar(fill="steelblue",col=NA),
     baseline=0,
     ylim=c(0,1),
     add_numbers=FALSE,
     axis_param=list(at=c(0,0.5,1),labels=c("0","0.5","1"),
                     labels_rot=0,
                     gp=grid::gpar(col="black",fontface="bold",fontsize=8)),
     width=grid::unit(BAR_WIDTH_CM,"cm")
  ),
  show_annotation_name=FALSE
)
hm <- ComplexHeatmap::Heatmap(
  as.matrix(M), name="AUC",right_annotation=row_ha,top_annotation=col_ha,
  col=c("#4195C1","#FFFFFF","#CB5746"), # ORIGINAL reference palette
  heatmap_legend_param=list(
    labels_gp=grid::gpar(col="black",fontface="bold",fontsize=9),
    title_gp=grid::gpar(col="black",fontface="bold",fontsize=10)
  ),
  rect_gp=grid::gpar(col="#172E3D",lwd=HEAT_CELL_LINEWIDTH),
  cluster_columns=FALSE,cluster_rows=FALSE,
  show_column_names=FALSE,show_row_names=TRUE,row_names_side="left",
  row_names_gp=grid::gpar(fontsize=if(nrow(M)>70)8.0 else 9.0,
                          col="black",fontface="bold"),
  width=grid::unit(ncol(M)+2,"cm"),
  height=grid::unit(0.50*nrow(M),"cm"),
  column_split=factor(colnames(M),levels=colnames(M)),column_title=NULL,
  layer_fun=function(j,i,x,y,w,h,fill) {
     grid::grid.text(sprintf("%.3f",M[cbind(i,j)]),x=x,y=y,
       gp=grid::gpar(fontsize=if(nrow(M)>70)HEAT_LABEL_FONT_PT else 9.5,
                     col="black",fontface="bold"))
  }
)

# Called once per graphics device (PDF and PNG), so both contain bar labels.
draw_I_polished <- function() {
  ComplexHeatmap::draw(hm,heatmap_legend_side="right",
                       annotation_legend_side="right")
  # Annotation-viewport y is bottom-up; matrix rows are top-down.
  # Use each bar midpoint and CENTER justification so labels never hug the right edge.
  n <- nrow(M)
  # x is the MIDPOINT of each bar. Keep the 0-1 AUC axis unchanged.
  xx <- as.numeric(avg_AUC) / 2
  yy <- (n - seq_len(n) + 0.5) / n
  # The values are DRAWN formatted, not used for ranking/selection.
  ComplexHeatmap::decorate_annotation("bar", {
    grid::grid.text(sprintf("%.3f",as.numeric(avg_AUC)),
      x=grid::unit(xx,"npc"), y=grid::unit(yy,"npc"),
      just="center",
      gp=grid::gpar(col="black",fontface="bold",
                    fontsize=if(n>80)BAR_LABEL_FONT_PT else 9.0))
  })
}
# More output height avoids cropped rows/labels; proportions match original source.
W <- 9.2
H <- max(11, (.50*nrow(M))/2.54 + 2.5)
grDevices::pdf(file.path(OUT,"Figure1I_model_AUC_OriginalStyle_CENTERED_BOLD_BLACK.pdf"),width=W,height=H,
               useDingbats=FALSE)
draw_I_polished()
grDevices::dev.off()
if (requireNamespace("ragg",quietly=TRUE)) {
 ragg::agg_png(file.path(OUT,"Figure1I_model_AUC_OriginalStyle_CENTERED_BOLD_BLACK.png"),
               width=W,height=H,units="in",res=400,background="white")
 draw_I_polished()
 grDevices::dev.off()
}
cat("Figure1I finished: ",nrow(M)," methods with valid external predictions\n",sep="")

} else {
  cat("Figure1I is preserved; regenerating Figure1J/K only.\n")
}

# -------- FIGURE J: original confusion matrix formatting --------------------
train_file <- file.path(RUN_DIR,"04_selected_model_train_OOF_predictions.tsv")
test_file <- file.path(RUN_DIR,"08_SELECTED_model_GSE18995_predictions.tsv")
if (!all(file.exists(c(train_file,test_file)))) {
 warning(" CV I J/K ")
} else {
 if (!requireNamespace("ggplot2",quietly=TRUE)) stop(" ggplot2")
 if (!requireNamespace("pROC",quietly=TRUE)) stop(" pROC")
 if (!requireNamespace("ragg",quietly=TRUE)) warning(" ragg, 600 dpi PNG")
 cv <- utils::read.delim(train_file,sep="\t",check.names=FALSE)
 ext <- utils::read.delim(test_file,sep="\t",check.names=FALSE)
 chosen <- r$Method[which(r$Status=="OK" & is.finite(r$CV_AUC))[1L]]
 if (!length(chosen)) stop("Check source data or configuration. ")
 COL <- list(blue="#4D9BC3",gold="#F2B45D",navy="#264653",
             teal="#2A9D8F",text="#000000",grid="#EAF0F4",
             reference="#A0AAB5",border="#D6DEE8")
 theme_original <- function(size=11.3) {
  ggplot2::theme_classic(base_size=size)+ggplot2::theme(
    text=ggplot2::element_text(face="bold",color="black"),
    plot.title=ggplot2::element_text(face="bold",hjust=.5,color=COL$text,
                                    size=size+1.9,margin=ggplot2::margin(b=12)),
    plot.tag=ggplot2::element_text(face="bold",size=size+5,color="black"),
    axis.text=ggplot2::element_text(color="black",face="bold",size=size-0.4),
    axis.title=ggplot2::element_text(color="black",face="bold",
                                    size=size+0.3),
    axis.line=ggplot2::element_line(color="#415566",linewidth=0.92),
    axis.ticks=ggplot2::element_line(color="#415566",linewidth=0.78),
    strip.background=ggplot2::element_rect(fill="#F8FAFC",
                 color="#C6D3DE",linewidth=0.85),
    strip.text=ggplot2::element_text(face="bold",size=10.3,color="black",
                                    margin=ggplot2::margin(t=6,b=6)),
    panel.border=ggplot2::element_rect(fill=NA,color="#AABBC9",linewidth=0.9),
    panel.grid=ggplot2::element_blank(),
    plot.background=ggplot2::element_rect(fill="white",color=NA),
    plot.margin=ggplot2::margin(12,14,12,14))
 }
 make_metrics <- function(d) {
  stopifnot(all(d$Type %in% c(0,1)),all(d$Predicted %in% c(0,1)))
  TP<-sum(d$Type==1 & d$Predicted==1)
  TN<-sum(d$Type==0 & d$Predicted==0)
  FN<-sum(d$Type==1 & d$Predicted==0)
  FP<-sum(d$Type==0 & d$Predicted==1)
  c(Accuracy=(TP+TN)/nrow(d),Sensitivity=TP/(TP+FN),Specificity=TN/(TN+FP))
 }
 # Keep identical columns before rbind: training contains Fold, external does not.
 # This fixes: Error in rbind(...): numbers of columns of arguments do not match.
 required_pred_cols <- c("Sample", "Donor", "Type", "Score_R2H", "Predicted")
 if (!all(required_pred_cols %in% names(cv)))
   stop("Check source data or configuration. ", paste(setdiff(required_pred_cols, names(cv)), collapse=", "))
 if (!all(required_pred_cols %in% names(ext)))
   stop("Check source data or configuration. ", paste(setdiff(required_pred_cols, names(ext)), collapse=", "))
 cv$Cohort <- "Train"
 ext$Cohort <- "GSE18995"
 cohort_order <- c("Train","GSE18995")
 plot_cols <- c("Cohort", required_pred_cols)
 combined <- rbind(cv[,plot_cols,drop=FALSE], ext[,plot_cols,drop=FALSE])
 combined$Actual <- factor(ifelse(combined$Type==1,"R2H","CIT"),levels=c("CIT","R2H"))
 combined$PredictedLabel <- factor(ifelse(combined$Predicted==1,"R2H","CIT"),levels=c("CIT","R2H"))
 combined$Cohort <- factor(combined$Cohort,levels=cohort_order)
 metric <- rbind(make_metrics(cv),make_metrics(ext))
 rownames(metric) <- cohort_order
 metric_table <- data.frame(Cohort=cohort_order,metric,check.names=FALSE)
 utils::write.table(metric_table,file.path(OUT,"Figure1J_classification_metrics.tsv"),
                    sep="\t",quote=FALSE,row.names=FALSE)
 labels <- setNames(sprintf("%s\n%s\nAccuracy=%.3f | Sensitivity=%.3f | Specificity=%.3f",
                            cohort_order,chosen,metric[,"Accuracy"],
                            metric[,"Sensitivity"],metric[,"Specificity"]),cohort_order)
 conf_df <- as.data.frame(table(Cohort=combined$Cohort,
                                 Actual=combined$Actual,Predicted=combined$PredictedLabel))
 names(conf_df)[4] <- "Count"
 conf_df$Correctness <- ifelse(as.character(conf_df$Actual)==as.character(conf_df$Predicted),
                               "Correct","Incorrect")
 # Same Y order as the original source: CIT appears above R2H.
 conf_df$PredictedPlot <- factor(as.character(conf_df$Predicted),levels=c("R2H","CIT"))
 pJ <- ggplot2::ggplot(conf_df,ggplot2::aes(x=Actual,y=PredictedPlot,fill=Correctness))+
   ggplot2::geom_tile(color="white",linewidth=CONFUSION_TILE_LINEWIDTH)+
   ggplot2::geom_text(ggplot2::aes(label=Count,color=Correctness),
                      size=5.5,fontface="bold")+
   ggplot2::scale_color_manual(values=c(Correct="black",Incorrect="black"),
                               guide="none")+
   ggplot2::facet_grid(Cohort~.,labeller=ggplot2::labeller(Cohort=labels))+
   ggplot2::scale_fill_manual(values=c(Correct=COL$blue,Incorrect=COL$gold))+
   ggplot2::coord_equal()+
   ggplot2::labs(title=paste0("Confusion matrices of the selected ",chosen," model"),
                 x="Actual class",y="Predicted class",fill=NULL,tag="J")+
   theme_original(11)+ggplot2::theme(legend.position="none",axis.ticks=ggplot2::element_blank())
 ggplot2::ggsave(file.path(OUT,"Figure1J_confusion_OriginalStyle_CENTERED_BOLD_BLACK.pdf"),pJ,
                 width=7.6,height=10.5,device=grDevices::cairo_pdf,bg="white")
 if (requireNamespace("ragg",quietly=TRUE))
 ggplot2::ggsave(file.path(OUT,"Figure1J_confusion_OriginalStyle_CENTERED_BOLD_BLACK.png"),pJ,
                 width=7.6,height=10.5,dpi=600,device=ragg::agg_png,bg="white")

 # -------- FIGURE K: preserve the original ROC presentation ----------------
 roc_df <- do.call(rbind,lapply(seq_along(cohort_order),function(i) {
    dat <- if (i==1L) cv else ext
    obj <- pROC::roc(dat$Type,dat$Score_R2H,levels=c(0,1),direction="<",quiet=TRUE)
    data.frame(FPR=1-obj$specificities,TPR=obj$sensitivities,Cohort=cohort_order[i])
 }))
 roc_df$Cohort <- factor(roc_df$Cohort,levels=cohort_order)
 roc_stat <- data.frame(Cohort=factor(cohort_order,levels=cohort_order),
   Label=c(sprintf("AUC = %.3f",as.numeric(pROC::auc(pROC::roc(cv$Type,cv$Score_R2H,levels=c(0,1),direction="<",quiet=TRUE)))),
           sprintf("AUC = %.3f",as.numeric(pROC::auc(pROC::roc(ext$Type,ext$Score_R2H,levels=c(0,1),direction="<",quiet=TRUE))))),
   X=.42,Y=.12)
 pK <- ggplot2::ggplot(roc_df,ggplot2::aes(FPR,TPR,color=Cohort))+
   ggplot2::geom_path(linewidth=ROC_LINEWIDTH,lineend="round",linejoin="round")+
   ggplot2::geom_abline(slope=1,intercept=0,linetype=2,
                         linewidth=0.90,color=COL$reference)+
   ggplot2::geom_label(data=roc_stat,ggplot2::aes(X,Y,label=Label),inherit.aes=FALSE,
                        color="black",fill="white",label.size=0.30,
                        size=3.3,fontface="bold",
                        label.padding=grid::unit(0.18,"lines"))+
   ggplot2::facet_grid(Cohort~.,labeller=ggplot2::labeller(Cohort=function(x)paste0(chosen,"\n",x)))+
   ggplot2::scale_color_manual(values=c(Train=COL$navy,GSE18995=COL$teal))+
   ggplot2::scale_x_continuous(breaks=seq(0,1,.2),limits=c(0,1))+
   ggplot2::scale_y_continuous(breaks=seq(0,1,.2),limits=c(0,1))+
   ggplot2::coord_equal()+ggplot2::labs(title="Receiver operating characteristic curves",
       x="1 - Specificity",y="Sensitivity",color=NULL,tag="K")+
   theme_original(11)+ggplot2::theme(legend.position="none")
 ggplot2::ggsave(file.path(OUT,"Figure1K_ROC_OriginalStyle_CENTERED_BOLD_BLACK.pdf"),pK,
                 width=7.6,height=10.5,device=grDevices::cairo_pdf,bg="white")
 if (requireNamespace("ragg",quietly=TRUE))
 ggplot2::ggsave(file.path(OUT,"Figure1K_ROC_OriginalStyle_CENTERED_BOLD_BLACK.png"),pK,
                 width=7.6,height=10.5,dpi=600,device=ragg::agg_png,bg="white")
 cat("Figure1J/K finished for selected model: ",chosen,"\n",sep="")
}
cat("Final figures in: ",normalizePath(OUT),"\n",sep="")



}

# =============================================================================
# STAGE figure2_shap_train
# =============================================================================
PIPELINE$figure2_shap_train <- function() {
# English consolidated source; frozen-model inference, no result replacement.
# ============================================================================
# RF+SVM | FROZEN-MODEL INTERVENTIONAL EXACT SHAP (GSE145989, 53 paired donors)
# 2026-10-21 | Compatible with ML_PAIR_FROM_SCRATCH_results output.
#
# RF was used for supervised FEATURE SELECTION. SVM is the fitted PREDICTOR.
# SHAP explains the FINAL SVM (sigmoid-transformed, class-1 decision margin).
# IMPORTANT: This SVM score is NOT a calibrated probability or clinical risk.
# There is NO retraining and NO external-based model/feature selection here.
# ============================================================================

options(stringsAsFactors = FALSE)
set.seed(123456)

# ------------------------------- CONFIG ------------------------------------
# Leave blank to select the EXACT input versions matching the frozen model.
# If you want to override, use TRAIN_FILE <- "data.train(3).txt", etc.
TRAIN_FILE <- repo_input("ml_train")
TEST_FILE <- repo_input("ml_test")
TRAIN_PAIR_FILE <- file.path(REPO_ROOT,"metadata/GSE145989_pairs.csv")
TEST_PAIR_FILE <- file.path(REPO_ROOT,"metadata/GSE18995_pairs.csv")
HELPER_FILE <- SCRIPT_FILE
RESULTS_ROOT <- repo_output("Figure1/ML")
RUN_DIR <- repo_ml_run()# Empty: auto-find one run. Or paste your run_* directory here.
# The GSE18995 SHAP folder is used ONLY to reuse identical BACKGROUND pairs,
# never to select the model or explain external samples.
MAX_BACKGROUND_PAIRS <- 20L  # Complete CIT+R2H training donor pairs.
EXACT_BATCH_MASKS <- 32L
BOOTSTRAP_B <- 1000L
PERMUTATION_B <- 100L
SEED <- 123456L
COL_CIT <- "#468CCB"
COL_R2H <- "#CC4F4E"
COL_OTHER <- "#91B9D2"
COL_DARK <- "#264653"
# ----------------------------------------------------------------------------

req <- c(TRAIN_PAIR_FILE,TEST_PAIR_FILE,HELPER_FILE)
missing <- req[!file.exists(req)]
if (length(missing)) stop("Missing project files: ",paste(missing,collapse=", "),
                          ". Set working directory to your original ML project.")
load_ml_helpers()
for (pkg in c("ggplot2","pROC"))
  if (!requireNamespace(pkg,quietly=TRUE)) stop("Install required package: ",pkg)

if (!nzchar(RUN_DIR)) {
  paths <- list.files(RESULTS_ROOT,pattern="^05_selected_frozen_model\\.rds$",
                      full.names=TRUE,recursive=TRUE)
  paths <- paths[basename(dirname(paths)) != "checkpoints_external"]
  if (length(paths)!=1L) stop("Detected ",length(paths)," frozen models. ",
     "Set RUN_DIR explicitly to the folder containing 05_selected_frozen_model.rds")
  RUN_DIR <- dirname(paths)
}
MODEL_FILE <- file.path(RUN_DIR,"05_selected_frozen_model.rds")
TRAIN_PRED_FILE <- file.path(RUN_DIR,"04_selected_model_train_OOF_predictions.tsv")
EXT_PRED_FILE <- file.path(RUN_DIR,"08_SELECTED_model_GSE18995_predictions.tsv")
for (f in c(MODEL_FILE,TRAIN_PRED_FILE,EXT_PRED_FILE))
  if (!file.exists(f)) stop("Saved model/prediction file missing: ",f)
OUT <- file.path(RUN_DIR,"RF_SVM_FROZEN_EXACT_SHAP_TRAIN_53PAIRS")
EXT_SHAP_DIR <- file.path(RUN_DIR,"RF_SVM_FROZEN_EXACT_SHAP_9GENES")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

read_expr <- function(path) {
  d <- utils::read.delim(path,header=TRUE,sep="\t",row.names=1,check.names=FALSE,
                         stringsAsFactors=FALSE,quote="",comment.char="")
  if (!"Type" %in% names(d) || anyDuplicated(rownames(d)))
    stop("Missing Type or duplicated rownames in ",path)
  ty <- tolower(trimws(as.character(d$Type)))
  y <- ifelse(ty %in% c("0","cit","control"),0L,
       ifelse(ty %in% c("1","r2h","treat","hr"),1L,NA_integer_))
  if (anyNA(y)) stop("Unrecognized class values")
  x <- as.data.frame(d[,setdiff(names(d),"Type"),drop=FALSE],check.names=FALSE)
  x[] <- lapply(x,function(t) as.numeric(as.character(t)))
  if (anyNA(x) || any(!is.finite(as.matrix(x)))) stop("Invalid numeric matrix: ",path)
  list(X=x,y=y,id=rownames(d))
}
# Auto-detect expression files: 21-gene and 22-gene variants can contain
# the same GSM IDs but different expression VALUES. Never silently mix them.
model <- readRDS(MODEL_FILE)
if (!nzchar(TRAIN_FILE)) {
  candidates <- c("data.train.txt","data.train(3).txt")
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("Neither data.train.txt nor data.train(3).txt exists")
  scaler <- model$scaler
  train_error <- vapply(candidates,function(path) {
    tryCatch({
      d <- read_expr(path)
      if (length(d$y)!=106L || !all(scaler$kept %in% colnames(d$X))) return(Inf)
      xx <- as.matrix(d$X[,scaler$kept,drop=FALSE])
      m <- colMeans(xx); sdv <- apply(xx,2L,stats::sd)
      max(c(abs(m-as.numeric(scaler$center[scaler$kept])),
            abs(sdv-as.numeric(scaler$scale[scaler$kept]))))
    },error=function(e) Inf)
  },numeric(1))
  if (all(!is.finite(train_error))) stop("No training matrix matches the saved model's predictors")
  TRAIN_FILE <- candidates[which.min(train_error)]
  message("Using training matrix ",TRAIN_FILE," (scaler discrepancy=",
          signif(min(train_error),4),")")
  if (min(train_error)>1e-5)
    stop("TRAIN_FILE does not match saved model training-scaler statistics. ",
         "Check model run and data versions.")
}
if (!nzchar(TEST_FILE)) {
  candidates <- c("data.test(3).txt","data.test.txt")
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("Neither data.test(3).txt nor data.test.txt exists")
  reference <- utils::read.delim(EXT_PRED_FILE,check.names=FALSE,stringsAsFactors=FALSE)
  test_error <- vapply(candidates,function(path) {
    tryCatch({
      d <- read_expr(path)
      if (length(d$y)!=32L || !setequal(d$id,as.character(reference$Sample))) return(Inf)
      p <- predict_ml_pipeline(model,d$X)
      max(abs(p-as.numeric(reference$Score_R2H[match(d$id,reference$Sample)])))
    },error=function(e) Inf)
  },numeric(1))
  if (all(!is.finite(test_error))) stop("No 32-sample GSE18995 matrix matches saved external scores")
  TEST_FILE <- candidates[which.min(test_error)]
  message("Using external matrix ",TEST_FILE," (saved-score discrepancy=",
          signif(min(test_error),4),")")
  if (min(test_error)>1e-6)
    stop("TEST_FILE cannot reproduce frozen GSE18995 predictions. ",
         "Check input matrices/helper versions.")
}
for (path in c(TRAIN_FILE,TEST_FILE))
  if (!file.exists(path)) stop("Selected expression file missing: ",path)
train <- read_expr(TRAIN_FILE)
test <- read_expr(TEST_FILE)
train_donor <- read_verified_pairs(TRAIN_PAIR_FILE,train$id,train$y,53L)
test_donor <- read_verified_pairs(TEST_PAIR_FILE,test$id,test$y,16L)
stopifnot(length(train$y)==106L,length(test$y)==32L)
if (!identical(model$label,"RF+SVM") || !identical(model$learner$method,"SVM")) {
  stop("Frozen model is NOT RF+SVM. Found: ",model$label,
       " (learner=",model$learner$method,"). Do not mislabel another classifier.")
}
genes <- as.character(model$selected_genes)
if (length(genes)!=9L || anyDuplicated(genes)) {
  stop("RF+SVM currently has ",length(genes)," selected genes, not the nine ",
       "in the supplied ranking. Confirm this RUN_DIR belongs to the same run.")
}
if (!all(genes %in% names(train$X)) || !all(genes %in% names(test$X)))
  stop("Selected genes missing from input expression matrices")
message("FINAL classifier: RF+SVM | selected genes: ",paste(genes,collapse=", "))

train_oof <- utils::read.delim(TRAIN_PRED_FILE,check.names=FALSE)
ext_oof <- utils::read.delim(EXT_PRED_FILE,check.names=FALSE)
for (d in list(train_oof,ext_oof)) {
  if (!all(c("Sample","Donor","Type","Score_R2H") %in% colnames(d)))
    stop("Prediction tables require Sample, Donor, Type, Score_R2H")
}
# Match columns to original input order rather than relying on saved row order.
align_pred <- function(pred,input,donor) {
  ii <- match(input$id,as.character(pred$Sample))
  if (anyNA(ii) || anyDuplicated(pred$Sample) || nrow(pred)!=length(input$id))
    stop("Saved prediction table does not match matrix samples")
  pred <- pred[ii,,drop=FALSE]
  if (any(as.integer(pred$Type)!=input$y) ||
      any(as.character(pred$Donor)!=as.character(donor)))
    stop("Saved predictions disagree with verified donor/case labels")
  pred
}
train_oof <- align_pred(train_oof,train,train_donor)
ext_oof <- align_pred(ext_oof,test,test_donor)
p_external <- predict_ml_pipeline(model,test$X)
if (max(abs(p_external-as.numeric(ext_oof$Score_R2H)))>1e-6)
  stop("Frozen model predictions do not reproduce saved GSE18995 scores. ",
       "Check RUN_DIR/input data/helper versions!")

calc_auc <- function(y,p) as.numeric(pROC::auc(
  pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE)))
auc_oof <- calc_auc(train$y,as.numeric(train_oof$Score_R2H))
auc_ext <- calc_auc(test$y,p_external)
cat("Verified AUC: grouped OOF=",sprintf("%.6f",auc_oof),
    ", external=",sprintf("%.6f",auc_ext),"\n",sep="")

# Paired (clustered) bootstrap: resample DONORS, retaining both biopsies.
set.seed(SEED+11L)
ci_oof <- paired_auc_bootstrap(train$y,as.numeric(train_oof$Score_R2H),
                               train_donor,B=BOOTSTRAP_B,seed=SEED+11L)
ci_ext <- paired_auc_bootstrap(test$y,p_external,test_donor,
                               B=BOOTSTRAP_B,seed=SEED+12L)
auc_summary <- data.frame(
  Dataset=c("Training donor-grouped 5-fold OOF","GSE18995 external"),
  Donor_pairs=c(53,16),AUC=c(auc_oof,auc_ext),
  CI_low=c(ci_oof[1],ci_ext[1]),CI_high=c(ci_oof[2],ci_ext[2]))
write.csv(auc_summary,file.path(OUT,"RF_SVM_AUC_paired_donor_bootstrap.csv"),row.names=FALSE)
write.csv(data.frame(Gene=genes,Source="RF selector, final full-training refit"),
          file.path(OUT,"RF_SVM_final_9_features.csv"),row.names=FALSE)

# Explain the frozen FULL-training SVM on training samples (NOT OOF SHAP).
# Panel A ROC still uses saved donor-grouped TRAINING OOF predictions and
# GSE18995 external predictions. SHAP B-F explains full-training frozen model.
explain <- train
explain_donor <- train_donor
p_explain <- predict_ml_pipeline(model,explain$X)
if (length(p_explain) != 106L || any(!is.finite(p_explain)))
  stop("Could not compute 106 full-training model scores")
auc_train_apparent <- calc_auc(explain$y,p_explain)
cat("Training SHAP: 53 pairs, 106 samples; frozen-model fitted AUC=",
    sprintf("%.3f",auc_train_apparent),"; grouped OOF AUC=",
    sprintf("%.3f",auc_oof),"\n",sep="")

# Keep exactly the SAME training-donor background as completed external SHAP.
# If an external background CSV exists, follow its sample IDs IN SAVED ORDER.
external_bg_csv <- file.path(EXT_SHAP_DIR,"SHAP_background_training_donor_pairs.csv")
if (file.exists(external_bg_csv)) {
  bg_tab <- utils::read.csv(external_bg_csv,stringsAsFactors=FALSE,check.names=FALSE)
  if (!all(c("Sample","Donor","Group") %in% names(bg_tab)) ||
      anyDuplicated(bg_tab$Sample) || nrow(bg_tab) %% 2L != 0L)
    stop("External SHAP background CSV is invalid")
  bg_idx <- match(as.character(bg_tab$Sample),train$id)
  if (anyNA(bg_idx) || anyDuplicated(bg_idx) ||
      any(as.character(bg_tab$Donor)!=as.character(train_donor[bg_idx])) ||
      any(as.character(bg_tab$Group)!=ifelse(train$y[bg_idx]==0,"CIT","R2H")))
    stop("Original external-SHAP background cannot be matched to this training matrix")
  chosen_bg_donors <- sort(unique(as.character(train_donor[bg_idx])))
  cat("Reusing EXTERNAL SHAP's original ",length(chosen_bg_donors),
      " training-background donor pairs\n",sep="")
} else {
  set.seed(SEED+21L)
  chosen_bg_donors <- sample(sort(unique(as.character(train_donor))),
                             min(MAX_BACKGROUND_PAIRS,length(unique(train_donor))))
  bg_idx <- which(as.character(train_donor) %in% chosen_bg_donors)
  warning("External SHAP background CSV not found; using deterministic same-seed selection. ",
          "For strict comparability, provide the earlier background CSV.")
}
if (length(bg_idx)!=2L*length(chosen_bg_donors) ||
    !all(table(train_donor[bg_idx],factor(train$y[bg_idx],levels=0:1))==1L))
  stop("Background does not contain complete CIT/R2H donor pairs")
# The frozen pipeline scales ALL kept training predictors BEFORE restricting
# to the RF-selected 9 genes. Therefore coalition data must contain the full
# scaler$kept set, even though only selected_genes enter the SHAP coalitions.
required_model_cols <- as.character(model$scaler$kept)
if (!length(required_model_cols) || anyNA(required_model_cols) ||
    anyDuplicated(required_model_cols) ||
    !all(required_model_cols %in% colnames(train$X)) ||
    !all(required_model_cols %in% colnames(test$X)))
  stop("Frozen model scaler$kept is incompatible with expression input columns")
if (!all(genes %in% required_model_cols))
  stop("RF-selected genes are not a subset of frozen scaler$kept")
background <- train$X[bg_idx,required_model_cols,drop=FALSE]
stopifnot(nrow(background)==2L*length(chosen_bg_donors))
write.csv(data.frame(Sample=train$id[bg_idx],Donor=train_donor[bg_idx],
                     Group=ifelse(train$y[bg_idx]==0,"CIT","R2H")),
          file.path(OUT,"SHAP_background_training_donor_pairs.csv"),row.names=FALSE)

# Ensure training-vs-external SHAP uses the same frozen model, background, and score scale.
# The interventional SHAP baseline equals the mean prediction on the background.
background_score <- mean(predict_ml_pipeline(model,background))
previous_ext_rds <- file.path(EXT_SHAP_DIR,"RF_SVM_exact_SHAP_result.rds")
if (file.exists(previous_ext_rds)) {
  old_ext <- readRDS(previous_ext_rds)
  if (!is.null(old_ext$baseline) &&
      (any(!is.finite(old_ext$baseline)) ||
       max(abs(old_ext$baseline-background_score))>1e-6)) {
    stop("Training SHAP background or expression matrix differs from external SHAP. ",
         "Check TRAIN_FILE and model versions before comparison.")
  }
  message("PASS: SHAP reference baseline matches previous external SHAP")
}

# EXACT interventional Shapley (2^9 coalitions).
# Value function: mean SVM 0-1 transformed margin after interventions on
# a training-donor-pair background. This is NOT a calibrated probability.
# Only selected gene columns are intervened; full scaler input is preserved.
# Inference of SHAP's values is descriptive; donors are the resampling units.
exact_one <- function(x,model,background,batch_masks=32L) {
  # Explain only the nine SELECTED genes; other scaler columns are retained
  # as required prediction inputs but never enumerated as SHAP features.
  genes <- as.character(model$selected_genes)
  if (!all(genes %in% colnames(background)))
    stop("Exact SHAP background is missing RF-selected genes")
  p <- length(genes)
  if (p>12L) stop("Exact SHAP intentionally limited to <=12 genes")
  x <- as.numeric(x[genes]); names(x) <- genes
  nmask <- 2L^p
  masks <- 0L:(nmask-1L)
  bits <- sapply(seq_len(p),function(j)
    bitwAnd(masks,bitwShiftL(1L,j-1L))!=0L)
  if (p==1L) bits <- matrix(bits,ncol=1L)
  K <- rowSums(bits)
  v <- numeric(nmask)
  B <- nrow(background)
  for (lo in seq.int(1L,nmask,by=as.integer(batch_masks))) {
    rows <- lo:min(nmask,lo+as.integer(batch_masks)-1L)
    data_b <- background[rep(seq_len(B),times=length(rows)),,drop=FALSE]
    for (h in seq_along(rows)) {
      use_cols <- which(bits[rows[h],])
      if (length(use_cols)) {
        rr <- ((h-1L)*B+1L):(h*B)
        # Match by gene NAME: selected genes need not be the first/adjacent
        # columns in the full scaler input matrix.
        for (j in use_cols) data_b[rr,genes[j]] <- x[j]
      }
    }
    pred <- predict_ml_pipeline(model,data_b)
    v[rows] <- colMeans(matrix(pred,nrow=B,ncol=length(rows)))
  }
  phi <- setNames(numeric(p),genes)
  for (j in seq_len(p)) {
    unmasked <- which(!bits[,j])
    w <- 1 / (p*choose(p-1,K[unmasked]))
    with_j <- masks[unmasked] + bitwShiftL(1L,j-1L) + 1L
    phi[j] <- sum(w*(v[with_j]-v[unmasked]))
  }
  list(shap=phi,baseline=v[1L],prediction=v[nmask],
       reconstruction=v[1L]+sum(phi))
}

# Resume SHAP by sample after interruption; refuse mismatched inputs/background.
signature <- c("TRAIN-53pairs-exact-shap-full-scaler-input-v3",
               unname(tools::md5sum(c(MODEL_FILE,TRAIN_FILE,HELPER_FILE))),
               paste(sort(bg_idx),collapse=","),paste(genes,collapse=";"),
               paste(required_model_cols,collapse=";"))
cache_file <- file.path(OUT,"CHECKPOINT_exact_SHAP.rds")
S <- matrix(NA_real_,nrow=length(explain$id),ncol=length(genes),
            dimnames=list(explain$id,genes))
baseline <- prediction <- rep(NA_real_,length(explain$id))
if (file.exists(cache_file)) {
  old <- tryCatch(readRDS(cache_file),error=function(e) NULL)
  if (!is.null(old) && identical(old$signature,signature) &&
      identical(dimnames(old$S),dimnames(S))) {
    S <- old$S; baseline <- old$baseline; prediction <- old$prediction
    cat("Resuming exact SHAP from checkpoint\n")
  } else stop("Existing SHAP checkpoint belongs to different data/model/background. ",
              "Rename/remove CHECKPOINT_exact_SHAP.rds before recomputing.")
}
for (i in seq_along(explain$id)) {
  if (all(is.finite(S[i,])) && is.finite(baseline[i])) next
  ans <- exact_one(explain$X[i,genes,drop=TRUE],model,background,
                   batch_masks=EXACT_BATCH_MASKS)
  if (abs(ans$prediction-p_explain[i])>1e-6 ||
      abs(ans$reconstruction-p_explain[i])>1e-6)
    stop("Exact SHAP additivity check failed at sample: ",explain$id[i])
  S[i,] <- ans$shap; baseline[i] <- ans$baseline; prediction[i] <- ans$prediction
  saveRDS(list(signature=signature,S=S,baseline=baseline,prediction=prediction),
          cache_file)
  cat(sprintf("Exact SHAP %2d/%2d: %s\n",i,length(explain$id),explain$id[i]))
}
if (max(abs(baseline+rowSums(S)-p_explain))>1e-6)
  stop("Final exact SHAP additivity mismatch")
saveRDS(list(S=S,baseline=baseline,prediction=prediction,
             model_label=model$label,selected_genes=genes,
             prediction_type="uncalibrated SVM margin transformed with sigmoid"),
        file.path(OUT,"RF_SVM_exact_SHAP_result.rds"))
write.csv(data.frame(Sample=explain$id,Donor=as.character(explain_donor),
                     Group=ifelse(explain$y==0,"CIT","R2H"),
                     Baseline=baseline,Predicted_Score=prediction,
                     Reconstruction=baseline+rowSums(S)),
          file.path(OUT,"RF_SVM_exact_SHAP_additivity.csv"),row.names=FALSE)

# Global importance: bootstrap COMPLETE TRAINING DONOR PAIRS (n=53).
set.seed(SEED+91L)
u <- sort(unique(as.character(explain_donor)))
imp_mat <- replicate(BOOTSTRAP_B,{
  donor_pick <- sample(u,length(u),replace=TRUE)
  ix <- unlist(lapply(donor_pick,function(d) which(explain_donor==d)),use.names=FALSE)
  colMeans(abs(S[ix,,drop=FALSE]))
})
imp <- data.frame(
  Gene=genes,MeanAbsSHAP=colMeans(abs(S)),
  CI_Lower=apply(imp_mat,1,stats::quantile,probs=.025,na.rm=TRUE),
  CI_Upper=apply(imp_mat,1,stats::quantile,probs=.975,na.rm=TRUE))
imp <- imp[order(-imp$MeanAbsSHAP),,drop=FALSE]
imp$Rank <- seq_len(nrow(imp))
write.csv(imp,file.path(OUT,"RF_SVM_SHAP_global_importance_53pair_CI.csv"),row.names=FALSE)

# Optional cross-cohort SHAP stability table (no refitting and no new inference).
ext_imp_file <- file.path(EXT_SHAP_DIR,"RF_SVM_SHAP_global_importance_16pair_CI.csv")
if (file.exists(ext_imp_file)) {
  ext_imp <- utils::read.csv(ext_imp_file,stringsAsFactors=FALSE,check.names=FALSE)
  cols <- c("Gene","MeanAbsSHAP","Rank","CI_Lower","CI_Upper")
  if (all(cols %in% colnames(ext_imp)) && setequal(ext_imp$Gene,imp$Gene)) {
    a <- imp[,cols,drop=FALSE]
    b <- ext_imp[,cols,drop=FALSE]
    names(a)[-1L] <- paste0("Train_",names(a)[-1L])
    names(b)[-1L] <- paste0("GSE18995_",names(b)[-1L])
    cmp <- merge(a,b,by="Gene",sort=FALSE)
    cmp <- cmp[order(cmp$Train_Rank),,drop=FALSE]
    write.csv(cmp,file.path(OUT,"RF_SVM_train_vs_external_SHAP_comparison.csv"),row.names=FALSE)
  }
}

# Each training sample appears once per selected gene in beeswarm/dependence data.
long <- do.call(rbind,lapply(genes,function(g) {
  raw <- as.numeric(explain$X[,g])
  ra <- range(raw)
  scaled <- if (diff(ra)==0) rep(.5,length(raw)) else (raw-ra[1])/diff(ra)
  data.frame(Sample=explain$id,Donor=as.character(explain_donor),
             Group=ifelse(explain$y==0,"CIT","R2H"),
             Gene=g,Expression=raw,ExpressionScaled=scaled,
             SHAP=S[,g],stringsAsFactors=FALSE)
}))
write.csv(long,file.path(OUT,"RF_SVM_SHAP_long_format.csv"),row.names=FALSE)

# Donor-paired gene expression analysis: 9 tests jointly BH corrected.
paired_vectors <- function(v,donor,y) {
  di <- sort(unique(as.character(donor)))
  cit <- r2h <- numeric(length(di))
  for (i in seq_along(di)) {
    a <- which(donor==di[i] & y==0L)
    b <- which(donor==di[i] & y==1L)
    if (length(a)!=1L || length(b)!=1L) stop("Unmatched donor: ",di[i])
    cit[i] <- v[a]; r2h[i] <- v[b]
  }
  list(donor=di,CIT=cit,R2H=r2h,Delta=r2h-cit)
}
p_raw <- sapply(genes,function(g) {
  q <- paired_vectors(as.numeric(explain$X[,g]),explain_donor,explain$y)
  if (all(q$Delta==0)) return(1)
  stats::wilcox.test(q$R2H,q$CIT,paired=TRUE,exact=FALSE)$p.value
})
p_adj <- stats::p.adjust(p_raw,method="BH")
paired_gene <- data.frame(Gene=genes,P_value=as.numeric(p_raw),
                          FDR_BH=as.numeric(p_adj))
write.csv(paired_gene,file.path(OUT,"GSE145989_9genes_paired_Wilcoxon_BH.csv"),row.names=FALSE)

# 9-gene paired coordinated program score from TRAINING mean/SD only.
explain_z <- transform_by_training(explain$X,model$scaler)[,genes,drop=FALSE]
score9 <- rowMeans(as.matrix(explain_z))
pair_score9 <- paired_vectors(score9,explain_donor,explain$y)
p_pair_score9 <- stats::wilcox.test(pair_score9$R2H,pair_score9$CIT,
                                    paired=TRUE,exact=FALSE)$p.value
write.csv(data.frame(Sample=explain$id,Donor=explain_donor,
                     Group=ifelse(explain$y==0,"CIT","R2H"),Score9=score9),
          file.path(OUT,"GSE145989_9gene_score.csv"),row.names=FALSE)

# Additional: paired DELTA expression correlation across 53 independent donors.
delta_expr <- sapply(genes,function(g)
  paired_vectors(as.numeric(explain$X[,g]),explain_donor,explain$y)$Delta)
colnames(delta_expr) <- genes
rho <- stats::cor(delta_expr,method="spearman",use="pairwise.complete.obs")
write.csv(rho,file.path(OUT,"GSE145989_9gene_paired_delta_Spearman.csv"))

# Training paired permutation importance (APPARENT, not out-of-fold):
# permute each gene across donor PAIRS
# (copy both CIT and R2H of the source donor together, by condition).
# This is an exploratory supplementary importance metric, not SHAP and not
# used for model selection.
set.seed(SEED+121L)
perm_drop <- matrix(NA_real_,nrow=PERMUTATION_B,ncol=length(genes),
                    dimnames=list(NULL,genes))
for (j in seq_along(genes)) {
  gene <- genes[j]
  for (b in seq_len(PERMUTATION_B)) {
    source_donors <- sample(u,length(u),replace=FALSE)
    xp <- explain$X
    for (t in seq_along(u)) {
      for (grp in 0:1) {
        ii <- which(explain_donor==u[t] & explain$y==grp)
        jj <- which(explain_donor==source_donors[t] & explain$y==grp)
        xp[ii,gene] <- explain$X[jj,gene]
      }
    }
    perm_drop[b,j] <- auc_train_apparent-calc_auc(explain$y,predict_ml_pipeline(model,xp))
  }
}
perm_imp <- data.frame(Gene=genes,MeanAUCdrop=colMeans(perm_drop),
                       Lower=apply(perm_drop,2,stats::quantile,.025),
                       Upper=apply(perm_drop,2,stats::quantile,.975))
perm_imp <- perm_imp[order(-perm_imp$MeanAUCdrop),]
write.csv(perm_imp,file.path(OUT,"RF_SVM_53pair_permutation_AUC_drop.csv"),row.names=FALSE)

# --------------------- Publication-style Figure 2 panels --------------------
library(ggplot2)
# Publication typography: all TEXT is bold black; original figure COLORS retained.
base_theme <- function(size=11) {
  theme_classic(base_size=size,base_family="sans") +
    theme(
      text=element_text(family="sans",face="bold",color="black"),
      plot.title=element_text(face="bold",color="black",hjust=.5,
                              size=size+1.8,margin=margin(b=9)),
      plot.subtitle=element_text(face="bold",color="black",hjust=.5,
                                 size=size-.4,margin=margin(b=8)),
      axis.title=element_text(face="bold",color="black",size=size+.4),
      axis.text=element_text(face="bold",color="black",size=size-.5),
      legend.title=element_text(face="bold",color="black",size=size-.3),
      legend.text=element_text(face="bold",color="black",size=size-.8),
      strip.text=element_text(face="bold",color="black"),
      plot.tag=element_text(face="bold",color="black"),
      axis.line=element_line(color="black",linewidth=.78),
      axis.ticks=element_line(color="black",linewidth=.65),
      panel.border=element_rect(color="black",fill=NA,linewidth=.72),
      legend.background=element_rect(fill="white",color=NA),
      legend.key=element_rect(fill="white",color=NA),
      plot.margin=margin(10,12,10,10))
}
save_panel <- function(p,name,width=7.3,height=5.5) {
  ggplot2::ggsave(file.path(OUT,paste0(name,".pdf")),p,
                  width=width,height=height,bg="white",useDingbats=FALSE)
  if (requireNamespace("ragg",quietly=TRUE)) {
    ggplot2::ggsave(file.path(OUT,paste0(name,".png")),p,
                    width=width,height=height,dpi=600,
                    device=ragg::agg_png,bg="white")
  } else {
    ggplot2::ggsave(file.path(OUT,paste0(name,".png")),p,
                    width=width,height=height,dpi=300,bg="white")
  }
}

roc_d <- function(y,p,name) {
  o <- pROC::roc(y,p,direction="<",levels=c(0,1),quiet=TRUE)
  data.frame(FPR=1-o$specificities,TPR=o$sensitivities,Curve=name)
}
roc_df <- rbind(roc_d(train$y,as.numeric(train_oof$Score_R2H),"Training grouped OOF"),
                roc_d(test$y,p_external,"GSE18995 external"))
leg <- c(sprintf("Training grouped OOF: AUC=%.3f",auc_oof),
         sprintf("GSE18995 external: AUC=%.3f",auc_ext))
names(leg) <- c("Training grouped OOF","GSE18995 external")
pA <- ggplot(roc_df,aes(FPR,TPR,color=Curve))+
  geom_path(linewidth=1.55,lineend="round",linejoin="round")+
  geom_abline(linetype=2,color="#999999",linewidth=.78)+
  scale_color_manual(values=c("Training grouped OOF"=COL_DARK,
                              "GSE18995 external"="#16A085"),labels=leg)+
  coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+
  labs(title="RF+SVM discriminatory performance",x="1 - Specificity",
       y="Sensitivity",color=NULL)+base_theme(10.5)+
  # Move the TWO legend rows from lower-right inward; retain line colors.
  theme(legend.position=c(.58,.23),legend.justification=c(.5,.5),
        legend.direction="vertical",
        legend.text=element_text(face="bold",color="black",size=9.0),
        legend.background=element_rect(fill=grDevices::adjustcolor("white",alpha.f=.92),color=NA))
save_panel(pA,"Panel_A_RF_SVM_ROC")

imp$Gene <- factor(imp$Gene,levels=rev(as.character(imp$Gene)))
pB <- ggplot(imp,aes(MeanAbsSHAP,Gene,fill=Gene=="SERPINE1"))+
  geom_col(width=.67)+
  geom_segment(aes(x=CI_Lower,xend=CI_Upper,y=Gene,yend=Gene),
               color=COL_DARK,linewidth=.65)+
  geom_point(shape=21,size=2,fill="white",color=COL_DARK)+
  scale_fill_manual(values=c(`TRUE`=COL_R2H,`FALSE`=COL_OTHER),guide="none")+
  scale_x_continuous(expand=expansion(mult=c(0,.22)))+
  geom_text(aes(x=CI_Upper,label=sprintf("%.3f",MeanAbsSHAP)),
            hjust=-.12,size=3.2,fontface="bold",color="black")+
  labs(title="Global exact SHAP importance",subtitle="GSE145989: 53-donor-pair bootstrap 95% CI",
       x="Mean |SHAP| on SVM decision score",y=NULL)+base_theme(10.5)
save_panel(pB,"Panel_B_RF_SVM_SHAP_global")

long$Gene <- factor(long$Gene,levels=rev(as.character(imp$Gene)))
if (requireNamespace("ggbeeswarm",quietly=TRUE)) {
  pt_layer <- ggbeeswarm::geom_quasirandom(width=.26,alpha=.82,size=1.9)
} else {
  set.seed(SEED)
  pt_layer <- geom_jitter(width=0,height=.18,alpha=.82,size=1.9)
}
pC <- ggplot(long,aes(SHAP,Gene,color=ExpressionScaled))+
  geom_vline(xintercept=0,linetype=2,color="#888888")+
  pt_layer+
  scale_color_gradientn(colors=c("#2166AC","#F7F7F7","#B2182B"),
                        limits=c(0,1),breaks=c(0,1),labels=c("Low","High"))+
  labs(title="Exact SHAP beeswarm",subtitle="Positive SHAP shifts score toward R2H",
       x="SVM score contribution (uncalibrated)",y=NULL,color="Expression")+
  base_theme(10.5)
save_panel(pC,"Panel_C_RF_SVM_SHAP_beeswarm")

ser <- long[as.character(long$Gene)=="SERPINE1",,drop=FALSE]
ser$Group <- factor(ser$Group,levels=c("CIT","R2H"))
pD <- ggplot(ser,aes(Expression,SHAP,color=Group))+
  geom_hline(yintercept=0,linetype=2,color="#999999")+
  geom_point(size=2.5,alpha=.85)+
  geom_smooth(aes(group=1),method="loess",formula=y~x,se=TRUE,
              inherit.aes=TRUE,color=COL_DARK,linewidth=.85)+
  scale_color_manual(labels=c(CIT="CIT",R2H="HR"),values=c(CIT=COL_CIT,R2H=COL_R2H))+
  labs(title="SERPINE1 SHAP dependence",subtitle="Exploratory model explanation",
       x="SERPINE1 expression",y="Exact SHAP (SVM decision score)",color=NULL)+
  base_theme(10.5)
save_panel(pD,"Panel_D_SERPINE1_dependence")

citser <- paired_vectors(as.numeric(explain$X[,"SERPINE1"]),explain_donor,explain$y)
paired_ser <- data.frame(Donor=rep(citser$donor,2L),
                         Group=factor(rep(c("CIT","R2H"),each=length(citser$donor)),
                                      levels=c("CIT","R2H")),
                         Expression=c(citser$CIT,citser$R2H))
ser_fdr <- paired_gene$FDR_BH[paired_gene$Gene=="SERPINE1"]
pE <- ggplot(paired_ser,aes(Group,Expression,group=Donor))+
  geom_line(color="#A9A9A9",linewidth=.5,alpha=.7)+
  geom_point(aes(color=Group),size=2.9,alpha=.93)+
  scale_color_manual(labels=c(CIT="CIT",R2H="HR"),values=c(CIT=COL_CIT,R2H=COL_R2H),guide="none")+
  labs(title="Training paired SERPINE1 expression",
       subtitle=sprintf("53 donor pairs | Wilcoxon BH-FDR=%.4g",ser_fdr),
       x=NULL,y="Normalized expression")+base_theme(10.5)
save_panel(pE,"Panel_E_GSE145989_SERPINE1_paired")

# Typical local case chosen by MEDIAN final-model score among 53 R2H TRAIN samples.
rr <- which(explain$y==1)
rr <- rr[order(abs(p_explain[rr]-stats::median(p_explain[rr])))]
k <- rr[1L]
phi <- S[k,]
phi <- phi[order(abs(phi),decreasing=TRUE)]
starts <- c(baseline[k],baseline[k]+head(cumsum(phi),-1L))
ends <- starts+phi
wf <- data.frame(Gene=factor(names(phi),levels=names(phi)),
                 i=seq_along(phi),Start=starts,End=ends,
                 Contribution=as.numeric(phi),
                 Direction=ifelse(phi>=0,"Toward HR","Toward CIT"))
pF <- ggplot(wf)+
  geom_rect(aes(xmin=i-.36,xmax=i+.36,
                ymin=pmin(Start,End),ymax=pmax(Start,End),fill=Direction))+
  geom_segment(data=utils::head(wf,-1L),
               aes(x=i+.36,xend=i+1-.36,y=End,yend=End),linewidth=.5,color=COL_DARK)+
  geom_hline(yintercept=baseline[k],linetype=2,color="#888888")+
  geom_hline(yintercept=p_explain[k],linetype=3,color=COL_DARK)+
  geom_text(aes(x=i,y=pmax(Start,End)+.02,label=sprintf("%+.3f",Contribution)),
            size=3,fontface="bold",color="black")+
  scale_x_continuous(breaks=wf$i,labels=as.character(wf$Gene))+
  scale_fill_manual(values=c("Toward HR"=COL_R2H,"Toward CIT"=COL_CIT),guide="none")+
  scale_y_continuous(limits=c(min(0,min(starts,ends)-.08),
                              max(1.05,max(starts,ends)+.12)))+
  labs(title="Local exact SHAP waterfall",
       subtitle=paste0("Training HR median-score sample: ",explain$id[k],
         " | baseline=",sprintf("%.3f",baseline[k]),
         ", score=",sprintf("%.3f",p_explain[k])),
       x=NULL,y="SVM decision score (uncalibrated)")+
  base_theme(10)+theme(axis.text.x=element_text(angle=50,hjust=1,
                                                     face="bold",color="black"))
save_panel(pF,"Panel_F_RF_SVM_local_waterfall",width=9,height=5.5)

# Supplement: 9-gene paired score.
sc_df <- data.frame(Donor=rep(pair_score9$donor,2),
                    Group=factor(rep(c("CIT","R2H"),each=length(pair_score9$donor)),
                                 levels=c("CIT","R2H")),
                    Score=c(pair_score9$CIT,pair_score9$R2H))
pS1 <- ggplot(sc_df,aes(Group,Score,group=Donor))+
  geom_line(color="grey70",alpha=.7)+
  geom_point(aes(color=Group),size=2.6)+
  scale_color_manual(labels=c(CIT="CIT",R2H="HR"),values=c(CIT=COL_CIT,R2H=COL_R2H),guide="none")+
  labs(title="RF-selected nine-gene coordinated score",
       subtitle=sprintf("GSE145989 donor-paired Wilcoxon P=%.4g",p_pair_score9),
       x=NULL,y="Mean training-referenced z-score")+base_theme(10.5)
save_panel(pS1,"Supplement_1_GSE145989_9gene_paired_score")

pS2 <- ggplot(perm_imp,aes(MeanAUCdrop,reorder(Gene,MeanAUCdrop),fill=Gene=="SERPINE1"))+
  geom_col(width=.7)+
  geom_segment(aes(x=Lower,xend=Upper,y=Gene,yend=Gene),linewidth=.7,color=COL_DARK)+
  scale_fill_manual(values=c(`TRUE`=COL_R2H,`FALSE`=COL_OTHER),guide="none")+
  labs(title="Training paired-donor permutation importance",
       subtitle="Exploratory: mean AUC reduction from donor-pair exchange",
       x="AUC decrease after paired permutation",y=NULL)+base_theme(10.5)
save_panel(pS2,"Supplement_2_GSE145989_paired_permutation_importance")

if (requireNamespace("patchwork",quietly=TRUE)) {
  fig <- (pA+pB+pC)/(pD+pE+pF)+
    patchwork::plot_layout(guides="keep",heights=c(1,1))
  save_panel(fig,"Figure2_RF_SVM_SHAP_A_to_F",width=19,height=12)
}
writeLines(c(
  "RF+SVM frozen exact interventional SHAP (9 final RF-selected genes)",
  "SVM output: sigmoid-transformed oriented decision margin; uncalibrated (NOT a risk probability)",
  "Main-model selection was based on GSE145989 donor-grouped OOF AUC; no external-based re-selection",
  sprintf("RF+SVM OOF AUC %.6f; GSE18995 external AUC %.6f",auc_oof,auc_ext),
  "GSE145989 SHAP values refer to the final full-training SVM, not to the fold-specific CV models",
  "SHAP interventional background: balanced CIT/R2H biopsies from complete training-donor pairs",
  "SHAP importance CI: bootstrap of 53 complete training donor pairs (conditional on frozen fit)",
  "D: SERPINE1 expression vs SHAP is an exploratory model dependence, not causal evidence",
  "CIT/R2H inference: training donor-paired Wilcoxon; nine individual gene tests BH adjusted",
  sprintf("Final-model fitted training AUC %.6f (SHAP explanation); do NOT report as OOF AUC",auc_train_apparent),
  "Potential limitation: if the upstream candidate gene panel was selected using all 53 training donors,",
  "the 5-fold machine-learning OOF AUC does not provide full pipeline-nested gene-discovery validation",
  "Training SHAP is an IN-SAMPLE final-model explanation; it is NOT grouped-OOF SHAP"),
  con=file.path(OUT,"READ_THIS_Methods_and_interpretation.txt"))
cat("\nSHAP finished; all output: ",normalizePath(OUT),"\n",sep="")

}

# =============================================================================
# STAGE figure2_shap_test
# =============================================================================
PIPELINE$figure2_shap_test <- function() {
# External frozen-model explanation; training-derived background only.
# English consolidated source; frozen-model inference, no result replacement.
# ============================================================================
# RF+SVM | FROZEN-MODEL INTERVENTIONAL EXACT SHAP (GSE145989, 53 paired donors)
# 2026-10-21 | Compatible with ML_PAIR_FROM_SCRATCH_results output.
#
# RF was used for supervised FEATURE SELECTION. SVM is the fitted PREDICTOR.
# SHAP explains the FINAL SVM (sigmoid-transformed, class-1 decision margin).
# IMPORTANT: This SVM score is NOT a calibrated probability or clinical risk.
# There is NO retraining and NO external-based model/feature selection here.
# ============================================================================

options(stringsAsFactors = FALSE)
set.seed(123456)

# ------------------------------- CONFIG ------------------------------------
# Leave blank to select the EXACT input versions matching the frozen model.
# If you want to override, use TRAIN_FILE <- "data.train(3).txt", etc.
TRAIN_FILE <- repo_input("ml_train")
TEST_FILE <- repo_input("ml_test")
TRAIN_PAIR_FILE <- file.path(REPO_ROOT,"metadata/GSE145989_pairs.csv")
TEST_PAIR_FILE <- file.path(REPO_ROOT,"metadata/GSE18995_pairs.csv")
HELPER_FILE <- SCRIPT_FILE
RESULTS_ROOT <- repo_output("Figure1/ML")
RUN_DIR <- repo_ml_run()# Empty: auto-find one run. Or paste your run_* directory here.
# The GSE18995 SHAP folder is used ONLY to reuse identical BACKGROUND pairs,
# never to select the model or explain external samples.
MAX_BACKGROUND_PAIRS <- 20L  # Complete CIT+R2H training donor pairs.
EXACT_BATCH_MASKS <- 32L
BOOTSTRAP_B <- 1000L
PERMUTATION_B <- 100L
SEED <- 123456L
COL_CIT <- "#468CCB"
COL_R2H <- "#CC4F4E"
COL_OTHER <- "#91B9D2"
COL_DARK <- "#264653"
# ----------------------------------------------------------------------------

req <- c(TRAIN_PAIR_FILE,TEST_PAIR_FILE,HELPER_FILE)
missing <- req[!file.exists(req)]
if (length(missing)) stop("Missing project files: ",paste(missing,collapse=", "),
                          ". Set working directory to your original ML project.")
load_ml_helpers()
for (pkg in c("ggplot2","pROC"))
  if (!requireNamespace(pkg,quietly=TRUE)) stop("Install required package: ",pkg)

if (!nzchar(RUN_DIR)) {
  paths <- list.files(RESULTS_ROOT,pattern="^05_selected_frozen_model\\.rds$",
                      full.names=TRUE,recursive=TRUE)
  paths <- paths[basename(dirname(paths)) != "checkpoints_external"]
  if (length(paths)!=1L) stop("Detected ",length(paths)," frozen models. ",
     "Set RUN_DIR explicitly to the folder containing 05_selected_frozen_model.rds")
  RUN_DIR <- dirname(paths)
}
MODEL_FILE <- file.path(RUN_DIR,"05_selected_frozen_model.rds")
TRAIN_PRED_FILE <- file.path(RUN_DIR,"04_selected_model_train_OOF_predictions.tsv")
EXT_PRED_FILE <- file.path(RUN_DIR,"08_SELECTED_model_GSE18995_predictions.tsv")
for (f in c(MODEL_FILE,TRAIN_PRED_FILE,EXT_PRED_FILE))
  if (!file.exists(f)) stop("Saved model/prediction file missing: ",f)
OUT <- file.path(RUN_DIR,"RF_SVM_FROZEN_EXACT_SHAP_9GENES")
EXT_SHAP_DIR <- file.path(RUN_DIR,"RF_SVM_FROZEN_EXACT_SHAP_9GENES")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

read_expr <- function(path) {
  d <- utils::read.delim(path,header=TRUE,sep="\t",row.names=1,check.names=FALSE,
                         stringsAsFactors=FALSE,quote="",comment.char="")
  if (!"Type" %in% names(d) || anyDuplicated(rownames(d)))
    stop("Missing Type or duplicated rownames in ",path)
  ty <- tolower(trimws(as.character(d$Type)))
  y <- ifelse(ty %in% c("0","cit","control"),0L,
       ifelse(ty %in% c("1","r2h","treat","hr"),1L,NA_integer_))
  if (anyNA(y)) stop("Unrecognized class values")
  x <- as.data.frame(d[,setdiff(names(d),"Type"),drop=FALSE],check.names=FALSE)
  x[] <- lapply(x,function(t) as.numeric(as.character(t)))
  if (anyNA(x) || any(!is.finite(as.matrix(x)))) stop("Invalid numeric matrix: ",path)
  list(X=x,y=y,id=rownames(d))
}
# Auto-detect expression files: 21-gene and 22-gene variants can contain
# the same GSM IDs but different expression VALUES. Never silently mix them.
model <- readRDS(MODEL_FILE)
if (!nzchar(TRAIN_FILE)) {
  candidates <- c("data.train.txt","data.train(3).txt")
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("Neither data.train.txt nor data.train(3).txt exists")
  scaler <- model$scaler
  train_error <- vapply(candidates,function(path) {
    tryCatch({
      d <- read_expr(path)
      if (length(d$y)!=106L || !all(scaler$kept %in% colnames(d$X))) return(Inf)
      xx <- as.matrix(d$X[,scaler$kept,drop=FALSE])
      m <- colMeans(xx); sdv <- apply(xx,2L,stats::sd)
      max(c(abs(m-as.numeric(scaler$center[scaler$kept])),
            abs(sdv-as.numeric(scaler$scale[scaler$kept]))))
    },error=function(e) Inf)
  },numeric(1))
  if (all(!is.finite(train_error))) stop("No training matrix matches the saved model's predictors")
  TRAIN_FILE <- candidates[which.min(train_error)]
  message("Using training matrix ",TRAIN_FILE," (scaler discrepancy=",
          signif(min(train_error),4),")")
  if (min(train_error)>1e-5)
    stop("TRAIN_FILE does not match saved model training-scaler statistics. ",
         "Check model run and data versions.")
}
if (!nzchar(TEST_FILE)) {
  candidates <- c("data.test(3).txt","data.test.txt")
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("Neither data.test(3).txt nor data.test.txt exists")
  reference <- utils::read.delim(EXT_PRED_FILE,check.names=FALSE,stringsAsFactors=FALSE)
  test_error <- vapply(candidates,function(path) {
    tryCatch({
      d <- read_expr(path)
      if (length(d$y)!=32L || !setequal(d$id,as.character(reference$Sample))) return(Inf)
      p <- predict_ml_pipeline(model,d$X)
      max(abs(p-as.numeric(reference$Score_R2H[match(d$id,reference$Sample)])))
    },error=function(e) Inf)
  },numeric(1))
  if (all(!is.finite(test_error))) stop("No 32-sample GSE18995 matrix matches saved external scores")
  TEST_FILE <- candidates[which.min(test_error)]
  message("Using external matrix ",TEST_FILE," (saved-score discrepancy=",
          signif(min(test_error),4),")")
  if (min(test_error)>1e-6)
    stop("TEST_FILE cannot reproduce frozen GSE18995 predictions. ",
         "Check input matrices/helper versions.")
}
for (path in c(TRAIN_FILE,TEST_FILE))
  if (!file.exists(path)) stop("Selected expression file missing: ",path)
train <- read_expr(TRAIN_FILE)
test <- read_expr(TEST_FILE)
train_donor <- read_verified_pairs(TRAIN_PAIR_FILE,train$id,train$y,53L)
test_donor <- read_verified_pairs(TEST_PAIR_FILE,test$id,test$y,16L)
stopifnot(length(train$y)==106L,length(test$y)==32L)
if (!identical(model$label,"RF+SVM") || !identical(model$learner$method,"SVM")) {
  stop("Frozen model is NOT RF+SVM. Found: ",model$label,
       " (learner=",model$learner$method,"). Do not mislabel another classifier.")
}
genes <- as.character(model$selected_genes)
if (length(genes)!=9L || anyDuplicated(genes)) {
  stop("RF+SVM currently has ",length(genes)," selected genes, not the nine ",
       "in the supplied ranking. Confirm this RUN_DIR belongs to the same run.")
}
if (!all(genes %in% names(train$X)) || !all(genes %in% names(test$X)))
  stop("Selected genes missing from input expression matrices")
message("FINAL classifier: RF+SVM | selected genes: ",paste(genes,collapse=", "))

train_oof <- utils::read.delim(TRAIN_PRED_FILE,check.names=FALSE)
ext_oof <- utils::read.delim(EXT_PRED_FILE,check.names=FALSE)
for (d in list(train_oof,ext_oof)) {
  if (!all(c("Sample","Donor","Type","Score_R2H") %in% colnames(d)))
    stop("Prediction tables require Sample, Donor, Type, Score_R2H")
}
# Match columns to original input order rather than relying on saved row order.
align_pred <- function(pred,input,donor) {
  ii <- match(input$id,as.character(pred$Sample))
  if (anyNA(ii) || anyDuplicated(pred$Sample) || nrow(pred)!=length(input$id))
    stop("Saved prediction table does not match matrix samples")
  pred <- pred[ii,,drop=FALSE]
  if (any(as.integer(pred$Type)!=input$y) ||
      any(as.character(pred$Donor)!=as.character(donor)))
    stop("Saved predictions disagree with verified donor/case labels")
  pred
}
train_oof <- align_pred(train_oof,train,train_donor)
ext_oof <- align_pred(ext_oof,test,test_donor)
p_external <- predict_ml_pipeline(model,test$X)
if (max(abs(p_external-as.numeric(ext_oof$Score_R2H)))>1e-6)
  stop("Frozen model predictions do not reproduce saved GSE18995 scores. ",
       "Check RUN_DIR/input data/helper versions!")

calc_auc <- function(y,p) as.numeric(pROC::auc(
  pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE)))
auc_oof <- calc_auc(train$y,as.numeric(train_oof$Score_R2H))
auc_ext <- calc_auc(test$y,p_external)
cat("Verified AUC: grouped OOF=",sprintf("%.6f",auc_oof),
    ", external=",sprintf("%.6f",auc_ext),"\n",sep="")

# Paired (clustered) bootstrap: resample DONORS, retaining both biopsies.
set.seed(SEED+11L)
ci_oof <- paired_auc_bootstrap(train$y,as.numeric(train_oof$Score_R2H),
                               train_donor,B=BOOTSTRAP_B,seed=SEED+11L)
ci_ext <- paired_auc_bootstrap(test$y,p_external,test_donor,
                               B=BOOTSTRAP_B,seed=SEED+12L)
auc_summary <- data.frame(
  Dataset=c("Training donor-grouped 5-fold OOF","GSE18995 external"),
  Donor_pairs=c(53,16),AUC=c(auc_oof,auc_ext),
  CI_low=c(ci_oof[1],ci_ext[1]),CI_high=c(ci_oof[2],ci_ext[2]))
write.csv(auc_summary,file.path(OUT,"RF_SVM_AUC_paired_donor_bootstrap.csv"),row.names=FALSE)
write.csv(data.frame(Gene=genes,Source="RF selector, final full-training refit"),
          file.path(OUT,"RF_SVM_final_9_features.csv"),row.names=FALSE)

# Explain the frozen FULL-training SVM on training samples (NOT OOF SHAP).
# Panel A ROC still uses saved donor-grouped TRAINING OOF predictions and
# GSE18995 external predictions. SHAP B-F explains full-training frozen model.
explain <- test
explain_donor <- test_donor
p_explain <- predict_ml_pipeline(model,explain$X)
if (length(p_explain) != 32L || any(!is.finite(p_explain)))
  stop("Could not compute 32 frozen-model test scores")
auc_train_apparent <- calc_auc(explain$y,p_explain)
cat("Test SHAP: 16 pairs, 32 samples; frozen-model fitted AUC=",
    sprintf("%.3f",auc_train_apparent),"; grouped OOF AUC=",
    sprintf("%.3f",auc_oof),"\n",sep="")

# Keep exactly the SAME training-donor background as completed external SHAP.
# If an external background CSV exists, follow its sample IDs IN SAVED ORDER.
external_bg_csv <- file.path(RUN_DIR,"RF_SVM_FROZEN_EXACT_SHAP_TRAIN_53PAIRS","SHAP_background_training_donor_pairs.csv")
if (file.exists(external_bg_csv)) {
  bg_tab <- utils::read.csv(external_bg_csv,stringsAsFactors=FALSE,check.names=FALSE)
  if (!all(c("Sample","Donor","Group") %in% names(bg_tab)) ||
      anyDuplicated(bg_tab$Sample) || nrow(bg_tab) %% 2L != 0L)
    stop("External SHAP background CSV is invalid")
  bg_idx <- match(as.character(bg_tab$Sample),train$id)
  if (anyNA(bg_idx) || anyDuplicated(bg_idx) ||
      any(as.character(bg_tab$Donor)!=as.character(train_donor[bg_idx])) ||
      any(as.character(bg_tab$Group)!=ifelse(train$y[bg_idx]==0,"CIT","R2H")))
    stop("Original external-SHAP background cannot be matched to this training matrix")
  chosen_bg_donors <- sort(unique(as.character(train_donor[bg_idx])))
  cat("Reusing EXTERNAL SHAP's original ",length(chosen_bg_donors),
      " training-background donor pairs\n",sep="")
} else {
  set.seed(SEED+21L)
  chosen_bg_donors <- sample(sort(unique(as.character(train_donor))),
                             min(MAX_BACKGROUND_PAIRS,length(unique(train_donor))))
  bg_idx <- which(as.character(train_donor) %in% chosen_bg_donors)
  warning("External SHAP background CSV not found; using deterministic same-seed selection. ",
          "For strict comparability, provide the earlier background CSV.")
}
if (length(bg_idx)!=2L*length(chosen_bg_donors) ||
    !all(table(train_donor[bg_idx],factor(train$y[bg_idx],levels=0:1))==1L))
  stop("Background does not contain complete CIT/R2H donor pairs")
# The frozen pipeline scales ALL kept training predictors BEFORE restricting
# to the RF-selected 9 genes. Therefore coalition data must contain the full
# scaler$kept set, even though only selected_genes enter the SHAP coalitions.
required_model_cols <- as.character(model$scaler$kept)
if (!length(required_model_cols) || anyNA(required_model_cols) ||
    anyDuplicated(required_model_cols) ||
    !all(required_model_cols %in% colnames(train$X)) ||
    !all(required_model_cols %in% colnames(test$X)))
  stop("Frozen model scaler$kept is incompatible with expression input columns")
if (!all(genes %in% required_model_cols))
  stop("RF-selected genes are not a subset of frozen scaler$kept")
background <- train$X[bg_idx,required_model_cols,drop=FALSE]
stopifnot(nrow(background)==2L*length(chosen_bg_donors))
write.csv(data.frame(Sample=train$id[bg_idx],Donor=train_donor[bg_idx],
                     Group=ifelse(train$y[bg_idx]==0,"CIT","R2H")),
          file.path(OUT,"SHAP_background_training_donor_pairs.csv"),row.names=FALSE)

# Ensure training-vs-external SHAP uses the same frozen model, background, and score scale.
# The interventional SHAP baseline equals the mean prediction on the background.
background_score <- mean(predict_ml_pipeline(model,background))
previous_ext_rds <- file.path(EXT_SHAP_DIR,"RF_SVM_exact_SHAP_result.rds")
if (file.exists(previous_ext_rds)) {
  old_ext <- readRDS(previous_ext_rds)
  if (!is.null(old_ext$baseline) &&
      (any(!is.finite(old_ext$baseline)) ||
       max(abs(old_ext$baseline-background_score))>1e-6)) {
    stop("Training SHAP background or expression matrix differs from external SHAP. ",
         "Check TRAIN_FILE and model versions before comparison.")
  }
  message("PASS: SHAP reference baseline matches previous external SHAP")
}

# EXACT interventional Shapley (2^9 coalitions).
# Value function: mean SVM 0-1 transformed margin after interventions on
# a training-donor-pair background. This is NOT a calibrated probability.
# Only selected gene columns are intervened; full scaler input is preserved.
# Inference of SHAP's values is descriptive; donors are the resampling units.
exact_one <- function(x,model,background,batch_masks=32L) {
  # Explain only the nine SELECTED genes; other scaler columns are retained
  # as required prediction inputs but never enumerated as SHAP features.
  genes <- as.character(model$selected_genes)
  if (!all(genes %in% colnames(background)))
    stop("Exact SHAP background is missing RF-selected genes")
  p <- length(genes)
  if (p>12L) stop("Exact SHAP intentionally limited to <=12 genes")
  x <- as.numeric(x[genes]); names(x) <- genes
  nmask <- 2L^p
  masks <- 0L:(nmask-1L)
  bits <- sapply(seq_len(p),function(j)
    bitwAnd(masks,bitwShiftL(1L,j-1L))!=0L)
  if (p==1L) bits <- matrix(bits,ncol=1L)
  K <- rowSums(bits)
  v <- numeric(nmask)
  B <- nrow(background)
  for (lo in seq.int(1L,nmask,by=as.integer(batch_masks))) {
    rows <- lo:min(nmask,lo+as.integer(batch_masks)-1L)
    data_b <- background[rep(seq_len(B),times=length(rows)),,drop=FALSE]
    for (h in seq_along(rows)) {
      use_cols <- which(bits[rows[h],])
      if (length(use_cols)) {
        rr <- ((h-1L)*B+1L):(h*B)
        # Match by gene NAME: selected genes need not be the first/adjacent
        # columns in the full scaler input matrix.
        for (j in use_cols) data_b[rr,genes[j]] <- x[j]
      }
    }
    pred <- predict_ml_pipeline(model,data_b)
    v[rows] <- colMeans(matrix(pred,nrow=B,ncol=length(rows)))
  }
  phi <- setNames(numeric(p),genes)
  for (j in seq_len(p)) {
    unmasked <- which(!bits[,j])
    w <- 1 / (p*choose(p-1,K[unmasked]))
    with_j <- masks[unmasked] + bitwShiftL(1L,j-1L) + 1L
    phi[j] <- sum(w*(v[with_j]-v[unmasked]))
  }
  list(shap=phi,baseline=v[1L],prediction=v[nmask],
       reconstruction=v[1L]+sum(phi))
}

# Resume SHAP by sample after interruption; refuse mismatched inputs/background.
signature <- c("TEST-16pairs-exact-shap-full-scaler-input-v3",
               unname(tools::md5sum(c(MODEL_FILE,TRAIN_FILE,HELPER_FILE))),
               paste(sort(bg_idx),collapse=","),paste(genes,collapse=";"),
               paste(required_model_cols,collapse=";"))
cache_file <- file.path(OUT,"CHECKPOINT_exact_SHAP.rds")
S <- matrix(NA_real_,nrow=length(explain$id),ncol=length(genes),
            dimnames=list(explain$id,genes))
baseline <- prediction <- rep(NA_real_,length(explain$id))
if (file.exists(cache_file)) {
  old <- tryCatch(readRDS(cache_file),error=function(e) NULL)
  if (!is.null(old) && identical(old$signature,signature) &&
      identical(dimnames(old$S),dimnames(S))) {
    S <- old$S; baseline <- old$baseline; prediction <- old$prediction
    cat("Resuming exact SHAP from checkpoint\n")
  } else stop("Existing SHAP checkpoint belongs to different data/model/background. ",
              "Rename/remove CHECKPOINT_exact_SHAP.rds before recomputing.")
}
for (i in seq_along(explain$id)) {
  if (all(is.finite(S[i,])) && is.finite(baseline[i])) next
  ans <- exact_one(explain$X[i,genes,drop=TRUE],model,background,
                   batch_masks=EXACT_BATCH_MASKS)
  if (abs(ans$prediction-p_explain[i])>1e-6 ||
      abs(ans$reconstruction-p_explain[i])>1e-6)
    stop("Exact SHAP additivity check failed at sample: ",explain$id[i])
  S[i,] <- ans$shap; baseline[i] <- ans$baseline; prediction[i] <- ans$prediction
  saveRDS(list(signature=signature,S=S,baseline=baseline,prediction=prediction),
          cache_file)
  cat(sprintf("Exact SHAP %2d/%2d: %s\n",i,length(explain$id),explain$id[i]))
}
if (max(abs(baseline+rowSums(S)-p_explain))>1e-6)
  stop("Final exact SHAP additivity mismatch")
saveRDS(list(S=S,baseline=baseline,prediction=prediction,
             model_label=model$label,selected_genes=genes,
             prediction_type="uncalibrated SVM margin transformed with sigmoid"),
        file.path(OUT,"RF_SVM_exact_SHAP_result.rds"))
write.csv(data.frame(Sample=explain$id,Donor=as.character(explain_donor),
                     Group=ifelse(explain$y==0,"CIT","R2H"),
                     Baseline=baseline,Predicted_Score=prediction,
                     Reconstruction=baseline+rowSums(S)),
          file.path(OUT,"RF_SVM_exact_SHAP_additivity.csv"),row.names=FALSE)

# Global importance: bootstrap COMPLETE TRAINING DONOR PAIRS (n=53).
set.seed(SEED+91L)
u <- sort(unique(as.character(explain_donor)))
imp_mat <- replicate(BOOTSTRAP_B,{
  donor_pick <- sample(u,length(u),replace=TRUE)
  ix <- unlist(lapply(donor_pick,function(d) which(explain_donor==d)),use.names=FALSE)
  colMeans(abs(S[ix,,drop=FALSE]))
})
imp <- data.frame(
  Gene=genes,MeanAbsSHAP=colMeans(abs(S)),
  CI_Lower=apply(imp_mat,1,stats::quantile,probs=.025,na.rm=TRUE),
  CI_Upper=apply(imp_mat,1,stats::quantile,probs=.975,na.rm=TRUE))
imp <- imp[order(-imp$MeanAbsSHAP),,drop=FALSE]
imp$Rank <- seq_len(nrow(imp))
write.csv(imp,file.path(OUT,"RF_SVM_SHAP_global_importance_53pair_CI.csv"),row.names=FALSE)

# Optional cross-cohort SHAP stability table (no refitting and no new inference).
ext_imp_file <- file.path(EXT_SHAP_DIR,"RF_SVM_SHAP_global_importance_16pair_CI.csv")
if (file.exists(ext_imp_file)) {
  ext_imp <- utils::read.csv(ext_imp_file,stringsAsFactors=FALSE,check.names=FALSE)
  cols <- c("Gene","MeanAbsSHAP","Rank","CI_Lower","CI_Upper")
  if (all(cols %in% colnames(ext_imp)) && setequal(ext_imp$Gene,imp$Gene)) {
    a <- imp[,cols,drop=FALSE]
    b <- ext_imp[,cols,drop=FALSE]
    names(a)[-1L] <- paste0("Train_",names(a)[-1L])
    names(b)[-1L] <- paste0("GSE18995_",names(b)[-1L])
    cmp <- merge(a,b,by="Gene",sort=FALSE)
    cmp <- cmp[order(cmp$Train_Rank),,drop=FALSE]
    write.csv(cmp,file.path(OUT,"RF_SVM_train_vs_external_SHAP_comparison.csv"),row.names=FALSE)
  }
}

# Each training sample appears once per selected gene in beeswarm/dependence data.
long <- do.call(rbind,lapply(genes,function(g) {
  raw <- as.numeric(explain$X[,g])
  ra <- range(raw)
  scaled <- if (diff(ra)==0) rep(.5,length(raw)) else (raw-ra[1])/diff(ra)
  data.frame(Sample=explain$id,Donor=as.character(explain_donor),
             Group=ifelse(explain$y==0,"CIT","R2H"),
             Gene=g,Expression=raw,ExpressionScaled=scaled,
             SHAP=S[,g],stringsAsFactors=FALSE)
}))
write.csv(long,file.path(OUT,"RF_SVM_SHAP_long_format.csv"),row.names=FALSE)

# Donor-paired gene expression analysis: 9 tests jointly BH corrected.
paired_vectors <- function(v,donor,y) {
  di <- sort(unique(as.character(donor)))
  cit <- r2h <- numeric(length(di))
  for (i in seq_along(di)) {
    a <- which(donor==di[i] & y==0L)
    b <- which(donor==di[i] & y==1L)
    if (length(a)!=1L || length(b)!=1L) stop("Unmatched donor: ",di[i])
    cit[i] <- v[a]; r2h[i] <- v[b]
  }
  list(donor=di,CIT=cit,R2H=r2h,Delta=r2h-cit)
}
p_raw <- sapply(genes,function(g) {
  q <- paired_vectors(as.numeric(explain$X[,g]),explain_donor,explain$y)
  if (all(q$Delta==0)) return(1)
  stats::wilcox.test(q$R2H,q$CIT,paired=TRUE,exact=FALSE)$p.value
})
p_adj <- stats::p.adjust(p_raw,method="BH")
paired_gene <- data.frame(Gene=genes,P_value=as.numeric(p_raw),
                          FDR_BH=as.numeric(p_adj))
write.csv(paired_gene,file.path(OUT,"GSE145989_9genes_paired_Wilcoxon_BH.csv"),row.names=FALSE)

# 9-gene paired coordinated program score from TRAINING mean/SD only.
explain_z <- transform_by_training(explain$X,model$scaler)[,genes,drop=FALSE]
score9 <- rowMeans(as.matrix(explain_z))
pair_score9 <- paired_vectors(score9,explain_donor,explain$y)
p_pair_score9 <- stats::wilcox.test(pair_score9$R2H,pair_score9$CIT,
                                    paired=TRUE,exact=FALSE)$p.value
write.csv(data.frame(Sample=explain$id,Donor=explain_donor,
                     Group=ifelse(explain$y==0,"CIT","R2H"),Score9=score9),
          file.path(OUT,"GSE18995_9gene_score.csv"),row.names=FALSE)

# Additional: paired DELTA expression correlation across 53 independent donors.
delta_expr <- sapply(genes,function(g)
  paired_vectors(as.numeric(explain$X[,g]),explain_donor,explain$y)$Delta)
colnames(delta_expr) <- genes
rho <- stats::cor(delta_expr,method="spearman",use="pairwise.complete.obs")
write.csv(rho,file.path(OUT,"GSE145989_9gene_paired_delta_Spearman.csv"))

# Training paired permutation importance (APPARENT, not out-of-fold):
# permute each gene across donor PAIRS
# (copy both CIT and R2H of the source donor together, by condition).
# This is an exploratory supplementary importance metric, not SHAP and not
# used for model selection.
set.seed(SEED+121L)
perm_drop <- matrix(NA_real_,nrow=PERMUTATION_B,ncol=length(genes),
                    dimnames=list(NULL,genes))
for (j in seq_along(genes)) {
  gene <- genes[j]
  for (b in seq_len(PERMUTATION_B)) {
    source_donors <- sample(u,length(u),replace=FALSE)
    xp <- explain$X
    for (t in seq_along(u)) {
      for (grp in 0:1) {
        ii <- which(explain_donor==u[t] & explain$y==grp)
        jj <- which(explain_donor==source_donors[t] & explain$y==grp)
        xp[ii,gene] <- explain$X[jj,gene]
      }
    }
    perm_drop[b,j] <- auc_train_apparent-calc_auc(explain$y,predict_ml_pipeline(model,xp))
  }
}
perm_imp <- data.frame(Gene=genes,MeanAUCdrop=colMeans(perm_drop),
                       Lower=apply(perm_drop,2,stats::quantile,.025),
                       Upper=apply(perm_drop,2,stats::quantile,.975))
perm_imp <- perm_imp[order(-perm_imp$MeanAUCdrop),]
write.csv(perm_imp,file.path(OUT,"RF_SVM_53pair_permutation_AUC_drop.csv"),row.names=FALSE)

# --------------------- Publication-style Figure 2 panels --------------------
library(ggplot2)
# Publication typography: all TEXT is bold black; original figure COLORS retained.
base_theme <- function(size=11) {
  theme_classic(base_size=size,base_family="sans") +
    theme(
      text=element_text(family="sans",face="bold",color="black"),
      plot.title=element_text(face="bold",color="black",hjust=.5,
                              size=size+1.8,margin=margin(b=9)),
      plot.subtitle=element_text(face="bold",color="black",hjust=.5,
                                 size=size-.4,margin=margin(b=8)),
      axis.title=element_text(face="bold",color="black",size=size+.4),
      axis.text=element_text(face="bold",color="black",size=size-.5),
      legend.title=element_text(face="bold",color="black",size=size-.3),
      legend.text=element_text(face="bold",color="black",size=size-.8),
      strip.text=element_text(face="bold",color="black"),
      plot.tag=element_text(face="bold",color="black"),
      axis.line=element_line(color="black",linewidth=.78),
      axis.ticks=element_line(color="black",linewidth=.65),
      panel.border=element_rect(color="black",fill=NA,linewidth=.72),
      legend.background=element_rect(fill="white",color=NA),
      legend.key=element_rect(fill="white",color=NA),
      plot.margin=margin(10,12,10,10))
}
save_panel <- function(p,name,width=7.3,height=5.5) {
  ggplot2::ggsave(file.path(OUT,paste0(name,".pdf")),p,
                  width=width,height=height,bg="white",useDingbats=FALSE)
  if (requireNamespace("ragg",quietly=TRUE)) {
    ggplot2::ggsave(file.path(OUT,paste0(name,".png")),p,
                    width=width,height=height,dpi=600,
                    device=ragg::agg_png,bg="white")
  } else {
    ggplot2::ggsave(file.path(OUT,paste0(name,".png")),p,
                    width=width,height=height,dpi=300,bg="white")
  }
}

roc_d <- function(y,p,name) {
  o <- pROC::roc(y,p,direction="<",levels=c(0,1),quiet=TRUE)
  data.frame(FPR=1-o$specificities,TPR=o$sensitivities,Curve=name)
}
roc_df <- rbind(roc_d(train$y,as.numeric(train_oof$Score_R2H),"Training grouped OOF"),
                roc_d(test$y,p_external,"GSE18995 external"))
leg <- c(sprintf("Training grouped OOF: AUC=%.3f",auc_oof),
         sprintf("GSE18995 external: AUC=%.3f",auc_ext))
names(leg) <- c("Training grouped OOF","GSE18995 external")
pA <- ggplot(roc_df,aes(FPR,TPR,color=Curve))+
  geom_path(linewidth=1.55,lineend="round",linejoin="round")+
  geom_abline(linetype=2,color="#999999",linewidth=.78)+
  scale_color_manual(values=c("Training grouped OOF"=COL_DARK,
                              "GSE18995 external"="#16A085"),labels=leg)+
  coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+
  labs(title="RF+SVM discriminatory performance",x="1 - Specificity",
       y="Sensitivity",color=NULL)+base_theme(10.5)+
  # Move the TWO legend rows from lower-right inward; retain line colors.
  theme(legend.position=c(.58,.23),legend.justification=c(.5,.5),
        legend.direction="vertical",
        legend.text=element_text(face="bold",color="black",size=9.0),
        legend.background=element_rect(fill=grDevices::adjustcolor("white",alpha.f=.92),color=NA))
save_panel(pA,"Panel_A_RF_SVM_ROC")

imp$Gene <- factor(imp$Gene,levels=rev(as.character(imp$Gene)))
pB <- ggplot(imp,aes(MeanAbsSHAP,Gene,fill=Gene=="SERPINE1"))+
  geom_col(width=.67)+
  geom_segment(aes(x=CI_Lower,xend=CI_Upper,y=Gene,yend=Gene),
               color=COL_DARK,linewidth=.65)+
  geom_point(shape=21,size=2,fill="white",color=COL_DARK)+
  scale_fill_manual(values=c(`TRUE`=COL_R2H,`FALSE`=COL_OTHER),guide="none")+
  scale_x_continuous(expand=expansion(mult=c(0,.22)))+
  geom_text(aes(x=CI_Upper,label=sprintf("%.3f",MeanAbsSHAP)),
            hjust=-.12,size=3.2,fontface="bold",color="black")+
  labs(title="Global exact SHAP importance",subtitle="GSE145989: 53-donor-pair bootstrap 95% CI",
       x="Mean |SHAP| on SVM decision score",y=NULL)+base_theme(10.5)
save_panel(pB,"Panel_B_RF_SVM_SHAP_global")

long$Gene <- factor(long$Gene,levels=rev(as.character(imp$Gene)))
if (requireNamespace("ggbeeswarm",quietly=TRUE)) {
  pt_layer <- ggbeeswarm::geom_quasirandom(width=.26,alpha=.82,size=1.9)
} else {
  set.seed(SEED)
  pt_layer <- geom_jitter(width=0,height=.18,alpha=.82,size=1.9)
}
pC <- ggplot(long,aes(SHAP,Gene,color=ExpressionScaled))+
  geom_vline(xintercept=0,linetype=2,color="#888888")+
  pt_layer+
  scale_color_gradientn(colors=c("#2166AC","#F7F7F7","#B2182B"),
                        limits=c(0,1),breaks=c(0,1),labels=c("Low","High"))+
  labs(title="Exact SHAP beeswarm",subtitle="Positive SHAP shifts score toward R2H",
       x="SVM score contribution (uncalibrated)",y=NULL,color="Expression")+
  base_theme(10.5)
save_panel(pC,"Panel_C_RF_SVM_SHAP_beeswarm")

ser <- long[as.character(long$Gene)=="SERPINE1",,drop=FALSE]
ser$Group <- factor(ser$Group,levels=c("CIT","R2H"))
pD <- ggplot(ser,aes(Expression,SHAP,color=Group))+
  geom_hline(yintercept=0,linetype=2,color="#999999")+
  geom_point(size=2.5,alpha=.85)+
  geom_smooth(aes(group=1),method="loess",formula=y~x,se=TRUE,
              inherit.aes=TRUE,color=COL_DARK,linewidth=.85)+
  scale_color_manual(labels=c(CIT="CIT",R2H="HR"),values=c(CIT=COL_CIT,R2H=COL_R2H))+
  labs(title="SERPINE1 SHAP dependence",subtitle="Exploratory model explanation",
       x="SERPINE1 expression",y="Exact SHAP (SVM decision score)",color=NULL)+
  base_theme(10.5)
save_panel(pD,"Panel_D_SERPINE1_dependence")

citser <- paired_vectors(as.numeric(explain$X[,"SERPINE1"]),explain_donor,explain$y)
paired_ser <- data.frame(Donor=rep(citser$donor,2L),
                         Group=factor(rep(c("CIT","R2H"),each=length(citser$donor)),
                                      levels=c("CIT","R2H")),
                         Expression=c(citser$CIT,citser$R2H))
ser_fdr <- paired_gene$FDR_BH[paired_gene$Gene=="SERPINE1"]
pE <- ggplot(paired_ser,aes(Group,Expression,group=Donor))+
  geom_line(color="#A9A9A9",linewidth=.5,alpha=.7)+
  geom_point(aes(color=Group),size=2.9,alpha=.93)+
  scale_color_manual(labels=c(CIT="CIT",R2H="HR"),values=c(CIT=COL_CIT,R2H=COL_R2H),guide="none")+
  labs(title="Training paired SERPINE1 expression",
       subtitle=sprintf("53 donor pairs | Wilcoxon BH-FDR=%.4g",ser_fdr),
       x=NULL,y="Normalized expression")+base_theme(10.5)
save_panel(pE,"Panel_E_GSE145989_SERPINE1_paired")

# Typical local case chosen by MEDIAN final-model score among 53 R2H TRAIN samples.
rr <- which(explain$y==1)
rr <- rr[order(abs(p_explain[rr]-stats::median(p_explain[rr])))]
k <- rr[1L]
phi <- S[k,]
phi <- phi[order(abs(phi),decreasing=TRUE)]
starts <- c(baseline[k],baseline[k]+head(cumsum(phi),-1L))
ends <- starts+phi
wf <- data.frame(Gene=factor(names(phi),levels=names(phi)),
                 i=seq_along(phi),Start=starts,End=ends,
                 Contribution=as.numeric(phi),
                 Direction=ifelse(phi>=0,"Toward HR","Toward CIT"))
pF <- ggplot(wf)+
  geom_rect(aes(xmin=i-.36,xmax=i+.36,
                ymin=pmin(Start,End),ymax=pmax(Start,End),fill=Direction))+
  geom_segment(data=utils::head(wf,-1L),
               aes(x=i+.36,xend=i+1-.36,y=End,yend=End),linewidth=.5,color=COL_DARK)+
  geom_hline(yintercept=baseline[k],linetype=2,color="#888888")+
  geom_hline(yintercept=p_explain[k],linetype=3,color=COL_DARK)+
  geom_text(aes(x=i,y=pmax(Start,End)+.02,label=sprintf("%+.3f",Contribution)),
            size=3,fontface="bold",color="black")+
  scale_x_continuous(breaks=wf$i,labels=as.character(wf$Gene))+
  scale_fill_manual(values=c("Toward HR"=COL_R2H,"Toward CIT"=COL_CIT),guide="none")+
  scale_y_continuous(limits=c(min(0,min(starts,ends)-.08),
                              max(1.05,max(starts,ends)+.12)))+
  labs(title="Local exact SHAP waterfall",
       subtitle=paste0("Training HR median-score sample: ",explain$id[k],
         " | baseline=",sprintf("%.3f",baseline[k]),
         ", score=",sprintf("%.3f",p_explain[k])),
       x=NULL,y="SVM decision score (uncalibrated)")+
  base_theme(10)+theme(axis.text.x=element_text(angle=50,hjust=1,
                                                     face="bold",color="black"))
save_panel(pF,"Panel_F_RF_SVM_local_waterfall",width=9,height=5.5)

# Supplement: 9-gene paired score.
sc_df <- data.frame(Donor=rep(pair_score9$donor,2),
                    Group=factor(rep(c("CIT","R2H"),each=length(pair_score9$donor)),
                                 levels=c("CIT","R2H")),
                    Score=c(pair_score9$CIT,pair_score9$R2H))
pS1 <- ggplot(sc_df,aes(Group,Score,group=Donor))+
  geom_line(color="grey70",alpha=.7)+
  geom_point(aes(color=Group),size=2.6)+
  scale_color_manual(labels=c(CIT="CIT",R2H="HR"),values=c(CIT=COL_CIT,R2H=COL_R2H),guide="none")+
  labs(title="RF-selected nine-gene coordinated score",
       subtitle=sprintf("GSE145989 donor-paired Wilcoxon P=%.4g",p_pair_score9),
       x=NULL,y="Mean training-referenced z-score")+base_theme(10.5)
save_panel(pS1,"Supplement_1_GSE145989_9gene_paired_score")

pS2 <- ggplot(perm_imp,aes(MeanAUCdrop,reorder(Gene,MeanAUCdrop),fill=Gene=="SERPINE1"))+
  geom_col(width=.7)+
  geom_segment(aes(x=Lower,xend=Upper,y=Gene,yend=Gene),linewidth=.7,color=COL_DARK)+
  scale_fill_manual(values=c(`TRUE`=COL_R2H,`FALSE`=COL_OTHER),guide="none")+
  labs(title="Training paired-donor permutation importance",
       subtitle="Exploratory: mean AUC reduction from donor-pair exchange",
       x="AUC decrease after paired permutation",y=NULL)+base_theme(10.5)
save_panel(pS2,"Supplement_2_GSE145989_paired_permutation_importance")

if (requireNamespace("patchwork",quietly=TRUE)) {
  fig <- (pA+pB+pC)/(pD+pE+pF)+
    patchwork::plot_layout(guides="keep",heights=c(1,1))
  save_panel(fig,"Figure2_RF_SVM_SHAP_A_to_F",width=19,height=12)
}
writeLines(c(
  "RF+SVM frozen exact interventional SHAP (9 final RF-selected genes)",
  "SVM output: sigmoid-transformed oriented decision margin; uncalibrated (NOT a risk probability)",
  "Main-model selection was based on GSE145989 donor-grouped OOF AUC; no external-based re-selection",
  sprintf("RF+SVM OOF AUC %.6f; GSE18995 external AUC %.6f",auc_oof,auc_ext),
  "GSE145989 SHAP values refer to the final full-training SVM, not to the fold-specific CV models",
  "SHAP interventional background: balanced CIT/R2H biopsies from complete training-donor pairs",
  "SHAP importance CI: bootstrap of 53 complete training donor pairs (conditional on frozen fit)",
  "D: SERPINE1 expression vs SHAP is an exploratory model dependence, not causal evidence",
  "CIT/R2H inference: training donor-paired Wilcoxon; nine individual gene tests BH adjusted",
  sprintf("Final-model fitted training AUC %.6f (SHAP explanation); do NOT report as OOF AUC",auc_train_apparent),
  "Potential limitation: if the upstream candidate gene panel was selected using all 53 training donors,",
  "the 5-fold machine-learning OOF AUC does not provide full pipeline-nested gene-discovery validation",
  "Training SHAP is an IN-SAMPLE final-model explanation; it is NOT grouped-OOF SHAP"),
  con=file.path(OUT,"READ_THIS_Methods_and_interpretation.txt"))
cat("\nSHAP finished; all output: ",normalizePath(OUT),"\n",sep="")

}

# =============================================================================
# STAGE figure2_nine_gene_scores
# =============================================================================
PIPELINE$figure2_nine_gene_scores <- function() {
# English consolidated source; frozen-model inference, no result replacement.
# =============================================================================
# RF-selected nine-gene coordinated score | PAIRED VIOLIN, TWO COHORTS
# Training: GSE145989 (53 donor pairs); external: GSE18995 (16 donor pairs)
# Score = mean of 9 gene z-scores standardized with TRAINING center and SD.
# NOT the SVM predicted probability / decision score / SHAP importance.
# Replot from completed SHAP score CSVs; if absent, reconstruct score with
# EXACT original frozen model scaler (never refit / never recompute SHAP).
# Uses the original CIT blue and R2H red, bold black typography.
# =============================================================================

options(stringsAsFactors = FALSE)

# ----- USER CONFIGURATION: No need to edit if directories remain unchanged -----
PROJECT_DIR <- file.path(REPO_ROOT,"metadata")
RUN_DIR <- repo_ml_run()
# If moved, change PROJECT_DIR / RUN_DIR above, NOT the SHAP score formula.
TRAIN_SCORE_FILE <- file.path(RUN_DIR,
  "RF_SVM_FROZEN_EXACT_SHAP_TRAIN_53PAIRS", "GSE145989_9gene_score.csv")
EXTERNAL_SCORE_FILE <- file.path(RUN_DIR,
  "RF_SVM_FROZEN_EXACT_SHAP_9GENES", "GSE18995_9gene_score.csv")
OUTPUT_DIR <- file.path(RUN_DIR, "RF_SVM_NINE_GENE_PAIRED_VIOLIN")
COL_CIT <- "#468CCB"
COL_R2H <- "#CC4F4E"
SEED <- 123456L
# -----------------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop('Please install ggplot2: install.packages("ggplot2")')
if (!dir.exists(RUN_DIR)) stop("RUN_DIR does not exist:\n", RUN_DIR)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Source sample-to-donor verification metadata, not an inferred GSM pairing.
TRAIN_META_FILE <- file.path(REPO_ROOT,"metadata/GSE145989_pairs.csv")
EXT_META_FILE <- file.path(REPO_ROOT,"metadata/GSE18995_pairs.csv")
MODEL_FILE <- file.path(RUN_DIR, "05_selected_frozen_model.rds")

read_meta <- function(path, n_expected) {
  if (!file.exists(path)) {
    message("Metadata absent; checking complete donor pairs from score CSV only: ", path)
    return(NULL)
  }
  m <- utils::read.csv(path, fileEncoding = "UTF-8-BOM", check.names = FALSE)
  wanted <- c("donor_id", "CIT_sample_ID", "R2H_sample_ID")
  if (!all(wanted %in% names(m)))
    stop("Metadata missing columns: ", paste(setdiff(wanted, names(m)), collapse = ", "))
  if (nrow(m) != n_expected || anyDuplicated(m$donor_id) ||
      anyDuplicated(c(m$CIT_sample_ID, m$R2H_sample_ID)))
    stop("Metadata is not the expected ", n_expected, " complete distinct pairs: ", path)
  m
}

meta_train <- read_meta(TRAIN_META_FILE, 53L)
meta_ext <- read_meta(EXT_META_FILE, 16L)

# Optional score reconstruction only if a completed score CSV does not exist.
# Strictly choose training matrix by matching frozen scaler mean/SD; NEVER
# assume the similarly named 21-/22-gene files are interchangeable.
read_expression <- function(path, n_expected) {
  d <- utils::read.delim(path, sep = "\t", header = TRUE, row.names = 1,
                         check.names = FALSE, quote = "", comment.char = "")
  if (nrow(d) != 2L * n_expected || !"Type" %in% names(d) ||
      anyDuplicated(rownames(d)))
    stop("Unexpected number of samples or Type column in ", path)
  ty <- tolower(trimws(as.character(d$Type)))
  group <- ifelse(ty %in% c("0", "cit", "control"), "CIT",
           ifelse(ty %in% c("1", "r2h", "treat", "hr"), "R2H", NA_character_))
  if (anyNA(group) || sum(group == "CIT") != n_expected ||
      sum(group == "R2H") != n_expected) stop("Invalid class assignment in ", path)
  X <- d[, setdiff(names(d), "Type"), drop = FALSE]
  X[] <- lapply(X, function(a) as.numeric(as.character(a)))
  if (any(!is.finite(as.matrix(X)))) stop("Invalid numerical expression in ", path)
  list(X = X, Sample = rownames(d), Group = group)
}

# Score CSV reconstruction needs original frozen model; do not search among
# multiple 05_selected_frozen_model.rds copies (previous ambiguous-model error).
rebuild_missing_scores <- function(need_train, need_ext) {
  if (!file.exists(MODEL_FILE))
    stop("Missing score CSV, and frozen model unavailable: ", MODEL_FILE,
         "\nPlease point RUN_DIR to the original RF+SVM run.")
  model <- readRDS(MODEL_FILE)
  genes <- as.character(model$selected_genes)
  sc <- model$scaler
  if (!identical(model$label, "RF+SVM") ||
      !identical(model$learner$method, "SVM") || length(genes) != 9L ||
      anyDuplicated(genes) || !all(genes %in% sc$kept))
    stop("The selected frozen model is not the original 9-gene RF+SVM model")
  if (!all(genes %in% names(sc$center)) ||
      !all(genes %in% names(sc$scale)))
    stop("Frozen model is missing training-scaled gene parameters")
  z_score <- function(a) {
    xx <- as.matrix(a$X[, genes, drop = FALSE]); storage.mode(xx) <- "double"
    z <- sweep(sweep(xx, 2L, sc$center[genes], "-"), 2L, sc$scale[genes], "/")
    if (any(!is.finite(z))) stop("Non-finite z-scores")
    rowMeans(z)
  }
  mk_df <- function(a, m, n) {
    if (is.null(m)) stop("Paired donor metadata required for reconstructing scores")
    e <- data.frame(Sample = c(m$CIT_sample_ID, m$R2H_sample_ID),
                    Donor = rep(as.character(m$donor_id), 2L),
                    Group = rep(c("CIT", "R2H"), each = n),
                    stringsAsFactors = FALSE)
    ix <- match(a$Sample, e$Sample)
    if (anyNA(ix) || anyDuplicated(a$Sample) ||
        !identical(as.character(a$Group), as.character(e$Group[ix])))
      stop("Expression file is inconsistent with previously verified sample pairs")
    data.frame(Sample = a$Sample, Donor = e$Donor[ix],
               Group = a$Group, Score9 = z_score(a))
  }
  # Choose one exact training matrix; if several match, lowest mean/SD error.
  train_candidates <- repo_input("ml_train")
  train_candidates <- train_candidates[file.exists(train_candidates)]
  scores <- vapply(train_candidates, function(path) {
    tryCatch({
      a <- read_expression(path, 53L)
      if (!all(sc$kept %in% names(a$X))) return(Inf)
      x <- as.matrix(a$X[, sc$kept, drop = FALSE])
      max(abs(colMeans(x) - sc$center[sc$kept]),
          abs(apply(x, 2L, stats::sd) - sc$scale[sc$kept]))
    }, error = function(e) Inf)
  }, numeric(1))
  if (!length(scores) || all(!is.finite(scores)))
    stop("No training expression matrix compatible with the frozen scaler")
  pick_train <- which.min(scores)
  if (scores[pick_train] > 1e-5)
    stop("Training matrix does not match the frozen model scaler. Best discrepancy: ",
         signif(scores[pick_train], 5), ". Do not mix score versions.")
  message("Score reconstruction: training matrix: ", train_candidates[pick_train])
  if (need_train) {
    a <- read_expression(train_candidates[pick_train], 53L)
    dat <- mk_df(a, meta_train, 53L)
    dir.create(dirname(TRAIN_SCORE_FILE), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(dat, TRAIN_SCORE_FILE, row.names = FALSE)
    message("Reconstructed missing TRAIN score CSV: ", TRAIN_SCORE_FILE)
  }
  if (need_ext) {
    test_candidates <- repo_input("ml_test")
    test_candidates <- test_candidates[file.exists(test_candidates)]
    ds <- lapply(test_candidates, function(f)
      tryCatch(read_expression(f, 16L), error = function(e) NULL))
    eligible <- which(vapply(ds, function(a) !is.null(a) &&
                        all(sc$kept %in% names(a$X)), logical(1)))
    if (!length(eligible)) stop("No eligible 32-sample external expression matrix")
    # If >1 possible test matrix, verify against original frozen predictions.
    chosen <- eligible
    pred_path <- file.path(RUN_DIR, "08_SELECTED_model_GSE18995_predictions.tsv")
    if (file.exists(pred_path)) {
      load_ml_helpers()
      old <- utils::read.delim(pred_path, check.names = FALSE)
      if (!all(c("Sample", "Score_R2H") %in% names(old)) ||
          anyDuplicated(old$Sample)) stop("Invalid original external predictions")
      diffs <- vapply(eligible, function(i) tryCatch({
        a <- ds[[i]]
        ix <- match(a$Sample, old$Sample)
        if (anyNA(ix) || length(ix) != nrow(old)) return(Inf)
        p <- predict_ml_pipeline(model, a$X)
        max(abs(p - as.numeric(old$Score_R2H[ix])))
      }, error = function(e) Inf), numeric(1))
      if (all(!is.finite(diffs)) || min(diffs) > 1e-6)
        stop("No test matrix reproduces the frozen GSE18995 predictions")
      chosen <- eligible[which.min(diffs)]
    } else if (length(eligible) != 1L) {
      stop("Multiple 32-sample external matrices; cannot choose without saved predictions")
    }
    a <- ds[[chosen[1]]]
    message("Score reconstruction: external matrix: ", test_candidates[chosen[1]])
    dat <- mk_df(a, meta_ext, 16L)
    dir.create(dirname(EXTERNAL_SCORE_FILE), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(dat, EXTERNAL_SCORE_FILE, row.names = FALSE)
    message("Reconstructed missing EXTERNAL score CSV: ", EXTERNAL_SCORE_FILE)
  }
}

need_train <- !file.exists(TRAIN_SCORE_FILE)
need_ext <- !file.exists(EXTERNAL_SCORE_FILE)
if (need_train || need_ext) {
  message("Missing saved score CSV(s). Reconstructing from the frozen TRAINING scaler.")
  rebuild_missing_scores(need_train, need_ext)
}

# Read and verify the actual saved scores. This is NOT a model prediction.
read_paired_score <- function(path, m, n, cohort) {
  d <- utils::read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  req <- c("Sample", "Donor", "Group", "Score9")
  if (!all(req %in% names(d)))
    stop(cohort, " score CSV missing: ", paste(setdiff(req, names(d)), collapse = ", "))
  if (nrow(d) != 2L * n || anyDuplicated(d$Sample) ||
      anyNA(d[, req]) || any(!is.finite(as.numeric(d$Score9))))
    stop(cohort, ": expected ", 2*n, " distinct finite score rows")
  d$Sample <- as.character(d$Sample)
  d$Donor <- as.character(d$Donor)
  d$Group <- as.character(d$Group)
  d$Score9 <- as.numeric(d$Score9)
  if (is.null(m)) {
    if (length(unique(d$Donor)) != n ||
        !all(table(d$Donor, factor(d$Group, levels = c("CIT", "R2H"))) == 1L))
      stop(cohort, " has incomplete / duplicated donor pairs")
  } else {
    expected <- data.frame(
      Sample = c(as.character(m$CIT_sample_ID), as.character(m$R2H_sample_ID)),
      Donor = rep(as.character(m$donor_id), 2L),
      Group = rep(c("CIT", "R2H"), each = n))
    ix <- match(d$Sample, expected$Sample)
    if (anyNA(ix) || !setequal(d$Sample, expected$Sample) ||
        !identical(d$Donor, as.character(expected$Donor[ix])) ||
        !identical(d$Group, as.character(expected$Group[ix])))
      stop(cohort, " score CSV does not match donor-pair provenance metadata")
  }
  if (length(unique(d$Donor)) != n ||
      any(table(d$Donor, factor(d$Group, levels = c("CIT", "R2H"))) != 1L))
    stop(cohort, " is not fully paired")
  d$Group <- factor(d$Group, levels = c("CIT", "R2H"))
  d$Cohort <- cohort
  d
}

train_df <- read_paired_score(TRAIN_SCORE_FILE, meta_train, 53L, "GSE145989")
ext_df <- read_paired_score(EXTERNAL_SCORE_FILE, meta_ext, 16L, "GSE18995")

pair_table <- function(d) {
  a <- d[d$Group == "CIT", c("Donor", "Score9")]
  b <- d[d$Group == "R2H", c("Donor", "Score9")]
  names(a)[2] <- "CIT"
  names(b)[2] <- "R2H"
  a <- merge(a, b, by = "Donor", sort = TRUE)
  a$Delta_R2H_minus_CIT <- a$R2H - a$CIT
  a
}
pair_train <- pair_table(train_df)
pair_ext <- pair_table(ext_df)
paired_p <- function(a) {
  if (all(a$Delta_R2H_minus_CIT == 0)) return(1)
  stats::wilcox.test(a$R2H, a$CIT, paired = TRUE, exact = FALSE)$p.value
}
p_train <- paired_p(pair_train)
p_ext <- paired_p(pair_ext)
# Explicit FDR family = TWO independent cohort-level comparisons,
# not the nine within-gene tests and not the Wilcoxon test on each sample.
p_fdr <- stats::p.adjust(c(p_train, p_ext), method = "BH")

summary_df <- data.frame(
  Dataset = c("GSE145989", "GSE18995"),
  Cohort = c("Training", "External validation"),
  Donor_pairs = c(53L, 16L),
  Paired_Wilcoxon_P = c(p_train, p_ext),
  BH_FDR_across_two_cohorts = p_fdr,
  Median_R2H_minus_CIT = c(median(pair_train$Delta_R2H_minus_CIT),
                            median(pair_ext$Delta_R2H_minus_CIT)),
  Mean_R2H_minus_CIT = c(mean(pair_train$Delta_R2H_minus_CIT),
                          mean(pair_ext$Delta_R2H_minus_CIT))
)
utils::write.csv(summary_df, file.path(OUTPUT_DIR, "NineGene_paired_score_stat_summary.csv"),
                 row.names = FALSE)
utils::write.csv(train_df, file.path(OUTPUT_DIR, "GSE145989_9gene_score_for_plot.csv"),
                 row.names = FALSE)
utils::write.csv(ext_df, file.path(OUTPUT_DIR, "GSE18995_9gene_score_for_plot.csv"),
                 row.names = FALSE)
utils::write.csv(pair_train, file.path(OUTPUT_DIR, "GSE145989_53paired_score_deltas.csv"),
                 row.names = FALSE)
utils::write.csv(pair_ext, file.path(OUTPUT_DIR, "GSE18995_16paired_score_deltas.csv"),
                 row.names = FALSE)

fmt_p <- function(x) {
  if (x < 0.001) format(x, scientific = TRUE, digits = 3)
  else sprintf("%.3f", x)
}

make_paired_violin <- function(d, pair_dat, raw_p, adj_p, subtitle, n_pairs) {
  d$Group_x <- ifelse(d$Group == "CIT", 1, 2)
  set.seed(SEED + n_pairs)
  jit <- stats::setNames(stats::runif(nrow(pair_dat), -.073, .073), pair_dat$Donor)
  d$Xpair <- d$Group_x + unname(jit[as.character(d$Donor)])
  pair_dat$Xcit <- 1 + unname(jit[as.character(pair_dat$Donor)])
  pair_dat$Xr2h <- 2 + unname(jit[as.character(pair_dat$Donor)])
  rng <- range(d$Score9)
  span <- diff(rng)
  if (!is.finite(span) || span <= 0) span <- max(0.5, abs(rng[2]) * .15)
  yline <- rng[2] + .080 * span
  ytick <- yline - .025 * span
  ytext <- rng[2] + .143 * span
  annotation <- paste0("Paired Wilcoxon P = ", fmt_p(raw_p),
                       "   |   BH-FDR = ", fmt_p(adj_p))
  alpha_pair <- if (n_pairs > 25L) .30 else .48

  ggplot2::ggplot() +
    ggplot2::geom_violin(data = d,
      ggplot2::aes(x = Group_x, y = Score9, group = Group, fill = Group),
      width = .79, trim = FALSE, adjust = 1.05,
      alpha = .50, color = "black", linewidth = .88) +
    ggplot2::geom_segment(data = pair_dat,
      ggplot2::aes(x = Xcit, xend = Xr2h, y = CIT, yend = R2H),
      inherit.aes = FALSE, linewidth = .50, alpha = alpha_pair,
      color = "#545454") +
    ggplot2::geom_boxplot(data = d,
      ggplot2::aes(x = Group_x, y = Score9, group = Group),
      width = .16, outlier.shape = NA, fill = "white", alpha = .87,
      color = "black", linewidth = .75) +
    ggplot2::geom_point(data = d,
      ggplot2::aes(x = Xpair, y = Score9, color = Group),
      alpha = if (n_pairs > 25L) .83 else .92,
      size = if (n_pairs > 25L) 1.8 else 2.25) +
    ggplot2::annotate("segment", x = 1, xend = 2,
      y = yline, yend = yline, color = "black", linewidth = .87) +
    ggplot2::annotate("segment", x = c(1, 2), xend = c(1, 2),
      y = c(ytick, ytick), yend = c(yline, yline),
      color = "black", linewidth = .87) +
    ggplot2::annotate("text", x = 1.5, y = ytext, label = annotation,
      color = "black", fontface = "bold", size = 3.35) +
    ggplot2::scale_x_discrete(labels=c(CIT="CIT",R2H="HR")) +
    ggplot2::scale_fill_manual(labels=c(CIT="CIT",R2H="HR"), values = c(CIT = COL_CIT, R2H = COL_R2H),
                               guide = "none") +
    ggplot2::scale_color_manual(labels=c(CIT="CIT",R2H="HR"), values = c(CIT = COL_CIT, R2H = COL_R2H),
                                guide = "none") +
    ggplot2::scale_x_continuous(breaks = c(1, 2), labels = c("CIT", "R2H"),
      limits = c(.48, 2.52), expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(limits = c(rng[1] - .12 * span,
                                            rng[2] + .225 * span),
      expand = ggplot2::expansion(mult = 0)) +
    ggplot2::labs(
      title = "RF-selected nine-gene coordinated score",
      subtitle = paste0(subtitle, "  |  ", n_pairs, " paired donors"),
      x = NULL, y = "Mean training-referenced z-score") +
    ggplot2::theme_classic(base_size = 11, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", face = "bold", color = "black"),
      plot.title = ggplot2::element_text(face = "bold", color = "black",
                          size = 13, hjust = .5, margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_text(face = "bold", color = "black",
                          size = 10.5, hjust = .5, margin = ggplot2::margin(b = 9)),
      axis.title = ggplot2::element_text(face = "bold", color = "black", size = 11.5),
      axis.text = ggplot2::element_text(face = "bold", color = "black", size = 11),
      axis.line = ggplot2::element_line(color = "black", linewidth = .85),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = .7),
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = .85),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(12, 18, 12, 14),
      legend.position = "none")
}

p_train_plot <- make_paired_violin(train_df, pair_train, p_train, p_fdr[1],
                                  "GSE145989 | Training", 53L)
p_ext_plot <- make_paired_violin(ext_df, pair_ext, p_ext, p_fdr[2],
                                "GSE18995 | External validation", 16L)

save_both <- function(p, stem) {
  f <- file.path(OUTPUT_DIR, stem)
  ggplot2::ggsave(paste0(f, ".pdf"), plot = p, width = 7.4, height = 5.6,
                  units = "in", bg = "white", useDingbats = FALSE)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(paste0(f, ".png"), plot = p, width = 7.4, height = 5.6,
                    units = "in", dpi = 600, device = ragg::agg_png, bg = "white")
  } else {
    ggplot2::ggsave(paste0(f, ".png"), plot = p, width = 7.4, height = 5.6,
                    units = "in", dpi = 600, bg = "white")
  }
  message("Saved PDF/PNG: ", f)
}
save_both(p_train_plot, "Training_GSE145989_53pairs_9gene_paired_violin")
save_both(p_ext_plot, "Validation_GSE18995_16pairs_9gene_paired_violin")

print(summary_df)
print(p_train_plot)
print(p_ext_plot)
message("COMPLETE: ", normalizePath(OUTPUT_DIR, winslash = "/"))
message("NOTE: BH-FDR adjusts the TWO cohort-wise paired-score comparisons only.")

}

# =============================================================================
# STAGE figure2_nine_gene_function
# =============================================================================
PIPELINE$figure2_nine_gene_function <- function() {
# English consolidated source; frozen-model inference, no result replacement.
# =============================================================================
# FIGURE F-G-H | GSE145989 TRAINING COHORT | FINAL RF+SVM 9-GENE PANEL
# Original layout/style adapted from Final_GSE145989_6Gene_FGH_VIVID_CIT_HR(4).R
#
# F = donor-paired 9-gene coordinated score (violin + box + true paired lines)
# G = exploratory GO BP and Reactome overrepresentation analysis (ORA)
# H = Spearman correlations of within-donor R2H-CIT expression changes
#
# USES ALREADY COMPUTED TRAINING SHAP TABLES, NOT geneexp(8).xlsx (6 genes).
# No model refitting, no recalculation of 9-gene scores, no SHAP recomputation.
# =============================================================================
options(stringsAsFactors = FALSE)
set.seed(20261021)

# --------------------------- 0. PATHS / STYLE -------------------------------
RUN_DIR <- repo_ml_run()

# Automatically try the user's renamed folder first, then the original folder.
TRAIN_DIR_CANDIDATES <- file.path(RUN_DIR, c(
  " shap", "RF_SVM_FROZEN_EXACT_SHAP_TRAIN_53PAIRS"
))
SCORE_BASENAME <- "GSE145989_9gene_score.csv"
LONG_BASENAME  <- "RF_SVM_SHAP_long_format.csv"
valid_dir <- vapply(TRAIN_DIR_CANDIDATES, function(d) {
  file.exists(file.path(d, SCORE_BASENAME)) &&
    file.exists(file.path(d, LONG_BASENAME))
}, logical(1))
if (!any(valid_dir)) {
  stop("Cannot find BOTH training-score and SHAP-long CSV files.\nChecked:\n",
       paste(TRAIN_DIR_CANDIDATES, collapse = "\n"),
       "\nExpected: ", SCORE_BASENAME, " AND ", LONG_BASENAME)
}
TRAIN_SHAP_DIR <- TRAIN_DIR_CANDIDATES[which(valid_dir)[1L]]
SCORE_FILE <- file.path(TRAIN_SHAP_DIR, SCORE_BASENAME)
LONG_FILE <- file.path(TRAIN_SHAP_DIR, LONG_BASENAME)
OUTPUT_DIR <- file.path(RUN_DIR, "FGH_TRAIN_GSE145989_RF9_PAIRED")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Match original uploaded VIVID F-G-H script: blue and red; black bold fonts.
COL_CIT <- "#2878D4"
COL_R2H <- "#E3322B"
BASE_FONT <- "sans"
EXPORT_WIDTH <- 18.5
EXPORT_HEIGHT <- 7.2
N_SHOW_GO <- 6L
N_SHOW_REACTOME <- 6L

# Exact selected 9 features from the frozen RF+SVM analysis.
SELECTED_GENES <- c(
  "ANGPTL4", "CCL20", "CXCL1", "CXCL3", "EREG",
  "ICAM1", "IL1A", "IL1B", "SERPINE1"
)

# This was the pre-selection gene universe in the 14-gene ML input matrix.
# Using all human genes as background after selecting from 14 candidates would
# make ORA less defensible; these nine genes were selected from this candidate set.
ML_CANDIDATE_UNIVERSE <- setdiff(names(read.delim(repo_input("ml_train"),row.names=1,check.names=FALSE)),"Type")
stopifnot(length(SELECTED_GENES) == 9L,
          length(ML_CANDIDATE_UNIVERSE) == 14L,
          all(SELECTED_GENES %in% ML_CANDIDATE_UNIVERSE))

for (pkg in c("ggplot2", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing package ", pkg, ". Install it: install.packages('", pkg, "')")
  }
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ------------------------- 1. STRICT FILE CHECKS ----------------------------
read_csv <- function(file) {
  utils::read.csv(file, check.names = FALSE, fileEncoding = "UTF-8-BOM",
                  stringsAsFactors = FALSE)
}
check_columns <- function(d, need, file) {
  missing <- setdiff(need, colnames(d))
  if (length(missing)) {
    stop("Missing column(s) in ", basename(file), ": ",
         paste(missing, collapse = ", "),
         "\nFound columns: ", paste(names(d), collapse = ", "))
  }
}
score <- read_csv(SCORE_FILE)
long <- read_csv(LONG_FILE)
check_columns(score, c("Sample", "Donor", "Group", "Score9"), SCORE_FILE)
check_columns(long, c("Sample", "Donor", "Group", "Gene", "Expression"), LONG_FILE)
score <- score[, c("Sample", "Donor", "Group", "Score9"), drop = FALSE]
long <- long[, c("Sample", "Donor", "Group", "Gene", "Expression"), drop = FALSE]

for (nm in c("Sample", "Donor", "Group")) {
  score[[nm]] <- trimws(as.character(score[[nm]]))
  long[[nm]] <- trimws(as.character(long[[nm]]))
}
score$Group <- toupper(score$Group)
long$Group <- toupper(long$Group)
long$Gene <- trimws(as.character(long$Gene))
score$Score9 <- suppressWarnings(as.numeric(score$Score9))
long$Expression <- suppressWarnings(as.numeric(long$Expression))

if (nrow(score) != 106L || length(unique(score$Donor)) != 53L) {
  stop("Training score table must have exactly 106 samples / 53 donors. Found: ",
       nrow(score), " rows and ", length(unique(score$Donor)), " donors")
}
if (anyNA(score) || any(!nzchar(score$Sample)) || any(!nzchar(score$Donor)) ||
    any(!is.finite(score$Score9)) || anyDuplicated(score$Sample)) {
  stop("Training score has missing/invalid or duplicate sample values")
}
if (!setequal(score$Group, c("CIT", "R2H"))) stop("Score groups must be CIT and R2H")
pair_count <- table(score$Donor, factor(score$Group, levels = c("CIT", "R2H")))
if (any(pair_count != 1L)) stop("Some training donors lack one CIT and one R2H sample")

if (nrow(long) != 106L * 9L || anyNA(long) ||
    any(!is.finite(long$Expression)) ||
    !setequal(unique(long$Gene), SELECTED_GENES)) {
  stop("SHAP long table must have 9 final features x 106 complete samples.\n",
       "Genes found: ", paste(unique(long$Gene), collapse = ", "))
}
key_long <- paste(long$Sample, long$Gene, sep = "::")
if (anyDuplicated(key_long)) stop("Repeated sample/gene in SHAP long table")
long_uniq <- unique(long[, c("Sample", "Donor", "Group")])
score_uniq <- score[, c("Sample", "Donor", "Group")]
row_key <- function(d) paste(d$Sample, d$Donor, d$Group, sep = "::")
if (nrow(long_uniq) != 106L || !setequal(row_key(long_uniq), row_key(score_uniq))) {
  stop("Training score CSV and SHAP long CSV sample/group/donor identities disagree")
}
if (!all(table(long$Sample) == 9L)) {
  stop("Every training sample needs exactly nine gene-expression records")
}
message("Confirmed: 53 paired donors / 106 training samples / 9 selected genes.")
message("Using existing score CSV, NO restandardization: ", SCORE_FILE)
message("Using expression values from: ", LONG_FILE)

# ------------------------ 2. COMMON PLOT THEME -------------------------------
plot_theme <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size, base_family = BASE_FONT) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black", face = "bold"),
      plot.title = ggplot2::element_text(
        color = "black", face = "bold", size = base_size + 1.6,
        hjust = 0.5, margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(
        color = "black", face = "bold", size = base_size - 0.4,
        hjust = 0.5, margin = ggplot2::margin(b = 7)),
      plot.caption = ggplot2::element_text(
        color = "black", face = "bold", size = base_size - 3),
      plot.tag = ggplot2::element_text(
        color = "black", face = "bold", size = base_size + 6),
      axis.title = ggplot2::element_text(
        color = "black", face = "bold", size = base_size),
      axis.text = ggplot2::element_text(
        color = "black", face = "bold", size = base_size - 1),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.70),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.56),
      strip.text = ggplot2::element_text(color = "black", face = "bold"),
      strip.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(color = "black", face = "bold"),
      legend.text = ggplot2::element_text(color = "black", face = "bold"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(9, 10, 9, 10))
}
format_p <- function(p) {
  if (!is.finite(p)) return("P = NA")
  if (p < 0.001) return("P < 0.001")
  sprintf("P = %.3f", p)
}

# --------------------- 3. PANEL F: PAIRED NINE-GENE SCORE -------------------
donors <- sort(unique(score$Donor))
ci <- match(donors, score$Donor[score$Group == "CIT"])
ri <- match(donors, score$Donor[score$Group == "R2H"])
cit <- score[score$Group == "CIT", , drop = FALSE]
r2h <- score[score$Group == "R2H", , drop = FALSE]
paired <- data.frame(
  Donor = donors,
  CIT = cit$Score9[ci],
  R2H = r2h$Score9[ri],
  stringsAsFactors = FALSE
)
paired$Delta_R2H_minus_CIT <- paired$R2H - paired$CIT
stopifnot(!anyNA(paired))
paired_test <- suppressWarnings(stats::wilcox.test(
  paired$R2H, paired$CIT, paired = TRUE, exact = FALSE))

score$Condition <- factor(score$Group, levels = c("CIT", "R2H"))
score$X <- ifelse(score$Group == "CIT", 1, 2)
offsets <- setNames(seq(-0.11, 0.11, length.out = length(donors)), donors)
score$PointX <- score$X + offsets[score$Donor]
lims <- range(score$Score9)
span <- max(diff(lims), 1e-3)
bracket_y <- lims[2] + 0.07 * span
text_y <- lims[2] + 0.17 * span
cols <- c(CIT = COL_CIT, R2H = COL_R2H)

pF <- ggplot2::ggplot(score, ggplot2::aes(x = X, y = Score9)) +
  ggplot2::geom_violin(
    ggplot2::aes(group = Condition, fill = Condition),
    width = 0.82, alpha = 0.52, trim = TRUE, color = NA,
    show.legend = FALSE) +
  ggplot2::geom_line(
    ggplot2::aes(x = PointX, group = Donor),
    color = "#7C838C", alpha = 0.49, linewidth = 0.45) +
  ggplot2::geom_boxplot(
    ggplot2::aes(group = Condition, fill = Condition, color = Condition),
    width = 0.27, alpha = 0.91, outlier.shape = NA,
    linewidth = 0.88, show.legend = FALSE) +
  ggplot2::geom_point(
    ggplot2::aes(x = PointX, color = Condition),
    size = 1.85, alpha = 0.96, show.legend = FALSE) +
  ggplot2::annotate("segment", x = 1, xend = 2,
                    y = bracket_y, yend = bracket_y,
                    color = "black", linewidth = 0.65) +
  ggplot2::annotate("segment", x = 1, xend = 1,
                    y = bracket_y - span * 0.025, yend = bracket_y,
                    color = "black", linewidth = 0.65) +
  ggplot2::annotate("segment", x = 2, xend = 2,
                    y = bracket_y - span * 0.025, yend = bracket_y,
                    color = "black", linewidth = 0.65) +
  ggplot2::annotate("text", x = 1.5, y = text_y,
                    label = paste0("Paired Wilcoxon  ", format_p(paired_test$p.value)),
                    color = "black", fontface = "bold", size = 3.5) +
  ggplot2::scale_x_discrete(labels=c(CIT="CIT",R2H="HR")) +
  ggplot2::scale_fill_manual(values = cols, labels=c(CIT="CIT",R2H="HR")) +
  ggplot2::scale_color_manual(values = cols, labels=c(CIT="CIT",R2H="HR")) +
  ggplot2::scale_x_continuous(breaks = 1:2, labels = c("CIT", "R2H"),
                              limits = c(0.47, 2.53)) +
  ggplot2::coord_cartesian(
    ylim = c(lims[1] - span * 0.09, lims[2] + span * 0.27), clip = "off") +
  ggplot2::labs(
    title = "RF-selected nine-gene coordinated score",
    subtitle = "GSE145989 training  |  53 matched donors",
    x = NULL, y = "Mean gene-wise training z-score",
    caption = "Grey lines connect CIT and R2H from the same donor") +
  plot_theme(10.7) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_line(
    color = "#EDEDED", linewidth = 0.30))

utils::write.csv(paired,
  file.path(OUTPUT_DIR, "Panel_F_53paired_9gene_score_and_delta.csv"),
  row.names = FALSE)
utils::write.csv(data.frame(
  Cohort = "GSE145989", Donors = 53L, Score = "Mean nine-gene training z-score",
  Test = "Paired two-sided Wilcoxon signed-rank",
  P = unname(paired_test$p.value),
  Median_Paired_Delta = stats::median(paired$Delta_R2H_minus_CIT)),
  file.path(OUTPUT_DIR, "Panel_F_paired_Wilcoxon_statistics.csv"), row.names = FALSE)

# ----------------- 4. PANEL G: EXPLORATORY GO / REACTOME ORA -----------------
# Importantly, background is the 22 pre-selection ML candidate genes rather
# than whole genome. G is exploratory, not independent pathway validation.
get_msig <- function(collection) {
  fm <- names(formals(msigdbr::msigdbr))
  collection_arg <- if ("collection" %in% fm) "collection" else "category"
  args <- list(species = "Homo sapiens")
  args[[collection_arg]] <- collection
  do.call(msigdbr::msigdbr, args)
}
ora_one <- function(msig, prefix, label, selected, universe) {
  required <- c("gs_name", "gene_symbol")
  if (!all(required %in% names(msig))) stop("Unexpected MSigDB column names")
  z <- unique(data.frame(
    term = as.character(msig$gs_name),
    gene = as.character(msig$gene_symbol), stringsAsFactors = FALSE))
  z <- z[grepl(prefix, z$term) & !is.na(z$gene) &
           nzchar(z$gene), , drop = FALSE]
  if (!nrow(z)) return(data.frame())
  annotated <- intersect(universe, unique(z$gene))
  selected <- intersect(selected, annotated)
  M <- length(annotated)
  N <- length(selected)
  if (M < 3L || N < 2L) return(data.frame())
  memberships <- split(z$gene, z$term)
  res <- lapply(names(memberships), function(term) {
    members <- intersect(unique(memberships[[term]]), annotated)
    K <- length(members)
    if (K < 2L || K >= M) return(NULL)
    k <- length(intersect(selected, members))
    data.frame(Term = term, Collection = label, Count = k,
               Background_Size = M, Query_Size = N,
               Set_In_Background = K,
               Fold = (k / N) / (K / M),
               P = stats::phyper(k - 1L, K, M - K, N, lower.tail = FALSE),
               stringsAsFactors = FALSE)
  })
  res <- Filter(Negate(is.null), res)
  if (!length(res)) return(data.frame())
  res <- do.call(rbind, res)
  # BH FDR over ALL tested terms in this collection, including zero overlaps.
  res$FDR <- stats::p.adjust(res$P, method = "BH")
  res[order(res$FDR, -res$Count, res$P), , drop = FALSE]
}
clean_term <- function(x, width = 27L) {
  x <- sub("^(GOBP_|REACTOME_)", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1))
}
ora_errors <- character(0)
go_res <- react_res <- data.frame()
if (requireNamespace("msigdbr", quietly = TRUE)) {
  go_ms <- tryCatch(get_msig("C5"), error = function(e) {
    ora_errors <<- c(ora_errors, paste("GO:", conditionMessage(e)))
    NULL
  })
  react_ms <- tryCatch(get_msig("C2"), error = function(e) {
    ora_errors <<- c(ora_errors, paste("Reactome:", conditionMessage(e)))
    NULL
  })
  if (!is.null(go_ms)) go_res <- ora_one(
    go_ms, "^GOBP_", "GO Biological Process",
    SELECTED_GENES, ML_CANDIDATE_UNIVERSE)
  if (!is.null(react_ms)) react_res <- ora_one(
    react_ms, "^REACTOME_", "Reactome",
    SELECTED_GENES, ML_CANDIDATE_UNIVERSE)
} else {
  ora_errors <- "Package msigdbr not installed"
}
if (nrow(go_res)) utils::write.csv(go_res,
  file.path(OUTPUT_DIR, "Panel_G_GO_BP_ORA_all_tested.csv"), row.names = FALSE)
if (nrow(react_res)) utils::write.csv(react_res,
  file.path(OUTPUT_DIR, "Panel_G_Reactome_ORA_all_tested.csv"), row.names = FALSE)

pick_top <- function(d, n) {
  if (!nrow(d)) return(d)
  d <- d[d$Count >= 2L, , drop = FALSE]
  head(d[order(d$FDR, -d$Count, d$P), , drop = FALSE], n)
}
chosen <- list(pick_top(go_res, N_SHOW_GO),
               pick_top(react_res, N_SHOW_REACTOME))
chosen <- chosen[vapply(chosen, nrow, integer(1)) > 0L]
ora_plot <- if (length(chosen)) do.call(rbind, chosen) else data.frame()
if (nrow(ora_plot)) {
  ora_plot$TermKey <- paste(ora_plot$Collection, ora_plot$Term, sep = "::")
  ora_plot$TermKey <- factor(ora_plot$TermKey,
                             levels = rev(unique(ora_plot$TermKey)))
  label_lookup <- stats::setNames(clean_term(ora_plot$Term),
                                   as.character(ora_plot$TermKey))
  ora_plot$NegLogFDR <- -log10(pmax(ora_plot$FDR, 1e-300))
  pG <- ggplot2::ggplot(ora_plot,
                        ggplot2::aes(x = Fold, y = TermKey)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = Fold, yend = TermKey),
      color = "#A7B1BB", linewidth = 0.82) +
    ggplot2::geom_point(
      ggplot2::aes(size = Count, color = NegLogFDR), alpha = 1.0) +
    ggplot2::facet_grid(Collection ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_y_discrete(labels = function(x) unname(label_lookup[x])) +
    ggplot2::scale_color_gradient(
      low = "#246FD2", high = "#E3292B", name = "-log10(FDR)") +
    ggplot2::scale_size_continuous(range = c(3.0, 6.8), name = "Overlap") +
    ggplot2::labs(
      title = "Functional convergence of nine genes",
      subtitle = "Exploratory GO BP / Reactome  |  14-gene background",
      x = "Fold enrichment", y = NULL,
      caption = "BH-FDR adjusted across all tested terms per collection") +
    plot_theme(10.1) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        color = "black", face = "bold", size = 8.2),
      strip.text.y = ggplot2::element_text(
        color = "black", face = "bold", angle = 0, size = 9.5),
      panel.grid.major.x = ggplot2::element_line(
        color = "#EEEEEE", linewidth = 0.32),
      legend.position = "right",
      legend.key.height = grid::unit(0.30, "cm"))
} else {
  reason <- if (length(ora_errors)) {
    "MSigDB unavailable\nSee R messages / install msigdbr"
  } else {
    "No GO / Reactome terms contain\nat least two of the nine genes\nunder the 14-gene background"
  }
  pG <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = reason,
                      color = "black", fontface = "bold", size = 3.7) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::labs(
      title = "Functional convergence of nine genes",
      subtitle = "Exploratory GO BP / Reactome  |  14-gene background") +
    ggplot2::theme_void(base_family = BASE_FONT) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = "black", face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = ggplot2::element_text(
        color = "black", face = "bold", hjust = 0.5, size = 10))
  if (length(ora_errors)) warning(paste(ora_errors, collapse = " | "))
}

# ---------------- 5. PANEL H: PAIRED EXPRESSION-DELTA SPEARMAN ---------------
# One expression change PER DONOR PER GENE: Delta = R2H - CIT.
# Long CSV values are the ORIGINAL expression, not SHAP contributions.
expr_changes <- matrix(NA_real_, nrow = length(donors),
                       ncol = length(SELECTED_GENES),
                       dimnames = list(donors, SELECTED_GENES))
for (j in seq_along(SELECTED_GENES)) {
  g <- SELECTED_GENES[j]
  dat <- long[long$Gene == g, , drop = FALSE]
  ca <- dat[dat$Group == "CIT", , drop = FALSE]
  ra <- dat[dat$Group == "R2H", , drop = FALSE]
  if (length(unique(dat$Donor)) != 53L ||
      anyDuplicated(ca$Donor) || anyDuplicated(ra$Donor) ||
      !setequal(ca$Donor, donors) || !setequal(ra$Donor, donors)) {
    stop("Mismatch when pairing gene ", g, " by donor")
  }
  expr_changes[, j] <- ra$Expression[match(donors, ra$Donor)] -
                       ca$Expression[match(donors, ca$Donor)]
}
if (anyNA(expr_changes)) stop("NA found in donor-level gene differences")
utils::write.csv(
  data.frame(Donor = donors, expr_changes, check.names = FALSE),
  file.path(OUTPUT_DIR, "Panel_H_53donors_9gene_R2H_minus_CIT_deltas.csv"),
  row.names = FALSE)

ngenes <- length(SELECTED_GENES)
rho <- matrix(NA_real_, ngenes, ngenes,
              dimnames = list(SELECTED_GENES, SELECTED_GENES))
p_mat <- rho
fdr_mat <- rho
diag(rho) <- 1
for (i in seq_len(ngenes - 1L)) {
  for (j in seq.int(i + 1L, ngenes)) {
    ct <- tryCatch(suppressWarnings(
      stats::cor.test(expr_changes[, i], expr_changes[, j],
                      method = "spearman", exact = FALSE)),
      error = function(e) NULL)
    if (!is.null(ct)) {
      rho[i, j] <- rho[j, i] <- unname(ct$estimate)
      p_mat[i, j] <- p_mat[j, i] <- ct$p.value
    }
  }
}
upper <- which(upper.tri(p_mat) & is.finite(p_mat), arr.ind = TRUE)
if (nrow(upper)) {
  fdr <- stats::p.adjust(p_mat[upper], method = "BH")
  for (k in seq_len(nrow(upper))) {
    i <- upper[k, 1]; j <- upper[k, 2]
    fdr_mat[i, j] <- fdr_mat[j, i] <- fdr[k]
  }
}
utils::write.csv(rho,
  file.path(OUTPUT_DIR, "Panel_H_9gene_delta_Spearman_rho.csv"))
utils::write.csv(p_mat,
  file.path(OUTPUT_DIR, "Panel_H_9gene_delta_Spearman_P.csv"))
utils::write.csv(fdr_mat,
  file.path(OUTPUT_DIR, "Panel_H_9gene_delta_Spearman_FDR_BH.csv"))

plot_order <- c("SERPINE1", setdiff(SELECTED_GENES, "SERPINE1"))
corrgrid <- expand.grid(GeneX = plot_order, GeneY = plot_order,
                        stringsAsFactors = FALSE)
corrgrid$Rho <- mapply(function(x, y) rho[y, x],
                       corrgrid$GeneX, corrgrid$GeneY)
corrgrid$FDR <- mapply(function(x, y) fdr_mat[y, x],
                       corrgrid$GeneX, corrgrid$GeneY)
corrgrid$Label <- ifelse(is.finite(corrgrid$Rho),
                          sprintf("%.2f", corrgrid$Rho), "NA")
corrgrid$Stars <- ifelse(is.na(corrgrid$FDR), "",
  ifelse(corrgrid$FDR < 0.001, "***",
    ifelse(corrgrid$FDR < 0.01, "**",
      ifelse(corrgrid$FDR < 0.05, "*", ""))))
corrgrid$Label <- paste0(corrgrid$Label, corrgrid$Stars)
corrgrid$GeneX <- factor(corrgrid$GeneX, levels = plot_order)
corrgrid$GeneY <- factor(corrgrid$GeneY, levels = rev(plot_order))

pH <- ggplot2::ggplot(corrgrid,
                       ggplot2::aes(x = GeneX, y = GeneY, fill = Rho)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.90) +
  ggplot2::geom_text(
    ggplot2::aes(label = Label), color = "black",
    fontface = "bold", size = 2.7) +
  ggplot2::scale_fill_gradient2(
    low = "#2476D9", mid = "#FFFFFF", high = "#EB2E30",
    midpoint = 0, limits = c(-1, 1), na.value = "#F1F1F1",
    name = "Spearman rho") +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Paired-change co-expression",
    subtitle = "GSE145989  |  53 donor-level R2H - CIT changes",
    x = NULL, y = NULL,
    caption = "Stars: BH-FDR <0.05 / <0.01 / <0.001 (36 tests)") +
  plot_theme(9.7) +
  ggplot2::theme(
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, color = "black", face = "bold", size = 8),
    axis.text.y = ggplot2::element_text(
      color = "black", face = "bold", size = 8),
    legend.position = "right",
    legend.key.height = grid::unit(0.33, "cm"))

# --------------------------- 6. EXPORT --------------------------------------
final_figure <- (pF | pG | pH) +
  patchwork::plot_layout(widths = c(1.10, 1.20, 1.08)) +
  patchwork::plot_annotation(
    tag_levels = list(c("F", "G", "H")),
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        color = "black", face = "bold", size = 19),
      plot.tag.position = c(0, 1),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)))

save_formats <- function(plot, base, width, height) {
  pdf_file <- paste0(base, ".pdf")
  tryCatch(
    ggplot2::ggsave(pdf_file, plot, width = width, height = height,
                    units = "in", device = grDevices::cairo_pdf, bg = "white"),
    error = function(e) {
      message("Cairo PDF unavailable; falling back to base PDF: ",
              conditionMessage(e))
      ggplot2::ggsave(pdf_file, plot, width = width, height = height,
                      units = "in", device = "pdf", useDingbats = FALSE,
                      bg = "white")
    })
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(paste0(base, ".png"), plot,
                    width = width, height = height, units = "in",
                    dpi = 600, device = ragg::agg_png, bg = "white")
    ggplot2::ggsave(paste0(base, ".tiff"), plot,
                    width = width, height = height, units = "in",
                    dpi = 600, device = ragg::agg_tiff, bg = "white")
  } else {
    ggplot2::ggsave(paste0(base, ".png"), plot,
                    width = width, height = height, units = "in",
                    dpi = 600, bg = "white")
    ggplot2::ggsave(paste0(base, ".tiff"), plot,
                    width = width, height = height, units = "in",
                    dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  }
  message("Saved: ", base, " [pdf / png / tiff]")
}

save_formats(final_figure,
  file.path(OUTPUT_DIR, "Figure_FGH_GSE145989_RF9_TRAIN_53pairs"),
  EXPORT_WIDTH, EXPORT_HEIGHT)
INDIV_DIR <- file.path(OUTPUT_DIR, "Individual_Panels")
dir.create(INDIV_DIR, recursive = TRUE, showWarnings = FALSE)
panels <- list(
  Panel_F_Training_9Gene_Paired_Violin = pF + ggplot2::labs(tag = "F"),
  Panel_G_Training_9Gene_Exploratory_ORA = pG + ggplot2::labs(tag = "G"),
  Panel_H_Training_9Gene_PairedDelta_Spearman = pH + ggplot2::labs(tag = "H")
)
sizes <- list(c(6.7, 5.9), c(8.6, 7.0), c(8.2, 7.1))
for (i in seq_along(panels)) {
  p <- panels[[i]] + ggplot2::theme(
    plot.tag = ggplot2::element_text(color = "black", face = "bold", size = 18),
    plot.tag.position = c(0, 1))
  save_formats(p, file.path(INDIV_DIR, names(panels)[i]),
               sizes[[i]][1], sizes[[i]][2])
}
print(final_figure)
message("DONE. Input: ", TRAIN_SHAP_DIR)
message("Output: ", normalizePath(OUTPUT_DIR, winslash = "/"))
message("Score paired Wilcoxon: ", format_p(paired_test$p.value))
message("9 genes: ", paste(SELECTED_GENES, collapse = ", "))

}

# =============================================================================
# STAGE figure3_upstream_qc
# =============================================================================
PIPELINE$figure3_upstream_qc <- function() {
# scripts/04_single_cell_qc.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

## 0. Load packages
## =========================
library(Seurat)
library(DoubletFinder)
library(harmony)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)

set.seed(1234)

dir.create("01_QC_DoubletFinder", showWarnings = FALSE)
dir.create("02_Global_Clustering", showWarnings = FALSE)
dir.create("03_FindAllMarkers", showWarnings = FALSE)
dir.create("03_FindAllMarkers/markers_by_cluster", showWarnings = FALSE)


## =========================
## 1. Read 10X data
## barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz
## =========================
HR1.data <- Read10X(file.path(repo_input("tenx_directory"),"HR1"))
HR2.data <- Read10X(file.path(repo_input("tenx_directory"),"HR2"))
HR3.data <- Read10X(file.path(repo_input("tenx_directory"),"HR3"))
HR4.data <- Read10X(file.path(repo_input("tenx_directory"),"HR4"))
HR5.data <- Read10X(file.path(repo_input("tenx_directory"),"HR5"))
HR6.data <- Read10X(file.path(repo_input("tenx_directory"),"HR6"))

CIT1.data <- Read10X(file.path(repo_input("tenx_directory"),"CIT1"))
CIT2.data <- Read10X(file.path(repo_input("tenx_directory"),"CIT2"))
CIT3.data <- Read10X(file.path(repo_input("tenx_directory"),"CIT3"))
CIT4.data <- Read10X(file.path(repo_input("tenx_directory"),"CIT4"))
CIT5.data <- Read10X(file.path(repo_input("tenx_directory"),"CIT5"))
CIT6.data <- Read10X(file.path(repo_input("tenx_directory"),"CIT6"))


## =========================
## 2. Create Seurat objects
## =========================
HR1 <- CreateSeuratObject(HR1.data, project = "HR1")
HR2 <- CreateSeuratObject(HR2.data, project = "HR2")
HR3 <- CreateSeuratObject(HR3.data, project = "HR3")
HR4 <- CreateSeuratObject(HR4.data, project = "HR4")
HR5 <- CreateSeuratObject(HR5.data, project = "HR5")
HR6 <- CreateSeuratObject(HR6.data, project = "HR6")

CIT1 <- CreateSeuratObject(CIT1.data, project = "CIT1")
CIT2 <- CreateSeuratObject(CIT2.data, project = "CIT2")
CIT3 <- CreateSeuratObject(CIT3.data, project = "CIT3")
CIT4 <- CreateSeuratObject(CIT4.data, project = "CIT4")
CIT5 <- CreateSeuratObject(CIT5.data, project = "CIT5")
CIT6 <- CreateSeuratObject(CIT6.data, project = "CIT6")

rm(
  HR1.data, HR2.data, HR3.data, HR4.data, HR5.data, HR6.data,
  CIT1.data, CIT2.data, CIT3.data, CIT4.data, CIT5.data, CIT6.data
)
gc()


## =========================
## 3. Define helper genes
## =========================
HB.genes <- c(
  "HBA1", "HBA2", "HBB",
  "HBD", "HBG1", "HBG2"
)


## =========================
## 4. DoubletFinder version compatibility
## =========================
run_param_sweep <- function(seu, pcs, sct = FALSE) {
  if (exists("paramSweep_v3", mode = "function")) {
    return(paramSweep_v3(seu, PCs = pcs, sct = sct))
  }
  if (exists("paramSweep", mode = "function")) {
    return(paramSweep(seu, PCs = pcs, sct = sct))
  }
  stop(" paramSweep_v3() paramSweep().")
}

run_doublet_finder <- function(
    seu, pcs, pN, pK, nExp,
    reuse.pANN = FALSE, sct = FALSE
) {
  if (exists("doubletFinder_v3", mode = "function")) {
    return(
      doubletFinder_v3(
        seu,
        PCs = pcs,
        pN = pN,
        pK = pK,
        nExp = nExp,
        reuse.pANN = reuse.pANN,
        sct = sct
      )
    )
  }
  if (exists("doubletFinder", mode = "function")) {
    return(
      doubletFinder(
        seu,
        PCs = pcs,
        pN = pN,
        pK = pK,
        nExp = nExp,
        reuse.pANN = reuse.pANN,
        sct = sct
      )
    )
  }
  stop(" doubletFinder_v3() doubletFinder().")
}


## =========================
## 5. Single-sample processing
## QC + DoubletFinder
## =========================
process_sample_with_DoubletFinder <- function(
    seu,
    sample_name,
    minGene = 300,
    maxGene = 7500,
    pctmt = 15,
    pctHB = 5,
    pcs = 1:15,
    resolution = 0.2,
    pN = 0.25
) {
  cat("\n============================\n")
  cat("Processing sample:", sample_name, "\n")
  cat("============================\n")

  # Source comment omitted during English-language consolidation.
  seu$sample_id <- sample_name
  seu$condition <- stringr::str_remove(sample_name, "[0-9]+$")

  ## QC metrics
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  hb_use <- intersect(HB.genes, rownames(seu))
  if (length(hb_use) > 0) {
    seu[["percent.HB"]] <- PercentageFeatureSet(seu, features = hb_use)
  } else {
    seu[["percent.HB"]] <- 0
  }

  ## QC plots before filtering
  p_qc_before <- VlnPlot(
    seu,
    features = c(
      "nFeature_RNA", "nCount_RNA",
      "percent.mt", "percent.HB"
    ),
    ncol = 4,
    pt.size = 0.01
  ) + plot_annotation(title = paste0(sample_name, " - before QC"))

  ggsave(
    filename = file.path(
      "01_QC_DoubletFinder",
      paste0(sample_name, "_QC_before.pdf")
    ),
    plot = p_qc_before,
    width = 12,
    height = 4
  )

  ## QC filtering
  before_n <- ncol(seu)

  seu <- subset(
    seu,
    subset =
      nFeature_RNA > minGene &
      nFeature_RNA < maxGene &
      percent.mt < pctmt &
      percent.HB < pctHB
  )

  after_qc_n <- ncol(seu)

  cat("Cells before QC:", before_n, "\n")
  cat("Cells after QC :", after_qc_n, "\n")

  if (after_qc_n < 100) {
    stop(
      "Check source data or configuration. ", sample_name,
      "Check source data or configuration. "
    )
  }

  ## QC plots after filtering
  p_qc_after <- VlnPlot(
    seu,
    features = c(
      "nFeature_RNA", "nCount_RNA",
      "percent.mt", "percent.HB"
    ),
    ncol = 4,
    pt.size = 0.01
  ) + plot_annotation(title = paste0(sample_name, " - after QC"))

  ggsave(
    filename = file.path(
      "01_QC_DoubletFinder",
      paste0(sample_name, "_QC_after.pdf")
    ),
    plot = p_qc_after,
    width = 12,
    height = 4
  )

  ## Preprocess for DoubletFinder
  seu <- NormalizeData(seu, verbose = FALSE) %>%
    FindVariableFeatures(
      selection.method = "vst",
      nfeatures = 2000,
      verbose = FALSE
    ) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE)

  seu <- FindNeighbors(
    seu,
    dims = pcs,
    verbose = FALSE
  ) %>%
    FindClusters(
      resolution = resolution,
      verbose = FALSE
    ) %>%
    RunUMAP(
      dims = pcs,
      verbose = FALSE
    )

  ## Parameter sweep
  set.seed(1234)

  sweep.res.list <- run_param_sweep(
    seu = seu,
    pcs = pcs,
    sct = FALSE
  )

  sweep.stats <- summarizeSweep(
    sweep.res.list,
    GT = FALSE
  )

  bcmvn <- find.pK(sweep.stats)

  write.csv(
    bcmvn,
    file = file.path(
      "01_QC_DoubletFinder",
      paste0(sample_name, "_DoubletFinder_pK.csv")
    ),
    row.names = FALSE
  )

  if (nrow(bcmvn) == 0 || all(is.na(bcmvn$BCmetric))) {
    stop("Check source data or configuration. ", sample_name, "  BCmetric.")
  }

  bcmetric_num <- as.numeric(
    as.character(bcmvn$BCmetric)
  )

  best_row <- which.max(bcmetric_num)

  best_pK <- as.numeric(
    as.character(bcmvn$pK[best_row])
  )

  if (length(best_pK) == 0 || is.na(best_pK)) {
    stop("Check source data or configuration. ", sample_name, "Check source data or configuration. ")
  }

  cat("Best pK =", best_pK, "\n")

  ## Estimate expected doublets
  homotypic.prop <- modelHomotypic(
    seu$seurat_clusters
  )

  # Source comment omitted during English-language consolidation.
  doublet_rate <- min(
    0.008 * ncol(seu) / 1000,
    0.10
  )

  nExp_poi <- max(
    round(doublet_rate * ncol(seu)),
    1
  )

  nExp_poi.adj <- max(
    round(nExp_poi * (1 - homotypic.prop)),
    1
  )

  cat("Recovered cells:", ncol(seu), "\n")
  cat(
    "Estimated doublet rate:",
    round(doublet_rate * 100, 2),
    "%\n"
  )
  cat("Homotypic proportion:", homotypic.prop, "\n")
  cat("Expected doublets (raw):", nExp_poi, "\n")
  cat("Expected doublets (adjusted):", nExp_poi.adj, "\n")

  df_parameters <- data.frame(
    sample_id = sample_name,
    cells_after_QC = ncol(seu),
    best_pK = best_pK,
    pN = pN,
    doublet_rate = doublet_rate,
    homotypic_prop = homotypic.prop,
    nExp_raw = nExp_poi,
    nExp_adjusted = nExp_poi.adj
  )

  write.csv(
    df_parameters,
    file = file.path(
      "01_QC_DoubletFinder",
      paste0(sample_name, "_DoubletFinder_parameters.csv")
    ),
    row.names = FALSE
  )

  ## Run DoubletFinder
  set.seed(1234)

  seu <- run_doublet_finder(
    seu = seu,
    pcs = pcs,
    pN = pN,
    pK = best_pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  ## Get classification column
  df_cols <- grep(
    "^DF.classifications",
    colnames(seu@meta.data),
    value = TRUE
  )

  if (length(df_cols) == 0) {
    stop(
      "DoubletFinder ",
      sample_name
    )
  }

  df_col <- tail(df_cols, 1)
  seu$doublet_status <- seu@meta.data[[df_col]]

  cat("Doublet classification table:\n")
  print(table(seu$doublet_status))

  ## Save summary
  df_summary <- as.data.frame(
    table(seu$doublet_status)
  )
  colnames(df_summary) <- c(
    "doublet_status",
    "cell_number"
  )
  df_summary$sample_id <- sample_name

  write.csv(
    df_summary,
    file = file.path(
      "01_QC_DoubletFinder",
      paste0(sample_name, "_DoubletFinder_summary.csv")
    ),
    row.names = FALSE
  )

  ## Doublet UMAP
  p_doublet <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "doublet_status"
  ) +
    ggtitle(paste0(sample_name, " - DoubletFinder"))

  ggsave(
    filename = file.path(
      "01_QC_DoubletFinder",
      paste0(sample_name, "_DoubletFinder_UMAP.pdf")
    ),
    plot = p_doublet,
    width = 6,
    height = 5
  )

  ## Keep singlets
  seu_clean <- subset(
    seu,
    subset = doublet_status == "Singlet"
  )

  after_doublet_n <- ncol(seu_clean)

  cat(
    "Cells after removing doublets:",
    after_doublet_n,
    "\n"
  )

  filtering_summary <- data.frame(
    sample_id = sample_name,
    cells_before_QC = before_n,
    cells_after_QC = after_qc_n,
    cells_after_DoubletFinder = after_doublet_n
  )

  write.csv(
    filtering_summary,
    file = file.path(
      "01_QC_DoubletFinder",
      paste0(sample_name, "_cell_filtering_summary.csv")
    ),
    row.names = FALSE
  )

  saveRDS(
    seu_clean,
    file = file.path(
      "01_QC_DoubletFinder",
      paste0(sample_name, "_clean.rds")
    )
  )

  return(seu_clean)
}


## =========================
## 6. Process each sample
## =========================
HR1 <- process_sample_with_DoubletFinder(HR1, "HR1")
HR2 <- process_sample_with_DoubletFinder(HR2, "HR2")
HR3 <- process_sample_with_DoubletFinder(HR3, "HR3")
HR4 <- process_sample_with_DoubletFinder(HR4, "HR4")
HR5 <- process_sample_with_DoubletFinder(HR5, "HR5")
HR6 <- process_sample_with_DoubletFinder(HR6, "HR6")

CIT1 <- process_sample_with_DoubletFinder(CIT1, "CIT1")
CIT2 <- process_sample_with_DoubletFinder(CIT2, "CIT2")
CIT3 <- process_sample_with_DoubletFinder(CIT3, "CIT3")
CIT4 <- process_sample_with_DoubletFinder(CIT4, "CIT4")
CIT5 <- process_sample_with_DoubletFinder(CIT5, "CIT5")
CIT6 <- process_sample_with_DoubletFinder(CIT6, "CIT6")

gc()


## =========================
## 7. Merge singlets
## =========================
scRNAlist <- merge(
  x = HR1,
  y = list(
    HR2, HR3, HR4, HR5, HR6,
    CIT1, CIT2, CIT3, CIT4, CIT5, CIT6
  ),
  add.cell.ids = c(
    "HR1", "HR2", "HR3", "HR4", "HR5", "HR6",
    "CIT1", "CIT2", "CIT3", "CIT4", "CIT5", "CIT6"
  ),
  project = "Lung_CIT_HR"
)

scRNAlist$tissue_type <- scRNAlist$condition

print(table(scRNAlist$orig.ident))
print(table(scRNAlist$sample_id))
print(table(scRNAlist$tissue_type))

write.csv(
  as.data.frame(
    table(
      sample_id = scRNAlist$sample_id,
      condition = scRNAlist$condition
    )
  ),
  file = file.path(
    "02_Global_Clustering",
    "Cell_numbers_after_merge.csv"
  ),
  row.names = FALSE
)


## =========================
## 8. Cell cycle scoring
## =========================
g2m_genes <- CaseMatch(
  search = cc.genes$g2m.genes,
  match = rownames(scRNAlist)
)

s_genes <- CaseMatch(
  search = cc.genes$s.genes,
  match = rownames(scRNAlist)
)

scRNAlist <- CellCycleScoring(
  object = scRNAlist,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = FALSE
)

p_cycle <- scRNAlist@meta.data %>%
  ggplot(
    aes(
      x = S.Score,
      y = G2M.Score,
      color = Phase
    )
  ) +
  geom_point(size = 0.2) +
  theme_minimal()

ggsave(
  filename = file.path(
    "02_Global_Clustering",
    "Cell_cycle_scores.pdf"
  ),
  plot = p_cycle,
  width = 6,
  height = 5
)


## =========================
## 9. Main preprocessing after merge
## =========================
DefaultAssay(scRNAlist) <- "RNA"

scRNAlist <- NormalizeData(scRNAlist) %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = 3000
  ) %>%
  ScaleData() %>%
  RunPCA(
    npcs = 30,
    verbose = TRUE
  )

## PCA by group / sample
p_pca_group <- DimPlot(
  scRNAlist,
  reduction = "pca",
  group.by = "tissue_type"
) +
  ggtitle("PCA by tissue_type")

p_pca_sample <- DimPlot(
  scRNAlist,
  reduction = "pca",
  group.by = "orig.ident"
) +
  ggtitle("PCA by orig.ident")

ggsave(
  filename = file.path(
    "02_Global_Clustering",
    "PCA_before_Harmony.pdf"
  ),
  plot = p_pca_group + p_pca_sample,
  width = 10,
  height = 8
)

## Variable genes visualization
top15 <- head(
  VariableFeatures(scRNAlist),
  15
)

plot1 <- VariableFeaturePlot(scRNAlist)

plot2 <- LabelPoints(
  plot = plot1,
  points = top15,
  repel = TRUE,
  size = 3
)

feat_15 <- plot1 + plot2

ggsave(
  filename = file.path(
    "02_Global_Clustering",
    "VariableFeatures_top15.pdf"
  ),
  plot = feat_15,
  width = 12,
  height = 6
)


## =========================
## 10. Harmony integration
## =========================
scRNA_harmony <- RunHarmony(
  object = scRNAlist,
  group.by.vars = "orig.ident",
  reduction = "pca"
)

p_harmony_group <- DimPlot(
  scRNA_harmony,
  reduction = "harmony",
  group.by = "tissue_type"
) +
  ggtitle("Harmony by tissue_type")

p_harmony_sample <- DimPlot(
  scRNA_harmony,
  reduction = "harmony",
  group.by = "orig.ident"
) +
  ggtitle("Harmony by orig.ident")

ggsave(
  filename = file.path(
    "02_Global_Clustering",
    "Harmony_embedding.pdf"
  ),
  plot = p_harmony_group + p_harmony_sample,
  width = 10,
  height = 8
)


## =========================
## 11. Clustering and dimensional reduction
## =========================
p_elbow <- ElbowPlot(
  scRNA_harmony,
  ndims = 30,
  reduction = "harmony"
)

ggsave(
  filename = file.path(
    "02_Global_Clustering",
    "ElbowPlot_harmony.pdf"
  ),
  plot = p_elbow,
  width = 6,
  height = 5
)

use_dims <- 1:20

scRNA_harmony <- FindNeighbors(
  scRNA_harmony,
  reduction = "harmony",
  dims = use_dims
) %>%
  FindClusters(
    resolution = 0.5
  )

scRNA_harmony <- RunTSNE(
  scRNA_harmony,
  reduction = "harmony",
  dims = use_dims
)

scRNA_harmony <- RunUMAP(
  scRNA_harmony,
  reduction = "harmony",
  dims = use_dims,
  n.neighbors = 10,
  min.dist = 0.5
)

## Plot clustering
tsne_integrated1 <- DimPlot(
  scRNA_harmony,
  reduction = "tsne",
  group.by = "tissue_type"
)

tsne_integrated2 <- DimPlot(
  scRNA_harmony,
  reduction = "tsne",
  label = TRUE
)

umap_integrated1 <- DimPlot(
  scRNA_harmony,
  reduction = "umap",
  group.by = "tissue_type"
)

umap_integrated2 <- DimPlot(
  scRNA_harmony,
  reduction = "umap",
  label = TRUE
)

umap_tsne_integrated <- (
  tsne_integrated1 + tsne_integrated2
) / (
  umap_integrated1 + umap_integrated2
)

ggsave(
  filename = file.path(
    "02_Global_Clustering",
    "UMAP_tSNE_integrated.pdf"
  ),
  plot = umap_tsne_integrated,
  width = 14,
  height = 10
)

## UMAP by sample
p_umap_sample <- DimPlot(
  scRNA_harmony,
  reduction = "umap",
  group.by = "orig.ident"
)

ggsave(
  filename = file.path(
    "02_Global_Clustering",
    "UMAP_by_sample.pdf"
  ),
  plot = p_umap_sample,
  width = 9,
  height = 7
)

## Split by group
p_umap_split <- DimPlot(
  scRNA_harmony,
  reduction = "umap",
  label = TRUE,
  split.by = "tissue_type"
)

ggsave(
  filename = file.path(
    "02_Global_Clustering",
    "UMAP_split_by_group.pdf"
  ),
  plot = p_umap_split,
  width = 12,
  height = 6
)

cluster_sample_table <- as.data.frame(
  table(
    cluster = scRNA_harmony$seurat_clusters,
    sample_id = scRNA_harmony$sample_id,
    condition = scRNA_harmony$condition
  )
)

write.csv(
  cluster_sample_table,
  file = file.path(
    "02_Global_Clustering",
    "Cluster_cell_numbers_by_sample.csv"
  ),
  row.names = FALSE
)


## =========================
## 12. Find cluster markers
## =========================
DefaultAssay(scRNA_harmony) <- "RNA"
Idents(scRNA_harmony) <- "seurat_clusters"

markers <- FindAllMarkers(
  object = scRNA_harmony,
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  return.thresh = 1
)

all.markers <- markers %>%
  dplyr::select(
    gene,
    everything()
  ) %>%
  dplyr::filter(
    p_val_adj < 0.05
  )

top10 <- all.markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

top20 <- all.markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  markers,
  file = file.path(
    "03_FindAllMarkers",
    "FindAllMarkers_all_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  all.markers,
  file = file.path(
    "03_FindAllMarkers",
    "FindAllMarkers_significant.csv"
  ),
  row.names = FALSE
)

write.csv(
  top10,
  file = file.path(
    "03_FindAllMarkers",
    "Top10_markers_each_cluster.csv"
  ),
  row.names = FALSE
)

write.csv(
  top20,
  file = file.path(
    "03_FindAllMarkers",
    "Top20_markers_each_cluster.csv"
  ),
  row.names = FALSE
)

cluster_ids <- sort(
  unique(all.markers$cluster)
)

for (cl in cluster_ids) {
  marker_cl <- all.markers %>%
    filter(cluster == cl) %>%
    arrange(desc(avg_log2FC))

  write.csv(
    marker_cl,
    file = file.path(
      "03_FindAllMarkers",
      "markers_by_cluster",
      paste0("Cluster_", cl, "_markers.csv")
    ),
    row.names = FALSE
  )
}


## =========================
## 13. Marker heatmap
## =========================
heatmap_genes <- unique(top10$gene)

heatmap_genes <- intersect(
  heatmap_genes,
  rownames(scRNA_harmony)
)

p_heat <- DoHeatmap(
  scRNA_harmony,
  features = heatmap_genes,
  assay = "RNA",
  slot = "scale.data",
  group.by = "seurat_clusters"
) +
  NoLegend()

ggsave(
  filename = file.path(
    "03_FindAllMarkers",
    "Top10_marker_heatmap.pdf"
  ),
  plot = p_heat,
  width = 12,
  height = 14,
  limitsize = FALSE
)


## =========================
## 14. Marker inspection for manual annotation
## =========================
genes_to_check <- c(
  ## B cell
  "EBF1", "CD79A", "CD79B", "MS4A1",
  "CD19", "CD37", "IGHM",

  ## T cell
  "CD3D", "CD3E", "TRBC1", "TRBC2",
  "IL7R", "TCF7", "CCR7",

  ## NK cell
  "NKG7", "GNLY", "GZMA", "GZMB",
  "PRF1", "KLRD1",

  ## Neutrophil
  "S100A8", "S100A9", "S100A12",
  "FCGR3B", "MMP8", "CSF3R",

  ## Macrophage / Monocyte
  "LYZ", "FCN1", "CCR2", "CSF1R",
  "MRC1", "MARCO", "C1QA", "C1QB",
  "C1QC", "CD68",

  ## Endothelial
  "PTPRB", "CDH5", "KDR", "CLDN5",
  "FLT1", "EPAS1", "PLVAP", "ESM1",
  "PECAM1", "VWF",

  ## Fibroblast
  "DCN", "COL1A1", "COL1A2",
  "COL3A1", "LUM", "PDGFRA",

  ## AT1
  "AGER", "CAV1", "CLIC5", "HOPX",

  ## AT2
  "SFTPA1", "SFTPA2", "SFTPB",
  "SFTPC", "NAPSA", "ABCA3",

  ## Other epithelial / Club
  "EPCAM", "KRT8", "KRT18", "KRT19",
  "SCGB1A1", "SCGB3A1",

  ## Ciliated
  "CAPS", "TPPP3", "FOXJ1", "PIFO",

  ## Mast
  "TPSAB1", "TPSB2", "KIT", "CPA3",

  ## Platelet
  "PF4", "PPBP",

  ## Cycling
  "MKI67", "TOP2A",

  ## Target gene
  "BAG3"
)

genes_to_check <- intersect(
  genes_to_check,
  rownames(scRNA_harmony)
)

p_dot <- DotPlot(
  scRNA_harmony,
  features = genes_to_check,
  assay = "RNA",
  group.by = "seurat_clusters"
) +
  coord_flip() +
  ggtitle("Marker inspection for manual annotation")

ggsave(
  filename = file.path(
    "03_FindAllMarkers",
    "Cluster_marker_dotplot.pdf"
  ),
  plot = p_dot,
  width = 10,
  height = 14,
  limitsize = FALSE
)


## =========================
## 15. Manual annotation template
## =========================
cluster_levels <- levels(
  Idents(scRNA_harmony)
)

if (is.null(cluster_levels)) {
  cluster_levels <- sort(
    unique(
      as.character(scRNA_harmony$seurat_clusters)
    )
  )
}

annotation_template <- data.frame(
  cluster = cluster_levels,
  celltype = "",
  comments = "",
  stringsAsFactors = FALSE
)

write.csv(
  annotation_template,
  file = file.path(
    "03_FindAllMarkers",
    "Manual_annotation_template.csv"
  ),
  row.names = FALSE
)


## =========================
## 16. Save object
## =========================
saveRDS(
  scRNA_harmony,
  file = file.path(
    "03_FindAllMarkers",
    "Seurat_object_after_FindAllMarkers.rds"
  )
)

save.image(
  file = file.path(
    "03_FindAllMarkers",
    "scRNA_workspace_after_FindAllMarkers.RData"
  )
)

cat("\n============================================\n")
cat(" FindAllMarkers.\n")
cat(" 03_FindAllMarkers marker \n")
cat("============================================\n")

saveRDS(scRNA_harmony,file.path(repo_output("Figure3/QC"),"global_unannotated.rds"))

}

# =============================================================================
# STAGE figure3_atlas
# =============================================================================
PIPELINE$figure3_atlas <- function() {
# scripts/05_atlas_panels.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.
sce <- load_sce()
# =============================================================================
#
#
#
#   hsa04218 Cellular senescence
#   hsa04066 HIF-1 signaling pathway
#   hsa04512 ECM-receptor interaction
#   hsa04520 Adherens junction
#
# =============================================================================

options(stringsAsFactors = FALSE)
set.seed(1234)

# =============================================================================
# =============================================================================

OUTDIR <- repo_output("Figure3/Atlas")
DIR_PANEL <- file.path(OUTDIR, "Individual_panels")
DIR_GSEA <- file.path(OUTDIR, "GSEA_individual")
DIR_TABLE <- file.path(OUTDIR, "Tables")

invisible(lapply(
  c(OUTDIR, DIR_PANEL, DIR_GSEA, DIR_TABLE),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

CELLTYPE_COL <- "celltype"
CONDITION_CANDIDATES <- c("tissue_type", "condition")
SAMPLE_CANDIDATES <- c("sample_id", "orig.ident")
ENDO_LABEL <- "Endothelial cell"

MIN_ENDO_CELLS_PER_SAMPLE <- 20
DE_FDR_CUTOFF <- 0.05
DE_LOGFC_CUTOFF <- 0.50
ORA_MIN_GENE_NUMBER <- 15
ORA_FALLBACK_TOP_N <- 250
GSEA_MIN_SIZE <- 5
GSEA_MAX_SIZE <- 500
RANDOM_SEED <- 1234

TARGET_KEGG <- data.frame(
  KEGG_ID = c(
    "hsa04218",
    "hsa04066",
    "hsa04512",
    "hsa04520"
  ),
  KEGG_Name = c(
    "Cellular senescence",
    "HIF-1 signaling pathway",
    "ECM-receptor interaction",
    "Adherens junction"
  ),
  Short_Name = c(
    "Cellular\nsenescence",
    "HIF-1 signaling\npathway",
    "ECM-receptor\ninteraction",
    "Adherens\njunction"
  ),
  stringsAsFactors = FALSE
)

PATHWAY_COLORS <- c(
  "hsa04218" = "#E67E22",
  "hsa04066" = "#2F6FB0",
  "hsa04512" = "#D63384",
  "hsa04520" = "#169C9C"
)

CONDITION_COLORS <- c(
  "CIT" = "#2166AC",
  "HR" = "#D73027"
)

CELLTYPE_COLORS <- c(
  "AT1" = "#168AAD",
  "AT2" = "#F57C20",
  "Other epithelial" = "#F46D43",
  "Ciliated cell" = "#313695",
  "Endothelial cell" = "#D73027",
  "Fibroblast" = "#4DAF4A",
  "Macrophage" = "#8064A2",
  "Monocyte" = "#9C755F",
  "Neutrophil" = "#F781BF",
  "T cell" = "#2CA6A4",
  "NK cell" = "#51438F",
  "B cell" = "#2468B4",
  "Mast cell" = "#C9A227",
  "Cycling cell" = "#E69F00",
  "Unknown" = "#8C8C8C",
  "NA" = "#8C8C8C"
)

SERPINE_PALETTE <- c(
  "#E5E9ED",
  "#F7C7B8",
  "#F47C68",
  "#D73027",
  "#8E1021"
)

# =============================================================================
# =============================================================================

cran_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "patchwork",
  "scales",
  "ragg"
)

bioc_packages <- c(
  "edgeR",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "KEGGREST",
  "fgsea"
)

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(ragg)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(KEGGREST)
  library(edgeR)
  library(fgsea)
})

# =============================================================================
# =============================================================================

if (!exists("sce", envir = .GlobalEnv, inherits = FALSE)) {
  stop(" R sce ?")
}

sce_fig4 <- get("sce", envir = .GlobalEnv)

if (!inherits(sce_fig4, "Seurat")) {
  stop("sce Seurat ?")
}

if (!CELLTYPE_COL %in% colnames(sce_fig4@meta.data)) {
  stop("sce@meta.data celltype ?")
}

condition_col <- CONDITION_CANDIDATES[CONDITION_CANDIDATES %in% colnames(sce_fig4@meta.data)][1]

sample_col <- SAMPLE_CANDIDATES[SAMPLE_CANDIDATES %in% colnames(sce_fig4@meta.data)][1]

if (is.na(condition_col)) {
  stop("sce@meta.data tissue_type condition ?")
}

if (is.na(sample_col)) {
  stop("sce@meta.data sample_id orig.ident ?")
}

if (!"umap" %in% Seurat::Reductions(sce_fig4)) {
  stop("sce umap ?")
}

Seurat::DefaultAssay(sce_fig4) <- "RNA"

if (!"SERPINE1" %in% rownames(sce_fig4)) {
  stop("RNA assay SERPINE1 ?")
}

condition_values <- unique(as.character(
  sce_fig4@meta.data[[condition_col]]
))

if (!all(c("CIT", "HR") %in% condition_values)) {
  stop(
    " CIT HR. ",
    paste(condition_values, collapse = ", ")
  )
}

# =============================================================================
# =============================================================================

fetch_data_compat <- function(object, vars, assay = "RNA") {
  Seurat::DefaultAssay(object) <- assay

  result <- tryCatch(
    Seurat::FetchData(
      object,
      vars = vars,
      layer = "data"
    ),
    error = function(e) NULL
  )

  if (is.null(result)) {
    result <- Seurat::FetchData(
      object,
      vars = vars,
      slot = "data"
    )
  }

  result
}

get_assay_matrix_compat <- function(
    object,
    assay = "RNA",
    layer_name = "counts"
) {
  result <- tryCatch(
    SeuratObject::GetAssayData(
      object,
      assay = assay,
      layer = layer_name
    ),
    error = function(e) NULL
  )

  if (is.null(result)) {
    result <- tryCatch(
      SeuratObject::GetAssayData(
        object,
        assay = assay,
        slot = layer_name
      ),
      error = function(e) NULL
    )
  }

  if (is.null(result) &&
      exists("JoinLayers", where = asNamespace("SeuratObject"),
             mode = "function")) {
    object_copy <- object
    object_copy[[assay]] <- SeuratObject::JoinLayers(
      object_copy[[assay]]
    )

    result <- tryCatch(
      SeuratObject::GetAssayData(
        object_copy,
        assay = assay,
        layer = layer_name
      ),
      error = function(e) NULL
    )
  }

  if (is.null(result)) {
    stop(
      "Check source data or configuration. ",
      assay,
      " assay ",
      layer_name,
      "Check source data or configuration. "
    )
  }

  result
}

format_probability <- function(x) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) {
    return("NA")
  }

  if (x < 0.001) {
    return("<0.001")
  }

  sprintf("%.3f", x)
}

theme_publication <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 2,
        hjust = 0.5,
        color = "black",
        margin = ggplot2::margin(b = 5)
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 0.5,
        hjust = 0.5,
        color = "black",
        margin = ggplot2::margin(b = 6)
      ),
      axis.title = ggplot2::element_text(
        face = "bold",
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        color = "black"
      ),
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = 0.8
      ),
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = 0.7
      ),
      panel.grid.major = ggplot2::element_line(
        color = "#E8EDF0",
        linewidth = 0.32
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      panel.border = ggplot2::element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.7
      ),
      strip.background = ggplot2::element_rect(
        fill = "#F1F5F7",
        color = "black",
        linewidth = 0.7
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        color = "black"
      ),
      legend.title = ggplot2::element_text(
        face = "bold",
        color = "black"
      ),
      legend.text = ggplot2::element_text(
        color = "black"
      ),
      plot.margin = ggplot2::margin(7, 8, 7, 7)
    )
}

theme_umap <- function(base_size = 10) {
  theme_publication(base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank()
    )
}

save_plot3 <- function(
    plot_object,
    filename,
    folder,
    width,
    height,
    dpi = 600
) {
  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggplot2::ggsave(
    filename = file.path(folder, paste0(filename, ".pdf")),
    plot = plot_object,
    width = width,
    height = height,
    device = pdf_device,
    bg = "white",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    filename = file.path(folder, paste0(filename, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = ragg::agg_png,
    bg = "white",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    filename = file.path(folder, paste0(filename, ".tiff")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = "tiff",
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
}

safe_wilcox <- function(data, value_col, group_col) {
  data <- data[is.finite(data[[value_col]]) & !is.na(data[[group_col]]), , drop = FALSE]

  if (length(unique(data[[group_col]])) != 2) {
    return(NA_real_)
  }

  if (any(table(data[[group_col]]) < 2)) {
    return(NA_real_)
  }

  tryCatch(
    stats::wilcox.test(
      data[[value_col]] ~ data[[group_col]],
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 5 ||
      length(unique(x)) < 3 ||
      length(unique(y)) < 3) {
    return(data.frame(
      Rho = NA_real_,
      P = NA_real_,
      N = length(x),
      stringsAsFactors = FALSE
    ))
  }

  result <- suppressWarnings(
    stats::cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )

  data.frame(
    Rho = unname(result$estimate),
    P = result$p.value,
    N = length(x),
    stringsAsFactors = FALSE
  )
}

safe_adjusted_lm <- function(data) {
  data <- data[is.finite(data$SERPINE1_Mean) & is.finite(data$PathwayScore_Mean) & !is.na(data$Condition), , drop = FALSE]

  if (nrow(data) < 8 ||
      stats::sd(data$SERPINE1_Mean) == 0 ||
      stats::sd(data$PathwayScore_Mean) == 0) {
    return(data.frame(
      Adjusted_Beta = NA_real_,
      Adjusted_P = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  data$SERPINE1_z <- as.numeric(scale(data$SERPINE1_Mean))
  data$Pathway_z <- as.numeric(scale(data$PathwayScore_Mean))

  fit <- tryCatch(
    stats::lm(
      Pathway_z ~ SERPINE1_z + Condition,
      data = data
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(data.frame(
      Adjusted_Beta = NA_real_,
      Adjusted_P = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  coef_table <- summary(fit)$coefficients

  if (!"SERPINE1_z" %in% rownames(coef_table)) {
    return(data.frame(
      Adjusted_Beta = NA_real_,
      Adjusted_P = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    Adjusted_Beta = coef_table["SERPINE1_z", "Estimate"],
    Adjusted_P = coef_table["SERPINE1_z", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}

ratio_to_numeric <- function(x) {
  vapply(
    strsplit(as.character(x), "/", fixed = TRUE),
    function(z) {
      if (length(z) != 2) {
        return(NA_real_)
      }

      as.numeric(z[1]) / as.numeric(z[2])
    },
    numeric(1)
  )
}

symbol_to_entrez <- function(symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]

  if (length(symbols) == 0) {
    return(character())
  }

  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )

  unique(unname(
    mapped[!is.na(mapped) & nzchar(mapped)]
  ))
}

get_official_kegg_symbols <- function(pathway_id, pathway_name) {
  entry <- tryCatch(
    KEGGREST::keggGet(pathway_id)[[1]],
    error = function(e) NULL
  )

  if (!is.null(entry)) {
    gene_field <- entry$GENE

    if (!is.null(gene_field) && length(gene_field) > 0) {
      if (!is.null(names(gene_field)) &&
          any(nzchar(names(gene_field)))) {
        entrez_ids <- names(gene_field)
      } else {
        entrez_ids <- gene_field[seq(1, length(gene_field), by = 2)]
      }

      entrez_ids <- sub("^hsa:", "", as.character(entrez_ids))
      entrez_ids <- sub("[^0-9].*$", "", entrez_ids)
      entrez_ids <- unique(
        entrez_ids[grepl("^[0-9]+$", entrez_ids)]
      )

      symbol_map <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = entrez_ids,
        keytype = "ENTREZID",
        column = "SYMBOL",
        multiVals = "first"
      )

      symbols <- unique(unname(
        symbol_map[!is.na(symbol_map) & nzchar(symbol_map)]
      ))

      if (length(symbols) >= GSEA_MIN_SIZE) {
        return(symbols)
      }
    }
  }

  # Source comment omitted during English-language consolidation.
  fallback_objects <- c(
    "official_kegg_genes",
    "Official_KEGG_pathway_genes"
  )

  for (object_name in fallback_objects) {
    if (exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
      object_i <- get(object_name, envir = .GlobalEnv)

      if (is.data.frame(object_i) &&
          all(c("KEGG_ID", "Gene") %in% colnames(object_i))) {
        genes_i <- unique(as.character(
          object_i$Gene[object_i$KEGG_ID == pathway_id]
        ))

        genes_i <- genes_i[!is.na(genes_i) & nzchar(genes_i)]

        if (length(genes_i) >= GSEA_MIN_SIZE) {
          return(genes_i)
        }
      }
    }
  }

  stop(
    " KEGG  ",
    pathway_id,
    "Check source data or configuration. ",
    pathway_name,
    "） KEGGREST ?"
  )
}

add_module_score_compat <- function(
    object,
    genes,
    score_name,
    seed = 1234
) {
  genes <- unique(intersect(genes, rownames(object)))

  if (length(genes) < GSEA_MIN_SIZE) {
    stop(
      score_name,
      "Check source data or configuration. ",
      GSEA_MIN_SIZE,
      "Check source data or configuration. "
    )
  }

  if (score_name %in% colnames(object@meta.data)) {
    object@meta.data[[score_name]] <- NULL
  }

  old_columns <- colnames(object@meta.data)

  result <- tryCatch(
    Seurat::AddModuleScore(
      object = object,
      features = list(genes),
      assay = "RNA",
      ctrl = 50,
      nbin = 24,
      seed = seed,
      name = paste0(score_name, "_tmp_")
    ),
    error = function(e) NULL
  )

  if (!is.null(result)) {
    new_columns <- setdiff(
      colnames(result@meta.data),
      old_columns
    )

    if (length(new_columns) == 1) {
      names(result@meta.data)[names(result@meta.data) == new_columns] <- score_name

      return(result)
    }
  }

  # Source comment omitted during English-language consolidation.
  normalized_matrix <- get_assay_matrix_compat(
    object,
    assay = "RNA",
    layer_name = "data"
  )

  object[[score_name]] <- Matrix::colMeans(
    normalized_matrix[genes, , drop = FALSE]
  )

  object
}

make_placeholder_plot <- function(title, message) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = message,
      size = 4,
      color = "black"
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      panel.border = ggplot2::element_rect(
        fill = NA,
        color = "black"
      )
    )
}

# =============================================================================
# =============================================================================

umap_matrix <- Seurat::Embeddings(
  sce_fig4,
  reduction = "umap"
)

umap_df <- data.frame(
  Cell = rownames(umap_matrix),
  UMAP_1 = as.numeric(umap_matrix[, 1]),
  UMAP_2 = as.numeric(umap_matrix[, 2]),
  stringsAsFactors = FALSE
)

serpine_df <- fetch_data_compat(
  sce_fig4,
  vars = "SERPINE1"
)

cell_df <- sce_fig4@meta.data %>%
  dplyr::mutate(
    Cell = rownames(sce_fig4@meta.data),
    Celltype = as.character(.data[[CELLTYPE_COL]]),
    Condition = as.character(.data[[condition_col]]),
    Sample = as.character(.data[[sample_col]])
  ) %>%
  dplyr::select(
    Cell,
    Celltype,
    Condition,
    Sample
  ) %>%
  dplyr::left_join(
    umap_df,
    by = "Cell"
  )

cell_df$SERPINE1 <- as.numeric(
  serpine_df[cell_df$Cell, "SERPINE1"]
)

cell_df$SERPINE1_Positive <- cell_df$SERPINE1 > 0

celltype_order_template <- c(
  "AT1",
  "AT2",
  "Other epithelial",
  "Ciliated cell",
  "Endothelial cell",
  "Fibroblast",
  "Macrophage",
  "Monocyte",
  "Neutrophil",
  "T cell",
  "NK cell",
  "B cell",
  "Mast cell",
  "Cycling cell",
  "Unknown",
  "NA"
)

celltype_levels <- c(
  intersect(
    celltype_order_template,
    unique(cell_df$Celltype)
  ),
  sort(setdiff(
    unique(cell_df$Celltype),
    celltype_order_template
  ))
)

cell_df$Celltype <- factor(
  cell_df$Celltype,
  levels = celltype_levels
)

cell_df$Condition <- factor(
  cell_df$Condition,
  levels = c("CIT", "HR")
)

celltype_colors_use <- CELLTYPE_COLORS[celltype_levels]

missing_color_types <- names(celltype_colors_use)[is.na(celltype_colors_use)]

if (length(missing_color_types) > 0) {
  extra_colors <- grDevices::hcl.colors(
    length(missing_color_types),
    palette = "Dark 3"
  )

  names(extra_colors) <- missing_color_types

  celltype_colors_use[missing_color_types] <- extra_colors
}

# =============================================================================
# =============================================================================

cell_centers <- cell_df %>%
  dplyr::group_by(Celltype) %>%
  dplyr::summarise(
    UMAP_1 = stats::median(UMAP_1, na.rm = TRUE),
    UMAP_2 = stats::median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

p4A <- ggplot2::ggplot(
  cell_df,
  ggplot2::aes(
    x = UMAP_1,
    y = UMAP_2,
    color = Celltype
  )
) +
  ggplot2::geom_point(
    size = 0.18,
    alpha = 0.78
  ) +
  ggrepel::geom_text_repel(
    data = cell_centers,
    ggplot2::aes(
      x = UMAP_1,
      y = UMAP_2,
      label = Celltype
    ),
    color = "black",
    fontface = "bold",
    size = 3,
    box.padding = 0.35,
    point.padding = 0.15,
    segment.color = "#7D8A93",
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = celltype_colors_use,
    drop = FALSE
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "All-cell UMAP",
    x = "UMAP_1",
    y = "UMAP_2",
    color = NULL
  ) +
  theme_umap(9.5) +
  ggplot2::theme(
    legend.position = "none"
  )

# =============================================================================
# =============================================================================

p4B <- ggplot2::ggplot(
  cell_df,
  ggplot2::aes(
    x = UMAP_1,
    y = UMAP_2,
    color = Celltype
  )
) +
  ggplot2::geom_point(
    size = 0.15,
    alpha = 0.78
  ) +
  ggplot2::facet_wrap(
    ~ Condition,
    nrow = 1
  ) +
  ggplot2::scale_color_manual(
    values = celltype_colors_use,
    drop = FALSE
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "CIT and HR UMAP",
    x = "UMAP_1",
    y = "UMAP_2",
    color = NULL
  ) +
  theme_umap(9.2) +
  ggplot2::theme(
    legend.position = "right",
    legend.text = ggplot2::element_text(size = 7.5)
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 2.6, alpha = 1),
      ncol = 1
    )
  )

# =============================================================================
# =============================================================================

marker_map <- c(
  "AT1" = "AGER",
  "AT2" = "SFTPC",
  "Other epithelial" = "EPCAM",
  "Ciliated cell" = "FOXJ1",
  "Endothelial cell" = "PECAM1",
  "Fibroblast" = "COL1A1",
  "Macrophage" = "C1QA",
  "Monocyte" = "FCN1",
  "Neutrophil" = "CSF3R",
  "T cell" = "CD3D",
  "NK cell" = "NKG7",
  "B cell" = "CD79A",
  "Mast cell" = "TPSAB1",
  "Cycling cell" = "MKI67"
)

marker_map <- marker_map[names(marker_map) %in% celltype_levels]

marker_map <- marker_map[marker_map %in% rownames(sce_fig4)]

marker_features <- unname(marker_map)

dot_object <- Seurat::DotPlot(
  object = sce_fig4,
  features = marker_features,
  group.by = CELLTYPE_COL,
  assay = "RNA",
  scale = TRUE,
  dot.scale = 7
)

marker_plot_df <- dot_object$data %>%
  dplyr::mutate(
    Gene = as.character(features.plot),
    Celltype = as.character(id)
  ) %>%
  dplyr::filter(
    Gene %in% marker_features,
    Celltype %in% celltype_levels
  )

marker_plot_df$Gene <- factor(
  marker_plot_df$Gene,
  levels = rev(marker_features)
)

marker_plot_df$Celltype <- factor(
  marker_plot_df$Celltype,
  levels = celltype_levels
)

p4C <- ggplot2::ggplot(
  marker_plot_df,
  ggplot2::aes(
    x = Celltype,
    y = Gene
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = pct.exp,
      color = avg.exp.scaled
    ),
    alpha = 0.96
  ) +
  ggplot2::scale_size_continuous(
    range = c(0.6, 6.5),
    limits = c(0, 100),
    breaks = c(25, 50, 75, 100),
    name = "Percent\nexpressing"
  ) +
  ggplot2::scale_color_gradient2(
    low = "#3B5BA5",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    name = "Average\nexpression"
  ) +
  ggplot2::labs(
    title = "Cell-type marker DotPlot",
    x = NULL,
    y = NULL
  ) +
  theme_publication(8.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 52,
      hjust = 1,
      vjust = 1,
      size = 7.5
    ),
    axis.text.y = ggplot2::element_text(
      face = "italic",
      size = 8.2
    ),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "right"
  )

# =============================================================================
# 9. Panel D:SERPINE1 FeaturePlot
# =============================================================================

serpine_upper <- as.numeric(stats::quantile(
  cell_df$SERPINE1,
  probs = 0.98,
  na.rm = TRUE
))

if (!is.finite(serpine_upper) || serpine_upper <= 0) {
  serpine_upper <- max(cell_df$SERPINE1, na.rm = TRUE)
}

if (!is.finite(serpine_upper) || serpine_upper <= 0) {
  serpine_upper <- 1
}

serpine_umap_df <- cell_df %>%
  dplyr::arrange(SERPINE1)

p4D <- ggplot2::ggplot(
  serpine_umap_df,
  ggplot2::aes(
    x = UMAP_1,
    y = UMAP_2
  )
) +
  ggplot2::geom_point(
    color = "#E2E6E9",
    size = 0.16,
    alpha = 0.52
  ) +
  ggplot2::geom_point(
    data = serpine_umap_df %>%
      dplyr::filter(SERPINE1 > 0),
    ggplot2::aes(color = SERPINE1),
    size = 0.22,
    alpha = 0.92
  ) +
  ggplot2::scale_color_gradientn(
    colors = SERPINE_PALETTE,
    trans = "sqrt",
    limits = c(0, serpine_upper),
    oob = scales::squish,
    name = "SERPINE1\nexpression"
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = expression(italic(SERPINE1) ~ "FeaturePlot"),
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_umap(9.5)

# =============================================================================
# =============================================================================

serpine_summary <- cell_df %>%
  dplyr::group_by(
    Celltype,
    Condition
  ) %>%
  dplyr::summarise(
    MeanExpression = mean(SERPINE1, na.rm = TRUE),
    PercentPositive = mean(
      SERPINE1_Positive,
      na.rm = TRUE
    ) * 100,
    CellNumber = dplyr::n(),
    .groups = "drop"
  )

serpine_summary$Celltype <- factor(
  serpine_summary$Celltype,
  levels = rev(celltype_levels)
)

serpine_summary$Condition <- factor(
  serpine_summary$Condition,
  levels = c("CIT", "HR")
)

serpine_mean_upper <- as.numeric(stats::quantile(
  serpine_summary$MeanExpression,
  0.98,
  na.rm = TRUE
))

if (!is.finite(serpine_mean_upper) ||
    serpine_mean_upper <= 0) {
  serpine_mean_upper <- max(
    serpine_summary$MeanExpression,
    na.rm = TRUE
  )
}

p4E <- ggplot2::ggplot(
  serpine_summary,
  ggplot2::aes(
    x = Condition,
    y = Celltype
  )
) +
  ggplot2::geom_tile(
    data = serpine_summary %>%
      dplyr::filter(
        as.character(Celltype) == ENDO_LABEL
      ),
    fill = "#FCE5E3",
    color = NA,
    width = 0.88,
    height = 0.78
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = PercentPositive,
      fill = MeanExpression
    ),
    shape = 21,
    color = "#263640",
    stroke = 0.45,
    alpha = 0.97
  ) +
  ggplot2::scale_size_continuous(
    range = c(1.5, 9.2),
    limits = c(0, 100),
    breaks = c(25, 50, 75, 100),
    name = "Percent\npositive"
  ) +
  ggplot2::scale_fill_gradientn(
    colors = SERPINE_PALETTE,
    trans = "sqrt",
    limits = c(0, serpine_mean_upper),
    oob = scales::squish,
    name = "Mean\nSERPINE1"
  ) +
  ggplot2::labs(
    title = expression(
      italic(SERPINE1) ~ "by cell type and condition"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_publication(9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      face = "bold",
      color = "black"
    ),
    axis.text.y = ggplot2::element_text(size = 8.3),
    axis.ticks = ggplot2::element_blank()
  )

# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================



md <- sce@meta.data
counts <- as.data.frame(table(Sample=md$orig.ident,Celltype=md$celltype))
names(counts)[3] <- "Cells"
counts$Condition <- sub("[1-6]$","",counts$Sample)
counts$Fraction <- counts$Cells/ave(counts$Cells,counts$Sample,FUN=sum)
p_composition <- ggplot2::ggplot(counts,ggplot2::aes(Sample,Fraction,fill=Celltype))+ggplot2::geom_col()+ggplot2::theme_classic()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1))+ggplot2::labs(x=NULL,y="Fraction within sample")
panels <- list(A=p4A,B=p4B,C=p4C,D=p_composition,E=p4D,F=p4E)
for (nm in names(panels)) save_publication_plot(panels[[nm]],file.path(repo_output("Figure3/Atlas"),paste0("Figure3",nm)),8,6)
saveRDS(panels,file.path(repo_output("Figure3/Atlas"),"Figure3_A_F_panels.rds"))


}

# =============================================================================
# STAGE figure3_paired_compute
# =============================================================================
PIPELINE$figure3_paired_compute <- function() {
# Figure 3 G-K: biological-sample summaries and donor-paired edgeR pseudobulk.
require_packages(c("SeuratObject","Matrix","edgeR","fgsea","msigdbr"))
sce<-load_sce();md<-sce@meta.data
if(!all(c("celltype","orig.ident")%in%names(md)))stop("Atlas metadata lack celltype/orig.ident.")
keep<-rownames(md)[as.character(md$celltype)==CONFIG$endothelial_label]
if(!length(keep))stop("No endothelial cells match CONFIG$endothelial_label.")
map<-read.csv(file.path(REPO_ROOT,"metadata/single_cell_sample_map.csv"))
get_rna<-function(layer)tryCatch(SeuratObject::GetAssayData(sce,assay="RNA",layer=layer),error=function(e)SeuratObject::GetAssayData(sce,assay="RNA",slot=layer))
counts<-get_rna("counts")[,keep,drop=FALSE];expr<-get_rna("data")[,keep,drop=FALSE]
sample<-as.character(md[keep,"orig.ident"])
if(!setequal(unique(sample),map$Sample))stop("Expected all 12 biological samples.")
G<-do.call(rbind,lapply(seq_len(nrow(map)),function(i){
 cells<-keep[sample==map$Sample[i]];v<-expr["SERPINE1",cells]
 data.frame(map[i,],MeanSERPINE1=mean(v),PositivePercent=100*mean(v>0))
}))
metrics<-c("MeanSERPINE1","PositivePercent")
GS<-do.call(rbind,lapply(metrics,function(m){
 cit<-G[G$Condition=="CIT",c("Donor",m)];hr<-G[G$Condition=="HR",c("Donor",m)]
 z<-merge(cit,hr,by="Donor",suffixes=c(".CIT",".HR"))
 data.frame(Metric=m,PairedWilcoxon_P=wilcox.test(z[[paste0(m,".HR")]],z[[paste0(m,".CIT")]],paired=TRUE,alternative="two.sided",exact=FALSE)$p.value)
}))
GS$BH_FDR_two_metrics<-p.adjust(GS$PairedWilcoxon_P,"BH")
indicator<-Matrix::sparseMatrix(i=seq_along(sample),j=match(sample,map$Sample),x=1,dims=c(length(sample),nrow(map)))
pb<-counts%*%indicator;colnames(pb)<-map$Sample
map$Condition<-factor(map$Condition,levels=c("CIT","HR"));map$Donor<-factor(map$Donor)
design<-model.matrix(~Donor+Condition,map)
y<-edgeR::DGEList(counts=pb);y<-y[edgeR::filterByExpr(y,design=design),,keep.lib.sizes=FALSE]
y<-edgeR::calcNormFactors(y);y<-edgeR::estimateDisp(y,design,robust=TRUE)
fit<-edgeR::glmQLFit(y,design,robust=TRUE);q<-edgeR::glmQLFTest(fit,coef="ConditionHR")
tab<-edgeR::topTags(q,n=Inf,sort.by="none")$table;tab$Gene<-rownames(tab)
tab$Signed_sqrt_QLF<-sign(tab$logFC)*sqrt(tab$F)
out<-repo_output("Figure3/Paired/Tables")
write.csv(tab,file.path(out,"05_edgeR_all_tested_genes.csv"),row.names=FALSE)
write.csv(tab[c("Gene","Signed_sqrt_QLF")],file.path(out,"06_paired_signed_rank_for_GSEA.csv"),row.names=FALSE)
write.csv(G,file.path(out,"07_Fig3G_EC_SERPINE1_6pairs.csv"),row.names=FALSE)
write.csv(GS,file.path(out,"08_Fig3G_SERPINE1_paired_P_FDR.csv"),row.names=FALSE)
groups<-list("Injury/senescence"=c("SERPINE1","CDKN1A","GDF15"),"Inflammatory adhesion"=c("ICAM1","VCAM1","SELE","SELP"),"Permeability/remodeling"=c("PLVAP","ANGPT2","FN1","COL4A1"),"Barrier homeostasis"=c("CLDN5","CDH5","TJP1","OCLN","TEK","KDR","ROBO4"))
annot<-data.frame(Gene=unlist(groups,use.names=FALSE),GeneGroup=rep(names(groups),lengths(groups)))
annot<-annot[annot$Gene%in%rownames(y),]
lc<-edgeR::cpm(y,log=TRUE,prior.count=2)
H<-do.call(rbind,lapply(seq_len(nrow(annot)),function(i){v<-lc[annot$Gene[i],];sdv<-sd(v);data.frame(Gene=annot$Gene[i],Sample=names(v),Expression=as.numeric(v),Condition=as.character(map$Condition),GeneGroup=annot$GeneGroup[i],Z=if(sdv>0)as.numeric(scale(v)) else 0)}))
write.csv(H,file.path(out,"09_Fig3H_paired_pseudobulk_logCPM_matrix.csv"),row.names=FALSE)
write.csv(merge(annot,tab[c("Gene","logFC","PValue","FDR")],by="Gene"),file.path(out,"10_Fig3H_gene_logFC_BH_FDR.csv"),row.names=FALSE)
# The exact gene-set release affects enrichment estimates. Save the complete
# annotation and result tables; do not overwrite calculated NES/FDR with text.
sets<-msigdbr::msigdbr(species="Homo sapiens")
write.csv(sets,file.path(out,"gene_set_annotation_snapshot.csv"),row.names=FALSE)
ranked<-setNames(tab$Signed_sqrt_QLF,tab$Gene);ranked<-sort(ranked,decreasing=TRUE)
run_gsea<-function(z){
 if(!nrow(z))stop("Required gene-set collection is absent from this msigdbr release.")
 p<-split(z$gene_symbol,z$gs_name);set.seed(1234)
 r<-as.data.frame(fgsea::fgseaMultilevel(pathways=p,stats=ranked,minSize=5,maxSize=500))
 r$leadingEdge<-vapply(r$leadingEdge,paste,collapse=";",FUN.VALUE=character(1));r
}
collection<-if("gs_cat"%in%names(sets))sets$gs_cat else sets$gs_collection
subcategory<-if("gs_subcat"%in%names(sets))sets$gs_subcat else sets$gs_subcollection
I<-run_gsea(sets[collection=="C2"&grepl("KEGG",subcategory),]);J<-run_gsea(sets[collection=="C5"&subcategory=="GO:BP",])
write.csv(I,file.path(out,"11_Fig3I_KEGG_paired_GSEA_all_tested.csv"),row.names=FALSE)
write.csv(J,file.path(out,"12_Fig3J_GO_BP_paired_GSEA_all_tested.csv"),row.names=FALSE)
targets<-c("CELLULAR_SENESCENCE","HIF_1_SIGNALING_PATHWAY","ECM_RECEPTOR_INTERACTION","ADHERENS_JUNCTION")
K<-do.call(rbind,lapply(targets,function(t){
 h<-grep(paste0(t,"$"),I$pathway);if(length(h)>1)stop("Ambiguous target pathway: ",t)
 data.frame(Pathway=gsub("_"," ",t),Found=length(h)==1,NES=if(length(h))I$NES[h] else NA_real_,P=if(length(h))I$pval[h] else NA_real_,FDR=if(length(h))I$padj[h] else NA_real_)
}))
write.csv(K,file.path(out,"13_Fig3K_four_Kegg_pathways_and_missing_status.csv"),row.names=FALSE)
saveRDS(list(y=y,design=design,fit=fit,test=q,results=tab),file.path(out,"paired_edgeR_checkpoint.rds"))

}

# =============================================================================
# STAGE figure3_paired_plots
# =============================================================================
PIPELINE$figure3_paired_plots <- function() {
# scripts/05_paired_endothelial_plots.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

# =============================================================================
# FIGURE 3 G-K | PREMIUM RE-PLOT FROM EXISTING PAIRED ANALYSIS TABLES
# R 4.3.2 | ggplot2 3.5.1 | patchwork 1.2.0 | ragg 1.3.2
# Preservation: same CIT blue, HR red, B/W/red expression gradient,
# bold BLACK text, donor paired design, HR-CIT direction.
# NO edgeR/limma/fgsea reruns, NO re-normalization, NO Seurat object needed.
# Input: Tables/ directory from FIG3_REVIEW_PAIRED_ONLY_G_to_K.zip
# Output: PDF (vector) + PNG (600 dpi) for G/H/I/J/K and a G-K composite.
# IMPORTANT: old full Figure3 A-F panels are NOT recalculated or overwritten.
# =============================================================================
options(stringsAsFactors = FALSE)
set.seed(20261022L)

# ---- 0. DATA LOCATION -------------------------------------------------------
# Use this if your results are already extracted somewhere:
TABLE_DIR <- repo_input("figure3_paired_tables")

# Find script directory when invoked via source(), otherwise working directory.
get_source_dir <- function() {
  fs <- sys.frames()
  for (j in rev(seq_along(fs))) {
    o <- fs[[j]]$ofile
    if (!is.null(o) && file.exists(o)) return(dirname(normalizePath(o)))
  }
  getwd()
}
SCRIPT_DIR <- get_source_dir()
if (!nzchar(TABLE_DIR)) {
  roots <- unique(c(SCRIPT_DIR, getwd()))
  cand <- unique(unlist(lapply(roots, function(r) c(file.path(r, "Tables"),
      file.path(r, "FIG3_REVIEW_PAIRED_ONLY_G_to_K", "Tables")))))
  have <- cand[file.exists(file.path(cand, "07_Fig3G_EC_SERPINE1_6pairs.csv"))]
  if (length(have)) TABLE_DIR <- have[1]
}
if (!nzchar(TABLE_DIR) || !dir.exists(TABLE_DIR)) {
  message("Please select 07_Fig3G_EC_SERPINE1_6pairs.csv in the original Tables folder.")
  selected_file <- file.choose()
  if (basename(selected_file) != "07_Fig3G_EC_SERPINE1_6pairs.csv")
    stop("Select the existing 07_Fig3G_EC_SERPINE1_6pairs.csv, not another file")
  TABLE_DIR <- dirname(selected_file)
}
TABLE_DIR <- normalizePath(TABLE_DIR, winslash = "/", mustWork = TRUE)
OUT <- file.path(dirname(TABLE_DIR), "FIG3_PREMIUM_G_K_REPLOT")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Please install patchwork")
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

# Original visual identity, stronger strokes and clearer black labels.
CIT <- "#2878D4"
HR  <- "#E3322B"
BLUE_LIGHT <- "#BDD9EE"
RED_LIGHT <- "#F3B7B4"
INK <- "#111111"
GRID <- "#E9EDF0"
WHITE <- "#FFFFFF"

read_table <- function(filename, required) {
  p <- file.path(TABLE_DIR, filename)
  if (!file.exists(p)) stop("Missing paired analysis table:\n", p)
  x <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE,
                       fileEncoding = "UTF-8-BOM")
  miss <- setdiff(required, colnames(x))
  if (length(miss)) stop(filename, " missing columns: ",paste(miss,collapse=", "))
  x
}

# Existing statistical results only (no test or FDR is recomputed).
Gdata <- read_table("07_Fig3G_EC_SERPINE1_6pairs.csv",
                    c("Sample", "Donor", "Condition", "MeanSERPINE1", "PositivePercent"))
Gstat <- read_table("08_Fig3G_SERPINE1_paired_P_FDR.csv",
                    c("Metric", "PairedWilcoxon_P", "BH_FDR_two_metrics"))
Hdata <- read_table("09_Fig3H_paired_pseudobulk_logCPM_matrix.csv",
                    c("Gene", "Sample", "Expression", "Condition", "GeneGroup", "Z"))
Hstat <- read_table("10_Fig3H_gene_logFC_BH_FDR.csv",
                    c("Gene", "logFC", "PValue", "FDR", "GeneGroup"))
Idata <- read_table("11_Fig3I_KEGG_paired_GSEA_all_tested.csv",
                    c("pathway", "NES", "padj", "size"))
Jdata <- read_table("12_Fig3J_GO_BP_paired_GSEA_all_tested.csv",
                    c("pathway", "NES", "padj", "size"))
Kdata <- read_table("13_Fig3K_four_Kegg_pathways_and_missing_status.csv",
                    c("Pathway", "Found", "NES", "P", "FDR"))

# Mathematical safety checks and locked samples.
expected <- c(paste0("CIT", 1:6), paste0("HR", 1:6))
stopifnot(nrow(Gdata) == 12L, nrow(Gstat) == 2L,
          setequal(as.character(Gdata$Sample),expected))
Gdata$Donor <- as.character(Gdata$Donor)
Gdata$Condition <- factor(as.character(Gdata$Condition), levels=c("CIT","HR"))
if (anyNA(Gdata$Condition) || length(unique(Gdata$Donor)) != 6L ||
    any(table(Gdata$Donor, Gdata$Condition) != 1L))
  stop("Figure G is not six complete donor pairs")
if (anyNA(Gdata[,c("MeanSERPINE1","PositivePercent")]) ||
    any(!is.finite(Gdata$MeanSERPINE1)) ||
    any(!is.finite(Gdata$PositivePercent))) stop("Nonfinite Figure G values")
if (!setequal(unique(Hdata$Sample), expected) ||
    anyDuplicated(paste(Hdata$Gene,Hdata$Sample)))
  stop("Figure H pseudobulk sample/gene IDs differ from expected 12 samples")
if (!setequal(Hdata$Gene,Hstat$Gene)) stop("Figure H gene statistics mismatch")
if (anyNA(Hstat$FDR) || any(!is.finite(Hstat$logFC))) stop("Invalid Figure H stats")
if (any(!is.finite(Gstat$PairedWilcoxon_P)) ||
    any(!is.finite(Gstat$BH_FDR_two_metrics))) stop("Invalid Figure G stats")

# ---- Publication typography -------------------------------------------------
font_theme <- function(base=11) {
  theme_classic(base_size=base, base_family="sans") +
    theme(
      text=element_text(family="sans",face="bold",color=INK),
      plot.title=element_text(size=base+1.5,face="bold",color=INK,
                              hjust=.5,margin=margin(b=6)),
      plot.subtitle=element_text(size=base-1,face="bold",color=INK,
                                 hjust=.5,margin=margin(b=9)),
      plot.caption=element_text(size=base-2.3,color=INK,face="bold",hjust=0,
                                margin=margin(t=7)),
      axis.title=element_text(size=base,face="bold",color=INK),
      axis.text=element_text(size=base-.9,face="bold",color=INK),
      axis.line=element_line(color=INK,linewidth=.62),
      axis.ticks=element_line(color=INK,linewidth=.53),
      panel.grid.major.y=element_line(color=GRID,linewidth=.32),
      panel.grid.major.x=element_blank(),
      panel.grid.minor=element_blank(),
      strip.text=element_text(color=INK,face="bold",size=base-.7),
      strip.background=element_rect(fill="#F4F6F8",color=NA),
      legend.text=element_text(face="bold",color=INK,size=base-1.3),
      legend.title=element_text(face="bold",color=INK,size=base-.8),
      legend.key.height=grid::unit(.42,"cm"),
      plot.tag=element_text(face="bold",color=INK,size=base+7),
      plot.margin=margin(11,12,11,12)
    )
}
fmt_q <- function(x) {
  ifelse(!is.finite(x),"NA",ifelse(x<0.001,formatC(x,format="e",digits=1),
                                        sprintf("%.3f",x)))
}
save_as <- function(p,name,w,h) {
  dest <- file.path(OUT,name)
  tryCatch(
    ggsave(paste0(dest,".pdf"),p,width=w,height=h,units="in",
           bg="white",device="pdf",useDingbats=FALSE),
    error=function(e) stop("PDF export failed for ",name,": ",conditionMessage(e))
  )
  if (requireNamespace("ragg",quietly=TRUE)) {
    ggsave(paste0(dest,".png"),p,width=w,height=h,units="in",
           dpi=600,bg="white",device=ragg::agg_png)
  } else {
    ggsave(paste0(dest,".png"),p,width=w,height=h,units="in",
           dpi=600,bg="white")
  }
  message("Saved: ",name," (.pdf, .png)")
}

# ---- G: TWO SIX-DONOR PAIRED BOXPLOTS ---------------------------------------
Gdata$xbase <- ifelse(Gdata$Condition=="CIT",1,2)
donors <- sort(unique(Gdata$Donor))
# Same horizontal offset for both samples of a donor (connections stay paired).
set.seed(126L)
donor_dx <- stats::setNames(stats::runif(length(donors),-.060,.060),donors)
Gdata$xj <- Gdata$xbase + unname(donor_dx[Gdata$Donor])
metric_plot <- function(metric,label,stat_idx,tag=NULL) {
  d <- Gdata
  d$value <- as.numeric(d[[metric]])
  value_min <- min(d$value);value_max <- max(d$value)
  ran <- value_max-value_min
  if(ran<=0) ran <- max(1,abs(value_max)*.2)
  line_y <- value_max + .092*ran
  text_y <- value_max + .163*ran
  ann <- sprintf("Paired P = %s   |   BH-FDR = %s",
                 fmt_q(Gstat$PairedWilcoxon_P[stat_idx]),
                 fmt_q(Gstat$BH_FDR_two_metrics[stat_idx]))
  ggplot(d,aes(y=value)) +
    geom_line(aes(x=xj,group=Donor),linewidth=.64,alpha=.45,color="#8C939B") +
    geom_boxplot(aes(x=xbase,group=Condition,fill=Condition),
                 width=.39,outlier.shape=NA,alpha=.35,
                 color=INK,linewidth=.98) +
    geom_point(aes(x=xj,color=Condition),size=2.65,alpha=.98) +
    annotate("segment",x=1,xend=2,y=line_y,yend=line_y,
             linewidth=.83,color=INK) +
    annotate("segment",x=c(1,2),xend=c(1,2),
             y=line_y-.029*ran,yend=line_y,linewidth=.83,color=INK) +
    annotate("text",x=1.5,y=text_y,label=ann,fontface="bold",
             color=INK,size=3.10) +
    scale_fill_manual(values=c(CIT=CIT,HR=HR),guide="none") +
    scale_color_manual(values=c(CIT=CIT,HR=HR),guide="none") +
    scale_x_continuous(breaks=c(1,2),labels=c("CIT","HR"),limits=c(.57,2.43)) +
    coord_cartesian(ylim=c(value_min-.10*ran,value_max+.24*ran),clip="off") +
    labs(title=label,x=NULL,y=label,tag=tag) +
    font_theme(11.2) + theme(panel.grid.major.x=element_blank())
}
pG1 <- metric_plot("MeanSERPINE1","Mean SERPINE1 expression",1,"G")
pG2 <- metric_plot("PositivePercent","SERPINE1-positive ECs (%)",2)
pG <- wrap_plots(pG1,pG2,nrow=1) +
  plot_annotation(title="Endothelial SERPINE1 response",
    subtitle="GSE220797  |  CIT vs HR  |  six matched donors",
    theme=theme(plot.title=element_text(face="bold",hjust=.5,color=INK,size=15),
                plot.subtitle=element_text(face="bold",hjust=.5,color=INK,size=10.5)))
save_as(pG,"Figure3G_SERPINE1_paired_PREMIUM",11.4,5.3)

# ---- H: HEATMAP + EDGE-R EFFECT PANEL --------------------------------------
ordered_groups <- c("Injury / senescence","Inflammatory activation",
                    "Permeability / remodeling","Barrier / homeostasis")
ordered_genes <- unique(Hdata$Gene)
if(!all(ordered_groups %in% unique(Hdata$GeneGroup)))
  stop("Unexpected Figure H gene group labels")
Hdata$GeneGroup <- factor(Hdata$GeneGroup, levels=ordered_groups)
Hstat$GeneGroup <- factor(Hstat$GeneGroup, levels=ordered_groups)
Hdata$Gene <- factor(Hdata$Gene,levels=rev(ordered_genes))
Hstat$Gene <- factor(Hstat$Gene,levels=rev(ordered_genes))
Hdata$Sample <- factor(Hdata$Sample,levels=expected)
Hdata$Zplot <- pmax(-2,pmin(2,Hdata$Z)) # visualization cap ONLY

pHheat <- ggplot(Hdata,aes(x=Sample,y=Gene,fill=Zplot))+
  geom_tile(color=WHITE,linewidth=.57,width=.99,height=.99)+
  geom_vline(xintercept=6.5,color=INK,linewidth=1.03)+
  facet_grid(GeneGroup~.,scales="free_y",space="free_y",switch="y")+
  scale_fill_gradient2(low="#377EBA",mid="white",high="#E35855",
                       midpoint=0,limits=c(-2,2),name="Row z-score")+
  labs(title="Paired endothelial expression",x=NULL,y=NULL,tag="H")+
  font_theme(9.1)+
  theme(panel.grid=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text.x=element_text(angle=58,hjust=1,vjust=1,size=8.3,face="bold"),
        strip.placement="outside",strip.text.y.left=element_text(angle=0,size=7.7),
        strip.background=element_rect(fill="#F5F6F7",color="#DCE1E5"),
        panel.spacing.y=grid::unit(.13,"cm"),
        legend.position="bottom",legend.key.width=grid::unit(.78,"cm"))

Hstat$Direction <- factor(ifelse(Hstat$logFC >= 0,"HR higher","CIT higher"),
                          levels=c("CIT higher","HR higher"))
# BH-FDR is for ALL genes in the original edgeR DEG table, per existing table.
Hstat$FDRLabel <- paste0("q=",fmt_q(Hstat$FDR))
minfc <- min(c(Hstat$logFC,0));maxfc <- max(c(Hstat$logFC,0))
spanfc <- max(maxfc-minfc,.1)
pHforest <- ggplot(Hstat,aes(y=Gene))+
  geom_vline(xintercept=0,color="#8C939B",linetype="dashed",linewidth=.7)+
  geom_segment(aes(x=0,xend=logFC,yend=Gene,color=Direction),linewidth=1.3,alpha=.82)+
  geom_point(aes(x=logFC,fill=Direction),shape=21,size=3.6,color=INK,stroke=.35)+
  geom_text(aes(x=maxfc+spanfc*.16,label=FDRLabel),color=INK,
            hjust=0,fontface="bold",size=2.65)+
  facet_grid(GeneGroup~.,scales="free_y",space="free_y")+
  scale_color_manual(values=c("CIT higher"=CIT,"HR higher"=HR),guide="none")+
  scale_fill_manual(values=c("CIT higher"=CIT,"HR higher"=HR),guide="none")+
  scale_x_continuous(limits=c(minfc-spanfc*.17,maxfc+spanfc*.83))+
  labs(title="Paired model effect",subtitle="edgeR QLF; q = genome-wide BH-FDR",
       x="log2 FC (HR - CIT)",y=NULL)+
  font_theme(9.1)+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),panel.spacing.y=grid::unit(.13,"cm"),
        strip.text.y=element_blank(),strip.background=element_blank(),
        panel.grid.major.x=element_line(color=GRID,linewidth=.3))
pH <- wrap_plots(pHheat,pHforest,nrow=1,widths=c(1.65,1.18))
save_as(pH,"Figure3H_EC_heatmap_and_paired_logFC_PREMIUM",15.6,8.5)

# ---- I + J: DIRECTIONAL GSEA LOLLIPOPS -------------------------------------
# Display prespecified top 6 up-/down-NES terms ranked by EXISTING padj.
# The full tested-table remains untouched and accompanies the plots.
pretty_term <- function(z,prefix,limit=36L) {
  z <- sub(paste0("^",prefix),"",z)
  z <- gsub("_"," ",z)
  z <- tools::toTitleCase(tolower(z))
  vapply(z,function(q) paste(strwrap(q,width=limit),collapse="\n"),character(1))
}
select_terms <- function(d,n=6L){
 d <- d[is.finite(d$NES)&is.finite(d$padj)&is.finite(d$size),,drop=FALSE]
 a <- head(d[which(d$NES>0)[order(d$padj[d$NES>0],-abs(d$NES[d$NES>0]))],],n)
 b <- head(d[which(d$NES<0)[order(d$padj[d$NES<0],-abs(d$NES[d$NES<0]))],],n)
 z <- rbind(a,b)
 z <- z[order(z$NES),,drop=FALSE]
 if(!nrow(z)) stop("No testable terms in GSEA table")
 z
}
make_enrichment <- function(d,tag,name,prefix){
 z <- select_terms(d,6L)
 z$Display <- pretty_term(z$pathway,prefix,33L)
 # A unique factor even if text wraps to identical display names.
 z$ItemID <- make.unique(z$Display)
 z$ItemID <- factor(z$ItemID,levels=z$ItemID)
 key <- stats::setNames(z$Display,as.character(z$ItemID))
 z$Direction <- ifelse(z$NES>0,"HR enriched","CIT enriched")
 z$Direction <- factor(z$Direction,levels=c("CIT enriched","HR enriched"))
 z$BHLabel <- paste0("q=",fmt_q(z$padj))
 xrange <- max(abs(z$NES))
 label_x <- xrange+.28
 p <- ggplot(z,aes(y=ItemID))+
   geom_vline(xintercept=0,linetype="dashed",color="#89939C",linewidth=.75)+
   geom_segment(aes(x=0,xend=NES,yend=ItemID,color=Direction),linewidth=1.18)+
   geom_point(aes(x=NES,fill=Direction,size=size),shape=21,
              color=INK,stroke=.44)+
   geom_text(aes(x=label_x,label=BHLabel),color=INK,fontface="bold",
             hjust=0,size=2.78)+
   scale_color_manual(values=c("CIT enriched"=CIT,"HR enriched"=HR),guide="none")+
   scale_fill_manual(values=c("CIT enriched"=CIT,"HR enriched"=HR),guide="none")+
   scale_size_continuous(range=c(2.9,6.2),name="Gene-set size")+
   scale_y_discrete(labels=function(k) unname(key[k]))+
   scale_x_continuous(limits=c(-xrange-.45,xrange+1.1),
                      breaks=pretty(c(-xrange,xrange),n=5))+
   labs(title=name,
        subtitle="Positive NES: HR  |  Negative NES: CIT",
        x="Normalized enrichment score (NES)",y=NULL,tag=tag,
        caption="Top 6 per direction by existing BH-FDR; all tested terms retained in input CSV")+
   font_theme(10.2)+
   theme(axis.text.y=element_text(size=8.4,lineheight=.96,color=INK),
         axis.text.x=element_text(size=9.1),
         panel.grid.major.y=element_blank(),
         panel.grid.major.x=element_line(color=GRID,linewidth=.35),
         legend.position="bottom",legend.direction="horizontal")
 p
}
pI <- make_enrichment(Idata,"I","KEGG pathway landscape","KEGG_")
pJ <- make_enrichment(Jdata,"J","GO biological-process landscape","GOBP_")
save_as(pI,"Figure3I_paired_KEGG_lollipop_PREMIUM",10.6,7.85)
save_as(pJ,"Figure3J_paired_GO_BP_lollipop_PREMIUM",10.6,7.85)

# ---- K: FOUR FOCUSED PATHWAYS, HONEST ABOUT UNTESTABLE TERMS ----------------
# The supplied paired analysis has two pathways with Found=FALSE. Without
# corresponding tested gene sets / running ES data, FOUR old GSEA curves
# CANNOT be recreated truthfully. Show a carefully labeled summary instead.
Kdata$Found <- as.logical(Kdata$Found)
Kdata$Pathway <- factor(Kdata$Pathway,levels=rev(Kdata$Pathway))
found <- Kdata[!is.na(Kdata$Found)&Kdata$Found,,drop=FALSE]
missing <- Kdata[is.na(Kdata$Found)|!Kdata$Found,,drop=FALSE]
if(nrow(found)) {
  if(any(!is.finite(found$NES))||any(!is.finite(found$FDR)))
    stop("K indicates a pathway is present, but NES or FDR missing")
  found$Direction <- factor(ifelse(found$NES>=0,"HR enriched","CIT enriched"),
                            levels=c("CIT enriched","HR enriched"))
  found$annotation <- paste0("NES ",sprintf("%.2f",found$NES),
                             "   |   q ",fmt_q(found$FDR))
}
xl <- max(2.2,if(nrow(found)) max(abs(found$NES))+.4 else 2.2)
pK <- ggplot()+
  geom_vline(xintercept=0,color="#9FA7AF",linewidth=.69,linetype="dashed")+
  geom_segment(data=found,aes(x=0,xend=NES,y=Pathway,yend=Pathway,
                              color=Direction),linewidth=1.7,alpha=.95)+
  geom_point(data=found,aes(x=NES,y=Pathway,fill=Direction),
             shape=21,color=INK,stroke=.5,size=4.2)+
  geom_text(data=found,aes(x=xl+.20,y=Pathway,label=annotation),
            hjust=0,size=3.3,color=INK,fontface="bold")+
  geom_text(data=missing,aes(x=.11,y=Pathway),
            label="Not available in tested KEGG set",
            hjust=0,size=3.1,color="#515A63",fontface="bold")+
  scale_fill_manual(values=c("CIT enriched"=CIT,"HR enriched"=HR),guide="none")+
  scale_color_manual(values=c("CIT enriched"=CIT,"HR enriched"=HR),guide="none")+
  scale_x_continuous(limits=c(-.45,xl+1.8),breaks=c(0,1,2))+
  labs(title="Focused endothelial KEGG pathways",
       subtitle="Paired pseudobulk  |  matched to the tested MSigDB KEGG collection",
       x="Normalized enrichment score (NES)",y=NULL,tag="K",
       caption="Unmatched pathways have no inferred NES/P/FDR; q is the original KEGG-wide BH-FDR")+
  font_theme(11.5)+
  theme(axis.text.y=element_text(size=10.2,color=INK),
        panel.grid.major.x=element_line(color=GRID,linewidth=.38))
save_as(pK,"Figure3K_focused_pathways_status_PREMIUM",14.4,4.7)

# ---- COMPOSITE: G–K only, preserve original A–F -----------------------------
# Authors should replace old G-K panels in the original Figure 3 layout using
# independent vector PDFs, without rerunning/reformatting A-F.
pTOP <- wrap_plots(pG,pH,widths=c(1.08,1.38),nrow=1)
pMID <- wrap_plots(pI,pJ,widths=c(1,1),nrow=1)
pCOMBINED <- wrap_plots(pTOP,pMID,pK,ncol=1,heights=c(1.10,1.16,.76))+
  plot_annotation(theme=theme(plot.background=element_rect(fill="white",color=NA)))
save_as(pCOMBINED,"Figure3_G-K_5panels_PREMIUM",23.5,22.8)

# Summary aids: sourced values, not recalculated inference.
write.csv(data.frame(Metric=Gstat$Metric,
                     paired_P=Gstat$PairedWilcoxon_P,
                     BH_FDR_two_metrics=Gstat$BH_FDR_two_metrics),
          file.path(OUT,"G_original_stats_no_recalculation.csv"),row.names=FALSE)
write.csv(Kdata,file.path(OUT,"K_pathway_availability_original_values.csv"),
          row.names=FALSE)
writeLines(c(
 "Figure 3 G-K: cosmetic replot ONLY (CSV input)",
 paste0("Data tables: ",TABLE_DIR),
 "CIT blue #2878D4; HR red #E3322B; all text bold and black.",
 "G paired P/FDR: exactly from Table 08, two metrics, n=6 donors.",
 "H expression Z values: from Table 09; edgeR HR-CIT logFC/FDR: Table 10.",
 "I/J full-GSEA values: Tables 11/12; top 6 per enrichment direction by FDR.",
 "K: two pathways unavailable in supplied tested KEGG set; never invent ES curves.",
 "A-F intentionally untouched.",
 capture.output(sessionInfo())
 ),file.path(OUT,"RENDER_PROVENANCE.txt"))
message("SUCCESS: output folder = ",normalizePath(OUT,winslash="/"))



}

# =============================================================================
# STAGE figure3_kegg_curves
# =============================================================================
PIPELINE$figure3_kegg_curves <- function() {
# scripts/05_targeted_kegg_curves.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

# =============================================================================
# Figure 3K — original Figure4 script's KEGG GSEA calculation and plotting
# GSE220797 | CIT1–HR1, ..., CIT6–HR6 | 6 DONOR-PAIRED pseudobulk samples
# R 4.3.2; fgsea 1.28.0; KEGGREST 1.42.0; ggplot2 3.5.1; patchwork 1.2.0
# 
# BASED ON original source Figure4_FINAL_WITH_CELL_COMPOSITION_NO_L(4).R
# Uses EXACT original sequence:
#   KEGGREST::keggGet -> AnnotationDbi::mapIds -> fgsea::fgsea ->
#   fgsea::plotEnrichment -> patchwork::wrap_plots(ncol=4)
# Only input replaced: 6-donor-paired edgeR pseudobulk SIGNED QLF rank.
# DOES NOT rerun edgeR, Seurat, QC, Harmony, clustering, or panels A-J.
# =============================================================================
options(stringsAsFactors=FALSE)
set.seed(1234)

# 1. FILES --------------------------------------------------------------------
# Keep Tables/ alongside this script, as in the supplied ZIP.
# Or set RANK_FILE to the exact saved paired-rank CSV path.
RANK_FILE <- repo_input("figure3_signed_rank")
RANK_BASENAME <- "06_paired_signed_rank_for_GSEA.csv"
if (!nzchar(RANK_FILE)) {
  roots <- unique(c(getwd(), file.path(getwd(), "Tables")))
  matches <- file.path(roots,RANK_BASENAME)
  matches <- matches[file.exists(matches)]
  if (length(matches)) RANK_FILE <- matches[[1L]]
  else {
    message("Select the ORIGINAL PAIRED rank CSV: ",RANK_BASENAME)
    RANK_FILE <- file.choose()
  }
}
if (!file.exists(RANK_FILE)) stop("Paired rank file not found: ",RANK_FILE)
RANK_FILE <- normalizePath(RANK_FILE,winslash="/",mustWork=TRUE)
if (basename(RANK_FILE) != RANK_BASENAME)
  stop("Incorrect rank table. Must be ",RANK_BASENAME,
       " (6 donor-paired edgeR QLF rank), not an unpaired file.")
OUTDIR <- file.path(dirname(RANK_FILE),"FIG3K_OFFICIAL_KEGG_ORIGINAL_STYLE")
GSEA_DIR <- file.path(OUTDIR,"GSEA_individual")
TABLE_DIR <- file.path(OUTDIR,"Tables")
for (d in c(OUTDIR,GSEA_DIR,TABLE_DIR)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

for (pkg in c("KEGGREST","AnnotationDbi","org.Hs.eg.db","fgsea","ggplot2","patchwork")) {
  if (!requireNamespace(pkg,quietly=TRUE))
    stop("Required package missing: ",pkg,". Install it in the existing R 4.3.2 environment.")
}

# 2. MATCH ORIGINAL PATHWAY IDS/COLORS --------------------------------------
TARGET_KEGG <- data.frame(
  KEGG_ID = c("hsa04218", "hsa04066", "hsa04512", "hsa04520"),
  KEGG_Name = c("Cellular senescence", "HIF-1 signaling pathway",
                "ECM-receptor interaction", "Adherens junction"),
  stringsAsFactors = FALSE
)
PATHWAY_COLORS <- c(
  "hsa04218" = "#E67E22", "hsa04066" = "#2F6FB0",
  "hsa04512" = "#D63384", "hsa04520" = "#169C9C"
)
GSEA_MIN_SIZE <- 5L
GSEA_MAX_SIZE <- 500L

# 3. RESTORE EXACT PAIRED PSEUDOBULK RANK -----------------------------------
pair_rank <- utils::read.csv(RANK_FILE,check.names=FALSE,
                             stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
required <- c("Gene","Signed_sqrt_QLF")
if (!all(required %in% names(pair_rank)))
  stop("Missing expected paired-rank columns. Found: ",paste(names(pair_rank),collapse=", "))
pair_rank$Gene <- trimws(as.character(pair_rank$Gene))
pair_rank$Signed_sqrt_QLF <- suppressWarnings(as.numeric(pair_rank$Signed_sqrt_QLF))
pair_rank <- pair_rank[!is.na(pair_rank$Gene) & nzchar(pair_rank$Gene) &
                         is.finite(pair_rank$Signed_sqrt_QLF),,drop=FALSE]
if (anyDuplicated(pair_rank$Gene))
  stop("Duplicate gene symbols in paired rank; inspect the original pseudobulk export.")
rank_score <- pair_rank$Signed_sqrt_QLF
names(rank_score) <- pair_rank$Gene
rank_score <- sort(rank_score,decreasing=TRUE)
# Preserve the tiny deterministic tie-break adjustment from the ORIGINAL GSEA code.
rank_score <- rank_score + seq(0,1e-10,length.out=length(rank_score))
if (length(rank_score)<500L) stop("Paired rank has too few genes: ",length(rank_score))
if (any(diff(rank_score)>0)) stop("Paired rank sorting error")
# positive = HR > CIT; negative = CIT > HR. This is a signed sqrt(QLF F).
message("Loaded ",length(rank_score)," paired pseudobulk ranked genes.")
message("Ranked left HR-enriched; ranked right CIT-enriched.")

# 4. ORIGINAL OFFICIAL KEGG FUNCTION -----------------------------------------
# This is the same get_official_kegg_symbols mechanism as the source script.
# The only addition is a frozen local cache for reproducibility.
get_official_kegg_symbols <- function(pathway_id, pathway_name) {
  entry <- tryCatch(KEGGREST::keggGet(pathway_id)[[1L]],error=function(e) NULL)
  if (!is.null(entry)) {
    gene_field <- entry$GENE
    if (!is.null(gene_field) && length(gene_field)>0L) {
      if (!is.null(names(gene_field)) && any(nzchar(names(gene_field))))
        entrez_ids <- names(gene_field)
      else
        entrez_ids <- gene_field[seq(1L,length(gene_field),by=2L)]
      entrez_ids <- sub("^hsa:","",as.character(entrez_ids))
      entrez_ids <- sub("[^0-9].*$","",entrez_ids)
      entrez_ids <- unique(entrez_ids[grepl("^[0-9]+$",entrez_ids)])
      if (!length(entrez_ids)) stop("KEGG GENE ID parsing failed: ",pathway_id)
      symbol_map <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys=entrez_ids,keytype="ENTREZID",column="SYMBOL",multiVals="first"
      )
      symbols <- unique(unname(as.character(symbol_map)))
      symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
      if (length(symbols)>=GSEA_MIN_SIZE) return(symbols)
    }
  }
  stop("Could not retrieve official KEGG ",pathway_id," (",pathway_name,
       "). Internet connection to KEGG REST or an existing membership cache is required.")
}

CACHE_FILE <- file.path(dirname(RANK_FILE),"Figure3K_official_KEGG_membership_cache.rds")
if (file.exists(CACHE_FILE)) {
  cache <- readRDS(CACHE_FILE)
  if (!is.list(cache) || !is.list(cache$pathways) ||
      !identical(names(cache$pathways),TARGET_KEGG$KEGG_ID))
    stop("KEGG membership cache does not match the four predefined official pathways")
  official_pathway_genes <- cache$pathways
  message("Reusing official KEGG membership cache: ",CACHE_FILE)
} else {
  message("Fetching four official KEGG pathways using ORIGINAL code mechanism...")
  official_pathway_genes <- stats::setNames(
    lapply(seq_len(nrow(TARGET_KEGG)),function(i) {
      get_official_kegg_symbols(TARGET_KEGG$KEGG_ID[i],TARGET_KEGG$KEGG_Name[i])
    }),TARGET_KEGG$KEGG_ID)
  cache <- list(pathways=official_pathway_genes,
                downloaded_UTC=format(Sys.time(),tz="UTC",usetz=TRUE),
                KEGGREST_version=as.character(utils::packageVersion("KEGGREST")),
                org_Hs_eg_db_version=as.character(utils::packageVersion("org.Hs.eg.db")))
  saveRDS(cache,CACHE_FILE)
  message("Cached official KEGG membership at: ",CACHE_FILE)
}

raw_members <- do.call(rbind,lapply(names(official_pathway_genes),function(id){
  data.frame(KEGG_ID=id,KEGG_Name=TARGET_KEGG$KEGG_Name[match(id,TARGET_KEGG$KEGG_ID)],
             Gene=official_pathway_genes[[id]],stringsAsFactors=FALSE)
}))
write.csv(raw_members,file.path(TABLE_DIR,"01_official_KEGG_gene_membership.csv"),row.names=FALSE)

# 5. ORIGINAL FGSEA ANALYSIS --------------------------------------------------
official_pathway_genes_ranked <- lapply(official_pathway_genes,function(genes_i) {
  intersect(genes_i,names(rank_score))
})
pathway_gene_summary <- data.frame(
  KEGG_ID=TARGET_KEGG$KEGG_ID,
  KEGG_Name=TARGET_KEGG$KEGG_Name,
  OfficialGeneNumber=lengths(official_pathway_genes[TARGET_KEGG$KEGG_ID]),
  MatchedRankGeneNumber=lengths(official_pathway_genes_ranked[TARGET_KEGG$KEGG_ID]),
  stringsAsFactors=FALSE
)
write.csv(pathway_gene_summary,
          file.path(TABLE_DIR,"02_Four_official_KEGG_gene_summary.csv"),row.names=FALSE)
valid_pathways <- official_pathway_genes_ranked[
  lengths(official_pathway_genes_ranked)>=GSEA_MIN_SIZE &
  lengths(official_pathway_genes_ranked)<=GSEA_MAX_SIZE]
if (length(valid_pathways)!=4L)
  stop("Cannot perform ORIGINAL four-panel GSEA: only ",length(valid_pathways),
       "/4 KEGG gene sets have 5-500 matched genes. Check 02_Four_official_KEGG_gene_summary.csv")

# EXACT same function and relevant arguments as original Figure4 script.
set.seed(1234)
gsea_result <- fgsea::fgsea(
  pathways=valid_pathways,stats=rank_score,
  minSize=GSEA_MIN_SIZE,maxSize=GSEA_MAX_SIZE,eps=0,nproc=1
)
gsea_result <- as.data.frame(gsea_result)
if (nrow(gsea_result)!=4L) stop("fgsea did not return all four official pathways")
if ("leadingEdge" %in% names(gsea_result)) {
  gsea_result$leadingEdge <- vapply(gsea_result$leadingEdge,
                                   paste,collapse=";",FUN.VALUE=character(1))
}
# As in source, padj from this FOUR-pathway fgsea test is BH-FDR across FOUR.
gsea_result$KEGG_ID <- as.character(gsea_result$pathway)
ix <- match(TARGET_KEGG$KEGG_ID,gsea_result$KEGG_ID)
if (anyNA(ix)) stop("One or more predefined pathways were not tested")
gsea_result <- gsea_result[ix,,drop=FALSE]
gsea_result$KEGG_Name <- TARGET_KEGG$KEGG_Name
gsea_result$FDR_4_pathways <- gsea_result$padj
gsea_result$EnrichedToward <- ifelse(gsea_result$NES>0,"HR","CIT")
gsea_result$MatchedGenes <- pathway_gene_summary$MatchedRankGeneNumber
if (anyNA(gsea_result$NES) || anyNA(gsea_result$pval))
  stop("Nonfinite GSEA statistics; inspect original paired-rank distribution")
write.csv(gsea_result,file.path(TABLE_DIR,"03_Four_official_KEGG_GSEA_paired_FDR4.csv"),row.names=FALSE)

# 6. ORIGINAL plotEnrichment() CLASSIC RUNNING-ES CURVES ---------------------
fmt <- function(x){
  if(!is.finite(x)) return("NA")
  if(x<0.001) format(x,digits=2,scientific=TRUE) else sprintf("%.3f",x)
}
original_style <- function(){
  ggplot2::theme_classic(base_size=10,base_family="sans")+
    ggplot2::theme(
      text=ggplot2::element_text(color="black",face="bold"),
      plot.title=ggplot2::element_text(color="black",face="bold",hjust=.5,size=12,
                                       margin=ggplot2::margin(b=6)),
      plot.subtitle=ggplot2::element_text(color="black",face="bold",hjust=.5,size=9,
                                          margin=ggplot2::margin(b=5)),
      axis.title=ggplot2::element_text(color="black",face="bold",size=9.5),
      axis.text=ggplot2::element_text(color="black",face="bold",size=9),
      axis.line=ggplot2::element_line(color="black",linewidth=.67),
      axis.ticks=ggplot2::element_line(color="black",linewidth=.55),
      panel.border=ggplot2::element_rect(fill=NA,color="black",linewidth=.55),
      panel.grid=ggplot2::element_blank(),
      plot.margin=ggplot2::margin(8,8,7,8))
}

plot_list <- vector("list",nrow(TARGET_KEGG))
for (i in seq_len(nrow(TARGET_KEGG))) {
  pathway_id <- TARGET_KEGG$KEGG_ID[i]
  result_i <- gsea_result[gsea_result$KEGG_ID==pathway_id,,drop=FALSE]
  # Exactly as in user's source: plotEnrichment includes running ES + gene ticks.
  plot_i <- fgsea::plotEnrichment(
    pathway=valid_pathways[[pathway_id]],stats=rank_score
  )+
    ggplot2::geom_hline(yintercept=0,linetype="dashed",linewidth=.35,
                        color="#89959E")+
    ggplot2::labs(
      title=TARGET_KEGG$KEGG_Name[i],
      subtitle=paste0("NES ",sprintf("%.2f",result_i$NES[1]),
                      "   P ",fmt(result_i$pval[1]),
                      "   FDR(4) ",fmt(result_i$FDR_4_pathways[1]),
                      "   ",result_i$EnrichedToward[1]),
      x="HR-enriched  ←  pseudobulk rank  →  CIT-enriched",
      y="Running enrichment score") + original_style()
  # Original four pathway line colors; ALL TEXT BLACK per user request.
  if (length(plot_i$layers)>=1L) {
    plot_i$layers[[1L]]$aes_params$colour <- unname(PATHWAY_COLORS[pathway_id])
    plot_i$layers[[1L]]$aes_params$linewidth <- 1.10
    plot_i$layers[[1L]]$aes_params$size <- NULL
  }
  plot_list[[i]] <- plot_i
  subdir_name <- paste0(pathway_id,"_",gsub("[^A-Za-z0-9]+","_",TARGET_KEGG$KEGG_Name[i]))
  ggplot2::ggsave(file.path(GSEA_DIR,paste0(subdir_name,".pdf")),plot_i,
                  width=6.7,height=4.4,bg="white",useDingbats=FALSE)
  if (requireNamespace("ragg",quietly=TRUE)) {
    ggplot2::ggsave(file.path(GSEA_DIR,paste0(subdir_name,".png")),plot_i,
                    width=6.7,height=4.4,dpi=600,bg="white",device=ragg::agg_png)
  } else {
    ggplot2::ggsave(file.path(GSEA_DIR,paste0(subdir_name,".png")),plot_i,
                    width=6.7,height=4.4,dpi=600,bg="white")
  }
}
names(plot_list) <- TARGET_KEGG$KEGG_ID

# Same 1 x 4 ordering as source: hsa04218, hsa04066, hsa04512, hsa04520.
figure3K_original <- patchwork::wrap_plots(plot_list,ncol=4) +
  patchwork::plot_annotation(
    title="Focused pathway GSEA (donor-paired endothelial pseudobulk)",
    theme=ggplot2::theme(
      plot.title=ggplot2::element_text(face="bold",color="black",hjust=.5,size=16)))

# Optional 2x2 version: identical plotting, larger typography in manuscript.
figure3K_2x2 <- patchwork::wrap_plots(plot_list,ncol=2) +
  patchwork::plot_annotation(
    title="Focused pathway GSEA (donor-paired endothelial pseudobulk)",
    theme=ggplot2::theme(
      plot.title=ggplot2::element_text(face="bold",color="black",hjust=.5,size=16)))

save_combined <- function(p,stem,width,height){
  ggplot2::ggsave(file.path(OUTDIR,paste0(stem,".pdf")),p,width=width,height=height,
                  bg="white",useDingbats=FALSE)
  if (requireNamespace("ragg",quietly=TRUE))
    ggplot2::ggsave(file.path(OUTDIR,paste0(stem,".png")),p,width=width,height=height,
                    dpi=600,bg="white",device=ragg::agg_png)
  else
    ggplot2::ggsave(file.path(OUTDIR,paste0(stem,".png")),p,width=width,height=height,
                    dpi=600,bg="white")
}
save_combined(figure3K_original,"Figure3K_four_official_KEGG_original_1x4",20,5.3)
save_combined(figure3K_2x2,"Figure3K_four_official_KEGG_alternative_2x2",13.8,9.4)

writeLines(c(
  "Based on original Figure4_FINAL_WITH_CELL_COMPOSITION_NO_L(4).R GSEA code.",
  paste("Paired 6-donor pseudobulk rank:",RANK_FILE),
  paste("Official KEGG membership cache:",CACHE_FILE),
  "Original workflow: KEGGREST::keggGet -> AnnotationDbi::mapIds ->",
  "fgsea::fgsea(4 selected pathways) -> fgsea::plotEnrichment.",
  "Positive signed rank = HR up vs CIT, displayed at LEFT.",
  "BH FDR(4) is for the four tested a priori KEGG pathways ONLY.",
  "This differs from the whole-KEGG FDR family in Fig 3I if using a different database.",
  "No rerun of edgeR, Seurat, Harmony or other figure panels."
),file.path(TABLE_DIR,"04_analysis_methods_and_statistical_scope.txt"))
writeLines(capture.output(sessionInfo()),file.path(TABLE_DIR,"05_sessionInfo.txt"))
print(gsea_result[,c("KEGG_ID","KEGG_Name","NES","pval","FDR_4_pathways","EnrichedToward","MatchedGenes")])
message("Figure 3K official KEGG GSEA finished: ",normalizePath(OUTDIR,winslash="/"))



}

# =============================================================================
# STAGE figure4_states
# =============================================================================
PIPELINE$figure4_states <- function() {
# scripts/06_endothelial_states.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

# =============================================================================
# Figure 4 reviewer revision: selective recomputation H, J, K, L, M
# Previously labelled 'Figure 5' in original source; manuscript Figure 4.
# Original source: Figure5_Step3_PREMIUM_A_to_M_FINAL_FIXED.R
# Source-declared package versions: R 4.3.2; Seurat 4.4.0; edgeR 4.0.16;
# ggplot2 3.5.1; patchwork 1.2.0; msigdbr 7.5.1.
# Re-execution on the final source object is required before submission.
# DO NOT rerun clustering, Harmony integration, marker detection or A-G / I.
# =============================================================================
options(stringsAsFactors = FALSE)
set.seed(1234)

# ---------------------------- CONFIGURATION ---------------------------------
RDS_FILE <- repo_input("endothelial_bundle")# full path if needed
OUTDIR <- repo_output("Figure4")
MIN_CELLS_PER_STATE_SAMPLE <- 10L   # H and M eligibility
MIN_CELLS_PER_HIGH_LOW <- 5L        # K: matched high/low per sample+state
MIN_PAIRS_PER_STATE <- 3L           # small samples: n reported, not assumed n=6
KEY_MIN_CELLS_PER_SAMPLE <- 10L     # J state pseudobulk
KEY_STATE <- 'Stress-responsive capillary EC'
SERPINE1_HIGH_QUANTILE <- 0.75      # fixed threshold from CIT positive cells
MAKE_MATCHED_UMAP <- TRUE           # may take time; only new UMAP, no Harmony rerun
RUN_PCA_CLUSTER_SENSITIVITY <- FALSE# NEW independent QC reclustering from saved PCA
RUN_KEY_LEAVE_ONE_DONOR_OUT <- FALSE# fit 6 small pseudobulk edgeR leave-one-pair-out
USE_EXISTING_PCA_AND_HARMONY <- TRUE
N_FEATURES_UMAP <- 20L
SEED <- 1234L
COLOR_CIT <- '#2166AC'; COLOR_HR <- '#D73027'
# ----------------------------------------------------------------------------

required <- c('Seurat','SeuratObject','Matrix','edgeR','ggplot2','patchwork',
              'msigdbr','clusterProfiler')
missing_pkg <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_pkg)) stop('Missing packages: ',paste(missing_pkg,collapse=', '),
  '. Install from CRAN/Bioconductor in your existing R environment first.')
if(!file.exists(RDS_FILE)) stop('Cannot find RDS_FILE: ',RDS_FILE)
dir.create(OUTDIR,showWarnings=FALSE,recursive=TRUE)
D <- file.path(OUTDIR,'Data_and_Audits'); P <- file.path(OUTDIR,'Revised_Panels')
S <- file.path(OUTDIR,'Supplementary_Stability')
for(dd in c(D,P,S)) dir.create(dd,showWarnings=FALSE,recursive=TRUE)
log_file <- file.path(OUTDIR,'EXECUTION_LOG.txt')
log_message <- function(...) { msg<-paste0(...);cat(msg,'\n');cat(msg,'\n',file=log_file,append=TRUE) }
log_message('Figure 4 selective reviewer revision; date: ',as.character(Sys.time()))
log_message('Original object: ',normalizePath(RDS_FILE,winslash='/',mustWork=TRUE))

obj <- readRDS(RDS_FILE)
if(!is.list(obj)||!inherits(obj$endo_clean,'Seurat')||
   !all(LETTERS[1:13] %in% names(obj$panels)))
  stop('RDS must be PREMIUM: obj$endo_clean is Seurat and obj$panels has A-M.')
endo <- obj$endo_clean
SeuratObject::DefaultAssay(endo) <- 'RNA'
md <- endo@meta.data
needed_md <- c('orig.ident','tissue_type','Endothelial_Annotation')
if(!all(needed_md %in% names(md))) stop('Missing metadata: ',
                                        paste(setdiff(needed_md,names(md)),collapse=', '))
if(!'SERPINE1' %in% rownames(endo[['RNA']])) stop('SERPINE1 missing in RNA assay')

# FIXED, VERIFIED 6-PAIR DESIGN. Never infer donor from row ordering.
md$Sample <- toupper(trimws(as.character(md$orig.ident)))
md$Condition <- toupper(trimws(as.character(md$tissue_type)))
expected <- c(paste0('CIT',1:6),paste0('HR',1:6))
if(!setequal(unique(md$Sample),expected))
 stop('Sample identifiers must be CIT1 ... CIT6 / HR1 ... HR6. Found: ',
       paste(sort(unique(md$Sample)),collapse=', '))
if(any(!md$Condition %in% c('CIT','HR')) ||
   any((grepl('^CIT',md$Sample)&md$Condition!='CIT') |
       (grepl('^HR',md$Sample)&md$Condition!='HR')))
 stop('tissue_type inconsistent with actual CIT/HR sample labels.')
md$Donor <- factor(sub('^(CIT|HR)','',md$Sample),levels=as.character(1:6))
md$Condition <- factor(md$Condition,levels=c('CIT','HR'))
md$State <- as.character(md$Endothelial_Annotation)
if(anyNA(md$Donor)||anyNA(md$State)) stop('Donor or state has missing values')
state_order <- unique(md$State)
if ('clean_cluster' %in% names(md)) {
  state_cluster_table <- table(md$clean_cluster, md$State)
} else {
  state_cluster_table <- NULL
}
meta_samp <- unique(md[,c('Sample','Condition','Donor'),drop=FALSE])
meta_samp <- meta_samp[match(expected,meta_samp$Sample),,drop=FALSE]
stopifnot(nrow(meta_samp)==12L, !anyNA(meta_samp$Sample),
          all(table(meta_samp$Donor)==2L))
write.csv(meta_samp,file.path(D,'01_verified_6pairs_CIT_HR_samples.csv'),row.names=FALSE)
log_message('Matched sample pairs verified: CIT1-HR1 ... CIT6-HR6; endothelial cells=',nrow(md))

# Original plot objects ONLY from RDS, no reanalysis of A-G and I
revised <- obj$panels
if(!all(vapply(revised,function(x) inherits(x,'ggplot')||inherits(x,'patchwork'),logical(1))))
 warning('Some original panels may be wrapped; composite export will be attempted.')
STATE_COLORS <- c(
  'Homeostatic capillary EC'='#43B1D0',
  'Transitional capillary EC'='#8ED3C7',
  'Stress-responsive capillary EC'='#E9B000',
  'Homeostatic arterial EC'='#6420A8',
  'Remodeling arterial EC'='#AC77CC',
  'Adhesion-activated venous EC'='#F09857',
  'Remodeling venous EC'='#C94E50',
  'Inflammatory venous EC'='#B90020',
  'Homeostatic lymphatic EC'='#079D7C',
  'Activated lymphatic EC'='#5A9E24')
extra <- setdiff(state_order,names(STATE_COLORS))
if(length(extra)) {newcols<-grDevices::hcl.colors(length(extra),'Dark 3');
                  names(newcols)<-extra;STATE_COLORS<-c(STATE_COLORS,newcols)}
STATE_COLORS <- STATE_COLORS[unique(c(state_order,names(STATE_COLORS)))]
COLORS <- c(CIT=COLOR_CIT,HR=COLOR_HR)

# -------- SEURAT V4 RAW + NORMALIZED MATRICES (NEVER integrated counts) ------
count_mat <- SeuratObject::GetAssayData(endo,assay='RNA',slot='counts')
expr_mat <- SeuratObject::GetAssayData(endo,assay='RNA',slot='data')
if(!identical(colnames(count_mat),rownames(md)) ||
   !identical(colnames(expr_mat),rownames(md)))
 stop('Seurat RNA counts/data columns do not match meta.data row order.')
ser <- as.numeric(expr_mat['SERPINE1',]); names(ser)<-colnames(expr_mat)
md$SERPINE1 <- ser[rownames(md)]
if(anyNA(md$SERPINE1)) stop('Missing SERPINE1 expression after metadata matching.')

black_theme <- function(size=10){
 ggplot2::theme_classic(base_size=size,base_family='Arial')+
 ggplot2::theme(text=ggplot2::element_text(color='black',face='bold'),
   plot.title=ggplot2::element_text(color='black',face='bold',size=size+1.5,hjust=.5),
   plot.tag=ggplot2::element_text(color='black',face='bold',size=size+4),
   plot.subtitle=ggplot2::element_text(color='black',face='bold',size=size-1,hjust=.5),
   axis.title=ggplot2::element_text(color='black',face='bold'),
   axis.text=ggplot2::element_text(color='black',face='bold'),
   axis.line=ggplot2::element_line(color='black',linewidth=.7),
   axis.ticks=ggplot2::element_line(color='black',linewidth=.5),
   strip.text=ggplot2::element_text(face='bold',color='black',size=size-1),
   strip.background=ggplot2::element_rect(fill='white',color='black',linewidth=.45),
   legend.title=ggplot2::element_text(face='bold',color='black'),
   legend.text=ggplot2::element_text(face='bold',color='black'),
   panel.border=ggplot2::element_rect(fill=NA,color='black',linewidth=.5),
   plot.margin=ggplot2::margin(7,10,7,9))}
fmt_p <- function(x) if(!is.finite(x)) 'NA' else if(x<.001) format(x,digits=2,scientific=TRUE) else sprintf('%.3f',x)
export <- function(p,name,dir=P,w=7,h=5){
  # separate files, original A-G and I untouched
  dev <- if(capabilities('cairo')) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(file.path(dir,paste0(name,'.pdf')),plot=p,width=w,height=h,
      units='in',device=dev,bg='white',limitsize=FALSE)
  if(requireNamespace('ragg',quietly=TRUE)) {
    ggplot2::ggsave(file.path(dir,paste0(name,'.png')),plot=p,width=w,height=h,
                    units='in',dpi=600,device=ragg::agg_png,bg='white',limitsize=FALSE)
  } else ggplot2::ggsave(file.path(dir,paste0(name,'.png')),plot=p,width=w,height=h,
                    units='in',dpi=600,bg='white',limitsize=FALSE)
}

# ---------------- H: 6 donor-paired sample-level expression -----------------
num <- aggregate(rep(1L,nrow(md)),
     by=list(Sample=md$Sample,State=md$State),FUN=sum)
names(num)[3]<-'Ncells'
agg_mean <- aggregate(md$SERPINE1,
     by=list(Sample=md$Sample,State=md$State),FUN=mean)
names(agg_mean)[3]<-'MeanSERPINE1'
sd <- merge(num,agg_mean,by=c('Sample','State'))
sd$Condition <- factor(ifelse(grepl('^CIT',sd$Sample),'CIT','HR'),levels=c('CIT','HR'))
sd$Donor <- sub('^(CIT|HR)','',sd$Sample)
write.csv(sd,file.path(D,'H_all_sample_state_mean_SERPINE1.csv'),row.names=FALSE)
valid_h <- sd[sd$Ncells >= MIN_CELLS_PER_STATE_SAMPLE,,drop=FALSE]

match_donor <- function(df,val){
  cit <- df[df$Condition=='CIT',c('Donor',val),drop=FALSE]
  hr <- df[df$Condition=='HR',c('Donor',val),drop=FALSE]
  names(cit)[2]<-'CIT';names(hr)[2]<-'HR'
  merge(cit,hr,by='Donor',sort=TRUE)
}
h_list <- lapply(sort(unique(md$State)),function(st){
  d<-valid_h[valid_h$State==st,,drop=FALSE]
  w<-match_donor(d,'MeanSERPINE1');n<-nrow(w)
  if (n >= MIN_PAIRS_PER_STATE && any(w$HR != w$CIT)) {
    p <- suppressWarnings(stats::wilcox.test(w$HR, w$CIT, paired=TRUE, exact=FALSE)$p.value)
  } else if (n >= MIN_PAIRS_PER_STATE) {
    p <- 1
  } else {
    p <- NA_real_
  }
  data.frame(State=st,Pairs=n,MeanPairedDelta=if(n)mean(w$HR-w$CIT) else NA_real_,
             MedianPairedDelta=if(n)median(w$HR-w$CIT) else NA_real_,P=p)
})
h_stats<-do.call(rbind,h_list);h_stats$FDR<-stats::p.adjust(h_stats$P,method='BH')
h_stats<-h_stats[order(-abs(h_stats$MeanPairedDelta),na.last=TRUE),]
write.csv(h_stats,file.path(D,'H_PAIRED_Wilcoxon_BH_10_states.csv'),row.names=FALSE)
show_h <- head(h_stats$State[is.finite(h_stats$P)],4L)
if(!length(show_h)) stop('No EC state has >=',MIN_PAIRS_PER_STATE,' valid complete donor pairs for H.')
hplot <- valid_h[valid_h$State %in% show_h,,drop=FALSE]
complete_h <- do.call(rbind,lapply(show_h,function(st){
  w<-match_donor(valid_h[valid_h$State==st,,drop=FALSE],'MeanSERPINE1')
  valid_h[valid_h$State==st & valid_h$Donor %in% w$Donor,,drop=FALSE]
}))
complete_h$State<-factor(complete_h$State,levels=show_h)
labels_h<-h_stats[match(show_h,h_stats$State),];labels_h$State<-factor(labels_h$State,levels=show_h)
labels_h$label<-paste0('BH-FDR=',vapply(labels_h$FDR,fmt_p,character(1)),' | n=',labels_h$Pairs)
labels_h$Y<-vapply(show_h,function(st){z<-complete_h$MeanSERPINE1[complete_h$State==st];max(z)+.17*max(diff(range(z)),.2)},numeric(1))
set.seed(SEED)
# Fixed horizontal offset per donor to maintain correct connections
jit<-stats::setNames(seq(-.075,.075,length.out=6),as.character(1:6))
complete_h$x<-ifelse(complete_h$Condition=='CIT',1,2)+jit[complete_h$Donor]
pH<-ggplot2::ggplot(complete_h,ggplot2::aes(y=MeanSERPINE1))+
 ggplot2::geom_line(ggplot2::aes(x=x,group=Donor),color='#8C8C8C',linewidth=.48,alpha=.7)+
 ggplot2::geom_boxplot(ggplot2::aes(x=ifelse(Condition=='CIT',1,2),group=Condition,fill=Condition),
                    width=.38,outlier.shape=NA,alpha=.18,linewidth=.7,color='black')+
 ggplot2::geom_point(ggplot2::aes(x=x,color=Condition),size=2.1)+
 ggplot2::geom_text(data=labels_h,ggplot2::aes(x=1.5,y=Y,label=label),
                    inherit.aes=FALSE,fontface='bold',color='black',size=2.8)+
 ggplot2::facet_wrap(~State,nrow=1,scales='free_y')+
 ggplot2::scale_color_manual(values=COLORS,guide='none')+
 ggplot2::scale_fill_manual(values=COLORS,guide='none')+
 ggplot2::scale_x_continuous(breaks=1:2,labels=c('CIT','HR'),limits=c(.65,2.35))+
 ggplot2::labs(title='Paired endothelial SERPINE1 expression',
             subtitle='Sample-level means; complete CIT/HR donor pairs only',
             x=NULL,y='Mean normalized SERPINE1')+black_theme(9)
pH <- pH + ggplot2::labs(tag='H')
revised$H<-pH;export(pH,'Figure4H_REVISED_paired_SERPINE1',w=11,h=4.5)
log_message('H: paired SERPINE1 comparisons for ',length(show_h),' displayed states; BH across tested states.')

# ---------------- J: paired pseudobulk key-state Hallmark GSEA --------------
key_md <- md[md$State==KEY_STATE,,drop=FALSE]
if(!nrow(key_md)) stop('KEY_STATE missing: ',KEY_STATE)
key_counts_by_sample<-table(factor(key_md$Sample,levels=expected))
write.csv(data.frame(Sample=expected,Cells=as.integer(key_counts_by_sample)),
          file.path(D,'J_key_state_cells_by_sample.csv'),row.names=FALSE)
valid_samples<-names(key_counts_by_sample)[key_counts_by_sample>=KEY_MIN_CELLS_PER_SAMPLE]
complete_pairs<-as.character(1:6)[vapply(as.character(1:6),function(id)
   all(c(paste0('CIT',id),paste0('HR',id)) %in% valid_samples),logical(1))]
if(length(complete_pairs)<3L) stop('J: insufficient complete key-state donor pairs after cell-count QC: ',
    length(complete_pairs),'. See J_key_state_cells_by_sample.csv.')
key_samples<-c(paste0('CIT',complete_pairs),paste0('HR',complete_pairs))
keep_cells<-rownames(key_md)[key_md$Sample %in% key_samples]
key_cell_samp<-factor(md[keep_cells,'Sample'],levels=key_samples)
A<-Matrix::sparse.model.matrix(~0+key_cell_samp);colnames(A)<-key_samples
key_pb<-count_mat[,keep_cells,drop=FALSE] %*% A
key_meta<-data.frame(Sample=key_samples,
   Donor=factor(rep(complete_pairs,2)),
   Condition=factor(rep(c('CIT','HR'),each=length(complete_pairs)),levels=c('CIT','HR')))
design_j<-stats::model.matrix(~Donor+Condition,data=key_meta)
if(qr(design_j)$rank!=ncol(design_j)) stop('J paired design rank deficient')
y<-edgeR::DGEList(counts=as.matrix(key_pb),samples=key_meta)
keep_gene<-edgeR::filterByExpr(y,design=design_j)
if(sum(keep_gene)<100) stop('J key-state pseudobulk retained <100 genes')
y<-y[keep_gene,,keep.lib.sizes=FALSE];y<-edgeR::calcNormFactors(y)
y<-edgeR::estimateDisp(y,design_j,robust=TRUE)
fit<-edgeR::glmQLFit(y,design_j,robust=TRUE)
test<-edgeR::glmQLFTest(fit,coef='ConditionHR')
key_de<-edgeR::topTags(test,n=Inf,sort.by='PValue')$table
key_de$Gene<-rownames(key_de);key_de$FDR<-p.adjust(key_de$PValue,'BH')
write.csv(key_de,file.path(D,'J_PAIRED_key_state_edgeR_DE.csv'),row.names=FALSE)
log_message('J: ',length(complete_pairs),' complete donors, ',nrow(key_de),' genes; HR-CIT paired edgeR.')

# ----------------------------- Hallmark GSEA --------------------------------
msig<-msigdbr::msigdbr(species='Homo sapiens',category='H')
if(!all(c('gs_name','gene_symbol')%in%names(msig)))stop('msigdbr H output lacks gene set fields')
t2g<-unique(data.frame(term=as.character(msig$gs_name),
                       gene=as.character(msig$gene_symbol)))
t2g<-t2g[!is.na(t2g$gene)&nzchar(t2g$gene),,drop=FALSE]
t2n<-data.frame(term=unique(t2g$term),stringsAsFactors=FALSE)
t2n$name<-tools::toTitleCase(tolower(gsub('_',' ',sub('^HALLMARK_','',t2n$term))))
run_hallmark <- function(de,exclude=character()){
 z<-de[is.finite(de$PValue)&is.finite(de$logFC)&!de$Gene%in%exclude,
       c('Gene','logFC','PValue'),drop=FALSE]
 z<-z[!duplicated(z$Gene),]
 rank<-sign(z$logFC)*(-log10(pmax(z$PValue,1e-300)))
 names(rank)<-z$Gene;rank<-sort(rank,decreasing=TRUE)
 if(length(rank)<100L) stop('GSEA rank contains fewer than 100 genes')
 out<-clusterProfiler::GSEA(geneList=rank,TERM2GENE=t2g,
       TERM2NAME=t2n,exponent=1,minGSSize=10,maxGSSize=500,
       pvalueCutoff=1,pAdjustMethod='BH',verbose=FALSE,
       seed=TRUE,by='fgsea')
 list(table=as.data.frame(out),rank=rank)
}
prep_gsea <- function(dat,max_each=5L){
 if(!nrow(dat)) return(dat)
 dat<-dat[is.finite(dat$NES)&is.finite(dat$p.adjust),,drop=FALSE]
 pos<-dat[dat$NES>0,,drop=FALSE];neg<-dat[dat$NES<0,,drop=FALSE]
 pos<-head(pos[order(pos$p.adjust,-abs(pos$NES)),,drop=FALSE],max_each)
 neg<-head(neg[order(neg$p.adjust,-abs(neg$NES)),,drop=FALSE],max_each)
 r<-rbind(neg,pos);r<-r[order(r$NES),,drop=FALSE]
 r$Description<-factor(r$Description,levels=unique(r$Description))
 r
}
plot_gsea <- function(data,title,negative,positive){
 if(!nrow(data)){
   return(ggplot2::ggplot()+ggplot2::annotate('text',x=0,y=0,label='No eligible Hallmark terms',
           fontface='bold',color='black',size=4)+ggplot2::theme_void()+
           ggplot2::labs(title=title)+black_theme(10))
 }
 data$Direction<-factor(ifelse(data$NES>0,positive,negative),levels=c(negative,positive))
 data$DotSize<-pmin(-log10(pmax(data$p.adjust,1e-12)),12)
 ggplot2::ggplot(data,ggplot2::aes(x=NES,y=Description))+
   ggplot2::geom_vline(xintercept=0,color='black',linetype='dashed',linewidth=.45)+
   ggplot2::geom_segment(ggplot2::aes(x=0,xend=NES,yend=Description,color=Direction),
                        linewidth=.9)+
   ggplot2::geom_point(ggplot2::aes(color=Direction,size=DotSize),alpha=.95)+
   ggplot2::scale_color_manual(values=stats::setNames(c(COLOR_CIT,COLOR_HR),c(negative,positive)))+
   ggplot2::scale_size_continuous(range=c(2.1,5.6),name='-log10(FDR)')+
   ggplot2::labs(title=title,x='Normalized enrichment score (NES)',y=NULL,color='Enriched toward')+
   black_theme(9)+ggplot2::theme(axis.text.y=ggplot2::element_text(size=7.3,color='black',face='bold'))
}
gsea_j<-run_hallmark(key_de)
write.csv(gsea_j$table,file.path(D,'J_PAIRED_key_state_Hallmark_GSEA.csv'),row.names=FALSE)
write.csv(data.frame(Gene=names(gsea_j$rank),Rank=gsea_j$rank),
          file.path(D,'J_PAIRED_key_state_Hallmark_rank.csv'),row.names=FALSE)
pJ<-plot_gsea(prep_gsea(gsea_j$table),
             'Stress-responsive capillary EC | paired Hallmark',
             'CIT','HR')
pJ <- pJ + ggplot2::labs(tag='J')
revised$J<-pJ;export(pJ,'Figure4J_REVISED_paired_Hallmark',w=7,h=5.1)

# --------------- K + L: balanced, sample-paired High vs Low ------------------
# One fixed threshold derived ONLY from CIT-positive expression, no ties.method
# = "first". Low = undetected (zero). High = >75th percentile CIT-positive.
# Primary: within each (sample,state), select equal number high and low cells;
# this balances state composition between High and Low *within every sample*.
# Biol. unit = DONOR (six donor-matched pairs). After matching cells within
# each sample+state, combine the 2 timepoints within donor for each High/Low arm.
# Model: ~ Donor + SERPINE1_State. This avoids treating 12 paired libraries
# as 12 independent biological donors.
cit_pos<-md$SERPINE1[md$Condition=='CIT' & md$SERPINE1>0]
if(length(cit_pos)<100L) stop('CIT-positive SERPINE1 count too low to define fixed threshold')
cutoff<-as.numeric(stats::quantile(cit_pos,probs=SERPINE1_HIGH_QUANTILE,na.rm=TRUE))
if(!is.finite(cutoff)||cutoff<=0) stop('Fixed positive SERPINE1 threshold is invalid')
log_message('Fixed High threshold (CIT positive q',SERPINE1_HIGH_QUANTILE,') = ',format(cutoff,digits=5))
md$HighFixed<-md$SERPINE1>cutoff
md$HL<-ifelse(md$SERPINE1>cutoff,'High',
              ifelse(md$SERPINE1==0,'Low','Middle'))
write.csv(data.frame(Method='CIT-positive reference quantile',Quantile=SERPINE1_HIGH_QUANTILE,
                     Threshold=cutoff,Low='expression == 0',High='expression > threshold'),
          file.path(D,'KL_fixed_high_low_definition.csv'),row.names=FALSE)
# Eligible sample/state strata, balance High/Low cell numbers
md$Stratum<-paste(md$Sample,md$State,sep='::')
byblock<-split(rownames(md),md$Stratum)
chosen<-list();audit<-vector('list',length(byblock));k<-0L
set.seed(SEED)
for(nm in names(byblock)){
 k<-k+1L; cells<-byblock[[nm]];d<-md[cells,,drop=FALSE]
 high<-rownames(d)[d$HL=='High'];low<-rownames(d)[d$HL=='Low']
 nmatch<-min(length(high),length(low))
 audit[[k]]<-data.frame(Sample=d$Sample[1],State=d$State[1],
       N_High=length(high),N_Low=length(low),N_PairMatched=nmatch,
       Included=(nmatch>=MIN_CELLS_PER_HIGH_LOW))
 if(nmatch<MIN_CELLS_PER_HIGH_LOW) next
 chosen[[nm]]<-c(sample(high,nmatch),sample(low,nmatch))
}
kl_audit<-do.call(rbind,audit)
write.csv(kl_audit,file.path(D,'KL_sample_state_high_low_balance_audit.csv'),row.names=FALSE)
kl_cells<-unlist(chosen,use.names=FALSE)
if(length(kl_cells)<2*MIN_CELLS_PER_HIGH_LOW*6) stop(
  'K: Too few balanced high/low cells. See KL_sample_state_high_low_balance_audit.csv')
kl_md<-md[kl_cells,,drop=FALSE]
kl_n<-as.data.frame(table(Sample=kl_md$Sample,HL=kl_md$HL))
kl_n<-kl_n[kl_n$HL%in%c('High','Low'),]
kl_mat<-reshape(kl_n,idvar='Sample',timevar='HL',direction='wide')
valid_kl_samples<-as.character(kl_mat$Sample[kl_mat$Freq.High>=MIN_CELLS_PER_HIGH_LOW &
                                 kl_mat$Freq.Low>=MIN_CELLS_PER_HIGH_LOW])
# Require both CIT and HR represented within a donor, before pooling counts.
valid_kl_donors<-as.character(1:6)[vapply(as.character(1:6),function(id)
  all(c(paste0('CIT',id),paste0('HR',id)) %in% valid_kl_samples),logical(1))]
if(length(valid_kl_donors)<3L)
 stop('K: fewer than 3 complete donors with BOTH CIT and HR high/low arms. See audit table.')
valid_kl_samples<-c(paste0('CIT',valid_kl_donors),paste0('HR',valid_kl_donors))
kl_md<-kl_md[kl_md$Sample%in%valid_kl_samples,,drop=FALSE]
log_message('K balanced high/low: ',length(valid_kl_donors),' DONOR pairs, ',
    length(valid_kl_samples),' samples; state composition matched per sample')
# Pair biological replicates at donor level (CIT and HR both contribute
# equally in cell counts to the high and low arms, by sample+state balancing).
kl_groups<-as.vector(rbind(paste0(valid_kl_donors,'__Low'),
                            paste0(valid_kl_donors,'__High')))
group_fac<-factor(paste(kl_md$Donor,kl_md$HL,sep='__'),levels=kl_groups)
B<-Matrix::sparse.model.matrix(~0+group_fac);colnames(B)<-kl_groups
kl_pb<-count_mat[,rownames(kl_md),drop=FALSE] %*% B
kl_meta<-data.frame(Group=kl_groups,
   Donor=factor(rep(valid_kl_donors,each=2),levels=valid_kl_donors),
   SERPINE1_State=factor(rep(c('Low','High'),times=length(valid_kl_donors)),
                         levels=c('Low','High')))
design_kl<-stats::model.matrix(~Donor+SERPINE1_State,data=kl_meta)
if(qr(design_kl)$rank!=ncol(design_kl)) stop('K paired DONOR high/low model is rank deficient')
k_y<-edgeR::DGEList(counts=as.matrix(kl_pb),samples=kl_meta)
kl_keep<-edgeR::filterByExpr(k_y,design=design_kl)
if(sum(kl_keep)<100) stop('K after filterByExpr: fewer than 100 genes')
k_y<-k_y[kl_keep,,keep.lib.sizes=FALSE];k_y<-edgeR::calcNormFactors(k_y)
k_y<-edgeR::estimateDisp(k_y,design_kl,robust=TRUE)
k_fit<-edgeR::glmQLFit(k_y,design_kl,robust=TRUE)
k_res<-edgeR::glmQLFTest(k_fit,coef='SERPINE1_StateHigh')
k_de<-edgeR::topTags(k_res,n=Inf,sort.by='PValue')$table
k_de$Gene<-rownames(k_de);k_de$FDR<-p.adjust(k_de$PValue,'BH')
write.csv(k_de,file.path(D,'K_matched_state_balanced_donor_paired_edgeR_DE_FULL.csv'),row.names=FALSE)
k_without<-k_de[k_de$Gene!='SERPINE1',,drop=FALSE]
write.csv(k_without,file.path(D,'K_matched_state_balanced_edgeR_DE_excl_SERPINE1.csv'),row.names=FALSE)
# Report complete gene-wise BH denominator; exclude SERPINE1 from display only
k_without$Class<-ifelse(k_without$FDR<.05 & k_without$logFC>=.5,'Higher in High',
                  ifelse(k_without$FDR<.05 & k_without$logFC<=-.5,'Higher in Low',
                         'Not significant'))
k_without$Y<-pmin(-log10(pmax(k_without$FDR,1e-300)),18)
label_k<-rbind(head(k_without[k_without$Class=='Higher in High',
  ][order(k_without$FDR[k_without$Class=='Higher in High']),,drop=FALSE],5L),
  head(k_without[k_without$Class=='Higher in Low',
  ][order(k_without$FDR[k_without$Class=='Higher in Low']),,drop=FALSE],5L))
pK<-ggplot2::ggplot(k_without,ggplot2::aes(x=logFC,y=Y))+
 ggplot2::geom_point(ggplot2::aes(color=Class),size=.82,alpha=.72)+
 ggplot2::geom_vline(xintercept=c(-.5,.5),linetype='dashed',color='#8C8C8C',linewidth=.4)+
 ggplot2::geom_hline(yintercept=-log10(.05),linetype='dashed',color='#8C8C8C',linewidth=.4)+
 ggplot2::scale_color_manual(values=c('Higher in High'=COLOR_HR,'Higher in Low'=COLOR_CIT,
                                       'Not significant'='#C2C7CC'))+
 ggplot2::labs(title='SERPINE1-high vs low | composition-matched',
       subtitle='Donor-paired pseudobulk; FDR shown; SERPINE1 excluded',
       x='log2 fold change (High - Low)',y='-log10(BH-FDR)',color='Differential expression')+
 black_theme(9)
if(nrow(label_k) && requireNamespace('ggrepel',quietly=TRUE))
 pK<-pK+ggrepel::geom_text_repel(data=label_k,ggplot2::aes(label=Gene),
       color='black',fontface='bold',size=2.6,max.overlaps=Inf,seed=SEED)
pK <- pK + ggplot2::labs(tag='K')
revised$K<-pK;export(pK,'Figure4K_REVISED_matched_high_low_volcano',w=7,h=5.3)

# ---------------- L: Hallmark GSEA, excluding SERPINE1 ----------------------
gsea_l<-run_hallmark(k_de,exclude='SERPINE1')
write.csv(gsea_l$table,file.path(D,'L_balanced_high_low_Hallmark_GSEA_excl_SERPINE1.csv'),row.names=FALSE)
write.csv(data.frame(Gene=names(gsea_l$rank),Rank=gsea_l$rank),
          file.path(D,'L_balanced_high_low_GSEA_rank_excl_SERPINE1.csv'),row.names=FALSE)
pL<-plot_gsea(prep_gsea(gsea_l$table),
   'SERPINE1-associated Hallmark programs | state-balanced',
   'SERPINE1-low','SERPINE1-high')
pL <- pL + ggplot2::labs(tag='L')
revised$L<-pL;export(pL,'Figure4L_REVISED_composition_matched_Hallmark',w=8.5,h=5)

# -------- M: donor-paired HR/CIT OR, with a full auditable sample/state grid --
# Critically, retain ALL sample x state combinations and the zero-count bins.
# No reshape() on a table with missing observations; do not silently lose donors.
# Per-donor log odds ratio = logit(HR high fraction) - logit(CIT high fraction),
# with Haldane-Anscombe 0.5 continuity correction at each time point.
# Mean within-donor logOR = pooled PAIR-AWARE estimate. A t CI on donor-level
# logORs needs >= MIN_PAIRS_PER_STATE complete donors and finite nonzero variance.
# Its estimand is a paired donor-level geometric mean odds ratio, NOT
# a fixed-cell-effect logit GLM and NOT an across-state enrichment OR.
# Fractions and exact donor eligibility are exported even if inference impossible.

m_ct <- as.data.frame(table(
  Sample=factor(md$Sample,levels=expected),
  State=factor(md$State,levels=state_order),
  High=factor(md$HighFixed,levels=c(FALSE,TRUE))))
colnames(m_ct) <- c('Sample','State','HighFlag','Cells')
m_ct$Sample <- as.character(m_ct$Sample)
m_ct$State <- as.character(m_ct$State)
m_ct$HighFlag <- as.character(m_ct$HighFlag)
m_ct$Cells <- as.integer(m_ct$Cells)

m_total <- stats::aggregate(Cells ~ Sample + State,data=m_ct,FUN=sum)
m_high <- m_ct[m_ct$HighFlag=='TRUE',c('Sample','State','Cells'),drop=FALSE]
names(m_total)[3] <- 'Total'
names(m_high)[3] <- 'High'
ct_wide <- merge(m_total,m_high,by=c('Sample','State'),all.x=TRUE,sort=FALSE)
ct_wide$High[is.na(ct_wide$High)] <- 0L
ct_wide$Low <- ct_wide$Total-ct_wide$High
ct_wide$Donor <- sub('^(CIT|HR)','',ct_wide$Sample)
ct_wide$Condition <- factor(ifelse(grepl('^CIT',ct_wide$Sample),'CIT','HR'),
                              levels=c('CIT','HR'))
ct_wide$HighFraction <- ifelse(ct_wide$Total>0,
                                ct_wide$High/ct_wide$Total,NA_real_)
ct_wide$EligibleSample <- ct_wide$Total>=MIN_CELLS_PER_STATE_SAMPLE
# 12 samples x each state; including combinations with ZERO state cells.
stopifnot(nrow(ct_wide)==12L*length(state_order),
          sum(ct_wide$Total)==nrow(md),
          sum(ct_wide$High)==sum(md$HighFixed),
          !anyDuplicated(paste(ct_wide$Sample,ct_wide$State)))
write.csv(ct_wide,file.path(D,'M_sample_state_fixed_high_fractions.csv'),row.names=FALSE)

m_donor_rows <- list()
m_result_rows <- list()
for(st in state_order){
  dd <- ct_wide[ct_wide$State==st,,drop=FALSE]
  donor_delta <- numeric(0)
  valid_fr_cit <- valid_fr_hr <- numeric(0)
  for(id in as.character(1:6)){
    c_row <- dd[dd$Sample==paste0('CIT',id),,drop=FALSE]
    h_row <- dd[dd$Sample==paste0('HR',id),,drop=FALSE]
    if(nrow(c_row)!=1L||nrow(h_row)!=1L)
      stop('M internal paired grid error at state ',st,' donor ',id)
    eligible <- c_row$EligibleSample && h_row$EligibleSample
    delta <- if(eligible){
      base::log((h_row$High+.5)/(h_row$Low+.5)) -
        base::log((c_row$High+.5)/(c_row$Low+.5))
    } else NA_real_
    if(length(delta)!=1L || (eligible && !is.finite(delta)))
      stop('M invalid paired log-OR: state=',st,', donor=',id)
    m_donor_rows[[length(m_donor_rows)+1L]] <- data.frame(
      State=st,Donor=id,CIT_n=c_row$Total,HR_n=h_row$Total,
      CIT_high=c_row$High,HR_high=h_row$High,
      CIT_fraction=c_row$HighFraction,HR_fraction=h_row$HighFraction,
      PairEligible=eligible,Paired_Log_OR=delta,
      ExclusionReason=if(eligible)'Included' else if(c_row$Total<MIN_CELLS_PER_STATE_SAMPLE &&
                                    h_row$Total<MIN_CELLS_PER_STATE_SAMPLE) 'Both sides too few cells'
                          else if(c_row$Total<MIN_CELLS_PER_STATE_SAMPLE) 'CIT too few cells'
                          else 'HR too few cells', stringsAsFactors=FALSE)
    if(eligible){
      donor_delta <- c(donor_delta,delta)
      valid_fr_cit <- c(valid_fr_cit,c_row$HighFraction)
      valid_fr_hr <- c(valid_fr_hr,h_row$HighFraction)
    }
  }
  n <- length(donor_delta)
  ratio <- if(n) exp(mean(donor_delta)) else NA_real_
  lo <- hi <- p <- wp <- NA_real_
  reason <- if(n<MIN_PAIRS_PER_STATE) 'Too few donor pairs for inference' else 'OK'
  if(n>=MIN_PAIRS_PER_STATE){
    s <- stats::sd(donor_delta)
    if(is.finite(s) && s>1e-12){
      se <- s/sqrt(n)
      tc <- stats::qt(.975,df=n-1L)
      lo <- exp(mean(donor_delta)-tc*se)
      hi <- exp(mean(donor_delta)+tc*se)
      p <- 2*stats::pt(-abs(mean(donor_delta)/se),df=n-1L)
    } else {
      reason <- 'Paired logOR has zero variance: CI and t-P not estimable'
    }
    if(all(abs(valid_fr_hr-valid_fr_cit)<1e-12)) {
      wp <- 1
    } else {
      wp <- suppressWarnings(stats::wilcox.test(valid_fr_hr,valid_fr_cit,
                                paired=TRUE,exact=FALSE)$p.value)
    }
  }
  m_result_rows[[length(m_result_rows)+1L]] <- data.frame(
     State=st,Pairs=n,OR=ratio,CI_low=lo,CI_high=hi,P=p,
     WilcoxonP=wp,Reason=reason,stringsAsFactors=FALSE)
}
m_donor_audit <- do.call(rbind,m_donor_rows)
write.csv(m_donor_audit,file.path(D,'M_DONOR_PAIR_AUDIT_ALL_STATES.csv'),row.names=FALSE)
m_stats <- do.call(rbind,m_result_rows)
m_stats$FDR <- stats::p.adjust(m_stats$P,method='BH')
m_stats$WilcoxonFDR <- stats::p.adjust(m_stats$WilcoxonP,method='BH')
write.csv(m_stats,file.path(D,'M_PAIRED_HR_vs_CIT_OR_fixed_high_fraction_BH.csv'),
          row.names=FALSE)
log_message('M diagnostic: states=',nrow(m_stats),'; pair counts: ',
    paste(m_stats$State,m_stats$Pairs,sep='=',collapse=' | '))

mm <- m_stats[is.finite(m_stats$OR)&m_stats$OR>0,,drop=FALSE]
if(nrow(mm)){
  mm$State <- factor(mm$State,levels=rev(state_order))
  mm$Label <- ifelse(is.finite(mm$FDR),
    paste0('n=',mm$Pairs,' | FDR=',vapply(mm$FDR,fmt_p,character(1))),
    paste0('n=',mm$Pairs,' | exploratory (CI unavailable)'))
  mm_ci <- mm[is.finite(mm$CI_low)&is.finite(mm$CI_high)&
                mm$CI_low>0&mm$CI_high>0,,drop=FALSE]
  pM <- ggplot2::ggplot(mm,ggplot2::aes(x=OR,y=State,color=State))+
    ggplot2::geom_vline(xintercept=1,linetype='dashed',color='black',linewidth=.55)+
    ggplot2::geom_segment(data=mm_ci,
      ggplot2::aes(x=CI_low,xend=CI_high,yend=State),linewidth=.9)+
    ggplot2::geom_point(shape=21,fill='white',size=2.6,stroke=.9)+
    ggplot2::geom_text(ggplot2::aes(label=Label),hjust=-.05,
      size=2.6,color='black',fontface='bold',show.legend=FALSE)+
    ggplot2::scale_x_log10(expand=ggplot2::expansion(mult=c(.08,.65)))+
    ggplot2::scale_color_manual(values=STATE_COLORS,guide='none')+
    ggplot2::labs(title='Paired SERPINE1-high fraction changes across EC states',
      subtitle=paste0('High: CIT-positive reference threshold=',signif(cutoff,3),
          '  |  paired donor OR; no CI when insufficient replication'),
      x='Geometric mean of donor-paired OR (HR/CIT; 95% CI where estimable)',
      y=NULL)+black_theme(8.3)+ggplot2::coord_cartesian(clip='off')+
    ggplot2::theme(axis.text.y=ggplot2::element_text(size=7.1),
      plot.margin=ggplot2::margin(9,30,9,9))
} else {
  # No estimable OR is a real result. Preserve the full figure without
  # fabricating values or altering prespecified minimum sample criteria.
  pM <- ggplot2::ggplot() +
    ggplot2::annotate('text',x=0,y=0,label=paste0(
       'No endothelial states have eligible paired ORs.\n',
       'Cell threshold: ',MIN_CELLS_PER_STATE_SAMPLE,
       ' per condition; min donors: ',MIN_PAIRS_PER_STATE,
       '\nSee M_DONOR_PAIR_AUDIT_ALL_STATES.csv'),
       color='black',fontface='bold',size=3.6)+
    ggplot2::theme_void(base_family='Arial')+
    ggplot2::labs(title='Paired SERPINE1-high fraction: not estimable')+
    ggplot2::theme(plot.title=ggplot2::element_text(
       color='black',face='bold',hjust=.5))
  log_message('M WARNING: zero states have eligible paired odds ratio; see diagnostic CSV')
}
pM <- pM + ggplot2::labs(tag='M')
revised$M <- pM
export(pM,'Figure4M_REVISED_PAIRED_HR_vs_CIT_OR',w=8.1,h=5.4)

# --------------------------- STABILITY AUDIT --------------------------------
# Sample contributions and pairwise fractions across all endothelial states
sample_state<-as.data.frame(table(Sample=factor(md$Sample,levels=expected),
                                   State=factor(md$State,levels=state_order)))
names(sample_state)[3]<-'Cells'
sample_state$Donor<-sub('^(CIT|HR)','',as.character(sample_state$Sample))
sample_state$Condition<-factor(ifelse(grepl('^CIT',sample_state$Sample),'CIT','HR'),
                                levels=c('CIT','HR'))
per_sample<-aggregate(Cells~Sample,sample_state,sum);names(per_sample)[2]<-'AllEC'
per_state<-aggregate(Cells~State,sample_state,sum);names(per_state)[2]<-'AllStateCells'
sample_state<-merge(sample_state,per_sample,by='Sample')
sample_state<-merge(sample_state,per_state,by='State')
sample_state$WithinSamplePct<-100*sample_state$Cells/sample_state$AllEC
sample_state$ShareOfStatePct<-100*sample_state$Cells/sample_state$AllStateCells
write.csv(sample_state,file.path(S,'S1_EC_state_contribution_12_samples.csv'),row.names=FALSE)
coverage<-do.call(rbind,lapply(state_order,function(st){
  dat<-sample_state[sample_state$State==st,,drop=FALSE]
  eligible<-dat$Sample[dat$Cells>=MIN_CELLS_PER_STATE_SAMPLE]
  dn<-as.character(1:6)
  paired<-dn[vapply(dn,function(n) all(c(paste0('CIT',n),paste0('HR',n)) %in% eligible),logical(1))]
  data.frame(State=st,TotalCells=sum(dat$Cells),MaxSingleSamplePercent=max(dat$ShareOfStatePct),
     SamplesWith10Plus=sum(dat$Cells>=MIN_CELLS_PER_STATE_SAMPLE),
     ValidCompleteDonorPairs=length(paired),
     Donors=paste(paired,collapse=';'))
}))
write.csv(coverage,file.path(S,'S2_subcluster_sample_coverage_audit.csv'),row.names=FALSE)
if(!is.null(state_cluster_table)) write.csv(as.data.frame.matrix(state_cluster_table),
        file.path(S,'S3_cluster_to_annotation_contingency.csv'))
sp<-sample_state
sp$State<-factor(sp$State,levels=state_order)
sp$Sample<-factor(sp$Sample,levels=expected)
p_comp<-ggplot2::ggplot(sp,ggplot2::aes(x=Sample,y=WithinSamplePct,fill=State))+
 ggplot2::geom_col(width=.8,color='white',linewidth=.15)+
 ggplot2::scale_fill_manual(values=STATE_COLORS)+
 ggplot2::labs(title='Endothelial state composition across 12 libraries',
       subtitle='Biological units: CIT1–HR1 through CIT6–HR6',
       x='Sample',y='Share of sample endothelial cells (%)',fill='State')+
 black_theme(10)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
export(p_comp,'Supplement_EC_state_composition_12samples',S,11,5.8)
p_heat<-ggplot2::ggplot(sp,ggplot2::aes(x=Sample,y=State,fill=ShareOfStatePct))+
 ggplot2::geom_tile(color='white',linewidth=.65)+
 ggplot2::geom_text(ggplot2::aes(label=ifelse(Cells==0,'0',as.character(Cells))),
        color='black',fontface='bold',size=2.6)+
 ggplot2::scale_fill_gradient(low='#F0F4F8',high='#3F86C6')+
 ggplot2::labs(title='Sample contributions to endothelial states',
     subtitle='Labels: cell counts; color: share of each state contributed by sample',
     x='Sample',y=NULL,fill='Share of state (%)')+
 black_theme(9)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
export(p_heat,'Supplement_EC_state_sample_contribution_heatmap',S,11,6.0)

# Paired state fractions: descriptive QC only, FDR over tested states
frac_summary<-lapply(state_order,function(st){
  d<-sp[sp$State==st,c('Donor','Condition','WithinSamplePct')]
  w<-match_donor(d,'WithinSamplePct')
  if (nrow(w) >= 3L && any(w$CIT != w$HR)) {
    p <- suppressWarnings(stats::wilcox.test(w$HR, w$CIT, paired=TRUE, exact=FALSE)$p.value)
  } else if (nrow(w) >= 3L) {
    p <- 1
  } else {
    p <- NA_real_
  }
  data.frame(State=st,Pairs=nrow(w),MeanPctDelta=mean(w$HR-w$CIT),P=p)
})
frac_sum<-do.call(rbind,frac_summary);frac_sum$FDR<-p.adjust(frac_sum$P,'BH')
write.csv(frac_sum,file.path(S,'S4_donor_paired_state_composition_test.csv'),row.names=FALSE)

# Capillary state donor-level marker validation (not independent annotation)
cap_states<-c('Homeostatic capillary EC','Transitional capillary EC')
markers<-intersect(c('CA4','RGCC','EMCN','FCN3','GPIHBP1','EDNRB','HPGD',
                    'SOSTDC1','PLVAP','ICAM1','SERPINE1'),rownames(expr_mat))
cap_md<-md[md$State %in% cap_states,,drop=FALSE]
if(nrow(cap_md) && length(markers)){
 marker_rows<-list();ii<-0L
 for(g in markers){
   v<-as.numeric(expr_mat[g,rownames(cap_md)])
   d<-data.frame(Sample=cap_md$Sample,State=cap_md$State,Value=v)
   tmp<-aggregate(Value~Sample+State,d,mean)
   tmp$Gene<-g
   marker_rows[[length(marker_rows)+1]]<-tmp
 }
 cap_marker<-do.call(rbind,marker_rows)
 write.csv(cap_marker,file.path(S,'S5_capillary_marker_means_by_sample.csv'),row.names=FALSE)
 pmark<-ggplot2::ggplot(cap_marker,
           ggplot2::aes(x=Sample,y=Value,color=State,group=State))+
    ggplot2::geom_point(size=2.3)+
    ggplot2::geom_line(alpha=.45,linewidth=.50)+
    ggplot2::facet_wrap(~Gene,scales='free_y',ncol=4)+
    ggplot2::scale_color_manual(values=STATE_COLORS)+
    ggplot2::labs(title='Capillary EC annotation: marker consistency across donors',
       subtitle='Mean normalized RNA expression per state per sample',
       x=NULL,y='Mean expression',color='Capillary state')+black_theme(9)+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=70,hjust=1,size=7))
 export(pmark,'Supplement_capillary_marker_consistency',S,14,9)
}

# Harmony integrity and sensitivity: saved PCA and Harmony embeddings.
# Compare same hyperparameters; no RunHarmony/FindClusters/FindAllMarkers.
reduc<-SeuratObject::Reductions(endo)
if(USE_EXISTING_PCA_AND_HARMONY && all(c('pca','harmony') %in% reduc)){
  emb1<-SeuratObject::Embeddings(endo,reduction='pca')
  emb2<-SeuratObject::Embeddings(endo,reduction='harmony')
  if(ncol(emb1)>=N_FEATURES_UMAP && ncol(emb2)>=N_FEATURES_UMAP){
    scatter_first2<-function(emb,label,groupby){
      d<-data.frame(X=emb[rownames(md),1],Y=emb[rownames(md),2],
                    Group=md[[groupby]])
      ggplot2::ggplot(d,ggplot2::aes(x=X,y=Y,color=Group))+
        ggplot2::geom_point(size=.16,alpha=.6)+
        ggplot2::labs(title=label,x='Embedding 1',y='Embedding 2',color=groupby)+
        black_theme(9)+ggplot2::theme(legend.position='right')
    }
    b <- scatter_first2(emb1,'PCA before Harmony | sample','Sample')
    a <- scatter_first2(emb2,'Harmony embedding | sample','Sample')
    ab <- b+a+patchwork::plot_layout(ncol=2)
    export(ab,'Supplement_PCA_vs_Harmony_embedding_by_sample',S,13,5.5)
    b2 <- scatter_first2(emb1,'PCA before Harmony | condition','Condition')+
      ggplot2::scale_color_manual(values=COLORS)
    a2 <- scatter_first2(emb2,'Harmony embedding | condition','Condition')+
      ggplot2::scale_color_manual(values=COLORS)
    export(b2+a2+patchwork::plot_layout(ncol=2),
           'Supplement_PCA_vs_Harmony_embedding_by_condition',S,13,5.5)
    if(MAKE_MATCHED_UMAP){
       set.seed(SEED)
       zz<-Seurat::RunUMAP(endo,reduction='pca',dims=1:N_FEATURES_UMAP,
                n.neighbors=20,min.dist=.15,spread=1.30,seed.use=SEED,
                reduction.name='umap_raw_review',reduction.key='RAW_',verbose=FALSE)
       set.seed(SEED)
       zz<-Seurat::RunUMAP(zz,reduction='harmony',dims=1:N_FEATURES_UMAP,
                n.neighbors=20,min.dist=.15,spread=1.30,seed.use=SEED,
                reduction.name='umap_harmony_review',reduction.key='HAR_',verbose=FALSE)
       eu1<-SeuratObject::Embeddings(zz,reduction='umap_raw_review')
       eu2<-SeuratObject::Embeddings(zz,reduction='umap_harmony_review')
       for(gr in c('Sample','Condition','State')){
         p1<-scatter_first2(eu1,paste('Before Harmony UMAP |',gr),gr)
         p2<-scatter_first2(eu2,paste('After Harmony UMAP |',gr),gr)
         if(gr=='Condition'){
           p1<-p1+ggplot2::scale_color_manual(values=COLORS)
           p2<-p2+ggplot2::scale_color_manual(values=COLORS)
         }
         export(p1+p2+patchwork::plot_layout(ncol=2),
           paste0('Supplement_matched_UMAP_pre_post_by_',gr),S,13,5.5)
       }
       rm(zz);invisible(gc())
    }
  } else warning('Stored PCA/Harmony reductions have fewer than 20 dimensions; skip stability embeddings')
} else warning('Stored PCA/Harmony reductions absent; no matched embedding QC performed')


# -------- REVIEWER ROBUSTNESS: PCA-only clustering and retained annotations ---
# A QC sensitivity run, not a replacement for existing Harmony-derived clusters.
ari <- function(x,y){
  tab<-table(x,y);cn2<-function(v) sum(v*(v-1)/2)
  numerator<-cn2(tab);rows<-cn2(rowSums(tab));cols<-cn2(colSums(tab));
  n<-sum(tab);total<-n*(n-1)/2
  expected<-rows*cols/total
  denom<-.5*(rows+cols)-expected
  if(denom==0) return(NA_real_)
  (numerator-expected)/denom
}
if(RUN_PCA_CLUSTER_SENSITIVITY && 'pca' %in% SeuratObject::Reductions(endo)){
  log_message('QC: PCA-only clustering (not replacing original Harmony results)')
  set.seed(SEED)
  endo_unadjusted<-Seurat::FindNeighbors(endo,reduction='pca',
        dims=1:N_FEATURES_UMAP,graph.name=c('review_pre_nn','review_pre_snn'),verbose=FALSE)
  endo_unadjusted<-Seurat::FindClusters(endo_unadjusted,
        graph.name='review_pre_snn',resolution=.6,random.seed=SEED,
        algorithm=1,verbose=FALSE)
  pre_clusters<-as.character(SeuratObject::Idents(endo_unadjusted))
  names(pre_clusters)<-colnames(endo_unadjusted)
  comparison<-table(Original_Annotated_EC_State=md$State,
                    PCA_only_cluster=pre_clusters[rownames(md)])
  write.csv(as.data.frame.matrix(comparison),
      file.path(S,'S6_unintegrated_PCA_only_clusters_vs_original_states.csv'))
  if ('clean_cluster' %in% names(md)) {
    ref_labels <- as.character(md$clean_cluster)
  } else {
    ref_labels <- as.character(md$State)
  }
  ari_value<-ari(ref_labels,pre_clusters[rownames(md)])
  write.csv(data.frame(Comparison='PCA-only clustering vs original Harmony state/cluster',
              AdjustedRandIndex=ari_value,Resolution=.6,Dimensions=N_FEATURES_UMAP),
      file.path(S,'S6_unintegrated_PCA_cluster_agreement_ARI.csv'),row.names=FALSE)
  log_message('PCA-only clustering vs original adjusted Rand index: ',fmt_p(ari_value))
  rm(endo_unadjusted);invisible(gc())
}

# J leave-one-donor-pair-out sensitivity: results are descriptive stability
# of chosen genes, NOT 6 independent external validations.
if(RUN_KEY_LEAVE_ONE_DONOR_OUT && length(complete_pairs)>=4L){
  chosen_genes <- intersect(c('SERPINE1','ICAM1','VCAM1','SELE','PLVAP',
                  'CDH5','CLDN5','JUN','FOS','CXCL8'),key_de$Gene)
  rows <- list()
  for(omit in complete_pairs){
    idx <- key_meta$Donor != omit
    yi <- edgeR::DGEList(counts=as.matrix(key_pb[,idx,drop=FALSE]))
    im <- droplevels(key_meta[idx,,drop=FALSE])
    di <- stats::model.matrix(~Donor+Condition,im)
    ki <- edgeR::filterByExpr(yi,design=di)
    if(sum(ki)<100L) next
    yi <- yi[ki,,keep.lib.sizes=FALSE]
    yi <- edgeR::calcNormFactors(yi)
    yi <- edgeR::estimateDisp(yi,di,robust=TRUE)
    fi <- edgeR::glmQLFit(yi,di,robust=TRUE)
    ti <- edgeR::glmQLFTest(fi,coef='ConditionHR')
    tt <- edgeR::topTags(ti,n=Inf,sort.by='none')$table
    g <- intersect(chosen_genes,rownames(tt))
    rows[[length(rows)+1L]] <- data.frame(OmittedDonor=omit,Gene=g,
        logFC=tt[g,'logFC'],P=tt[g,'PValue'])
  }
  if(length(rows)){
    loo <- do.call(rbind,rows)
    write.csv(loo,file.path(S,'S7_key_state_leave_one_donor_pair_out_DE.csv'),row.names=FALSE)
    longsummary <- do.call(rbind,lapply(split(loo,loo$Gene),function(d) data.frame(
      Gene=d$Gene[1],NumberFits=nrow(d),MeanLogFC=mean(d$logFC),
      FractionDirectionMatchesFull=mean(sign(d$logFC)==sign(key_de$logFC[
       match(d$Gene[1],key_de$Gene)])))))
    write.csv(longsummary,file.path(S,'S7_key_state_leave_one_pair_out_summary.csv'),row.names=FALSE)
    log_message('J leave-one-donor sensitivity generated for ',length(complete_pairs),' donors.')
  }
}

# ----------------------- FINAL A-M ASSEMBLY ---------------------------------
# PRECISE ORIGINAL row/column proportions from re-export source;
# originals A-G and I are references to obj$panels; no computation rerun.
wrap<-function(p) patchwork::wrap_elements(full=p)
p<-lapply(revised,wrap)
row1 <- (p$A+p$B+p$C)+patchwork::plot_layout(widths=c(1,1.34,1.36))
row2 <- (p$D+p$E+p$F)+patchwork::plot_layout(widths=c(1.16,1.2,1))
row3 <- (p$G+p$H)+patchwork::plot_layout(widths=c(.82,2.18))
row4 <- (p$I+p$J+p$K)+patchwork::plot_layout(widths=c(1.04,1.08,1.08))
row5 <- (p$L+p$M)+patchwork::plot_layout(widths=c(1.48,.92))
combo<-(row1/row2/row3/row4/row5)+
        patchwork::plot_layout(heights=c(1.06,1.06,.84,.98,.83))+
        patchwork::plot_annotation(title='Figure 4. Endothelial states and SERPINE1-associated programs',
           theme=ggplot2::theme(plot.title=ggplot2::element_text(
                  face='bold',color='black',hjust=.5,size=16)))
export(combo,'Figure4_REVISED_A_to_M',OUTDIR,18,22.2)
# Save NEW RDS without overwriting the original PREMIUM object
obj$panels<-revised
obj$key_state_DE<-key_de
obj$key_state_Hallmark_GSEA<-gsea_j$table
obj$SERPINE1_high_low_DE_full<-k_de
obj$SERPINE1_high_low_DE_SERPINE1_removed<-k_without
obj$SERPINE1_high_low_Hallmark_GSEA<-gsea_l$table
obj$SERPINE1_high_enrichment_OR<-m_stats
obj$Figure4<-combo
obj$Figure5<-NULL
obj$revision_info<-list(OriginalFile=normalizePath(RDS_FILE,mustWork=TRUE),
   time=Sys.time(),DonorPairs=paste0('CIT',1:6,'-HR',1:6),
   HighExpressionThreshold=cutoff,HighQuantile=SERPINE1_HIGH_QUANTILE,
   KL_StratumMatching='equal high/low cell counts within sample x EC state',
   K_Design='~ Donor + SERPINE1_State (CIT/HR pooled per donor after sample-state balancing)',
   J_Design='~ Donor + Condition',
   M_Design='mean within-donor continuity-corrected logOR; t CI across donor pairs',
   A_to_G_and_I='unaltered original objects',
   FDR='BH within distinct prespecified test families')
saveRDS(obj,file.path(OUTDIR,'Figure4_REVISED_PAIRED_ANALYSIS.rds'),compress=TRUE)
writeLines(capture.output(sessionInfo()),file.path(OUTDIR,'sessionInfo.txt'))
log_message('COMPLETE | New Figure4 A-M PDF+PNG and revised RDS in: ',normalizePath(OUTDIR,winslash='/'))
log_message('IMPORTANT: no A-G/I recomputation, no Seurat/Harmony clustering rerun.')
log_message('Statistical notes: H donor-paired Wilcoxon; J paired edgeR; ',
    'K donor-paired matched-state pseudobulk; L Hallmark; M paired log-odds OR.')



}

# =============================================================================
# STAGE figure5_trajectory_AF
# =============================================================================
PIPELINE$figure5_trajectory_AF <- function() {
# scripts/07_trajectory_panels.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

require_packages(c("monocle","Seurat","SeuratObject","Matrix","dplyr","tidyr","tibble","ggplot2","patchwork","ragg","clusterProfiler","org.Hs.eg.db","AnnotationDbi","ggrepel"))
suppressPackageStartupMessages({library(monocle);library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
restore_named_list(repo_input("trajectory_session"))

# =============================================================================
# Current-session continuation script
#
# IMPORTANT
#   - Does NOT reload packages.
#   - Does NOT reread RDS files.
#   - Does NOT rerun DDRTree/orderCells.
#   - Continues from the objects already present in the current R session.
#
# Required current-session objects
#   cds_fig6
#   trajectory_df
#   rna_data
#
# Recommended existing objects (used if present)
#   senescence_score_obj
#   senescence_genes
#   dynamic_test
#   serpine_module_genes
#   serpine_module_go
#
# Main figure logic
#   A  Capillary endothelial state continuum
#   B  DDRTree pseudotime topology with explicit Early/Late labels
#   C  SERPINE1 distribution on DDRTree
#   D  Reperfusion pathology programs along pseudotime:
#        p53/p21-associated cell-cycle arrest
#        oxidative stress
#        inflammatory/SASP (SERPINE1 excluded)
#        ECM remodeling
#        barrier homeostasis
#   E  Data-driven LATE-UP senescence-related genes across pseudotime
#   F  SERPINE1-associated dynamic module + GO BP
#
# Biological interpretation
#   "premature senescence-like" is treated as a coordinated stress phenotype,
#   not as a claim that the complete canonical senescence score must increase
#   monotonically from early to late pseudotime.
# =============================================================================

options(stringsAsFactors = FALSE)
set.seed(1234)

# =============================================================================
# 0. Input checks
# =============================================================================

required_objects <- c(
  "cds_fig6",
  "trajectory_df",
  "rna_data"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv,
    inherits = FALSE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Current R session is missing required objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

if (utils::packageVersion("monocle") < "2.0.0" || utils::packageVersion("monocle") >= "3.0.0") {
  stop("This continuation script requires Monocle 2.x.")
}

if (utils::packageVersion("igraph") != "1.6.0") {
  warning(
    "Recommended igraph for the current Monocle2 session is 1.6.0. Current: ",
    as.character(utils::packageVersion("igraph"))
  )
}

# =============================================================================
# 1. Output
# =============================================================================

OUTDIR_7 <- repo_output("Figure5/Trajectory")
PANEL_DIR_7 <- file.path(OUTDIR_7, "01_panels")
TABLE_DIR_7 <- file.path(OUTDIR_7, "02_tables")
OBJECT_DIR_7 <- file.path(OUTDIR_7, "03_objects")
FINAL_DIR_7 <- file.path(OUTDIR_7, "04_final")
SUPP_DIR_7 <- file.path(OUTDIR_7, "05_supplementary")

invisible(
  lapply(
    c(
      PANEL_DIR_7,
      TABLE_DIR_7,
      OBJECT_DIR_7,
      FINAL_DIR_7,
      SUPP_DIR_7
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

# =============================================================================
# 2. Parameters
# =============================================================================

CAPILLARY_STATES_7 <- c(
  "Homeostatic capillary EC",
  "Transitional capillary EC",
  "Stress-responsive capillary EC"
)

STRESS_STATE_7 <- "Stress-responsive capillary EC"

CONDITION_COLORS_7 <- c(
  "CIT" = "#2166AC",
  "HR" = "#D73027"
)

STATE_COLORS_7 <- c(
  "Homeostatic capillary EC" = "#43B1D0",
  "Transitional capillary EC" = "#8ED3C7",
  "Stress-responsive capillary EC" = "#E9B000"
)

PSEUDOTIME_COLORS_7 <- c(
  "#253494",
  "#2C7FB8",
  "#7FCDBB",
  "#FEC44F",
  "#E34A33"
)

PROGRAM_COLORS_7 <- c(
  "p53/p21-associated cell-cycle arrest" = "#7B3294",
  "Oxidative stress" = "#008DB8",
  "Inflammatory / SASP" = "#E7298A",
  "ECM remodeling" = "#FF6A00",
  "Barrier homeostasis" = "#009E73"
)

N_PROGRAM_BINS_7 <- 16
N_HEAT_BINS_7 <- 55
N_LATE_HEAT_GENES_7 <- 30
N_MODULE_GENES_7 <- 28
MIN_CELLS_PER_BIN_7 <- 5
MIN_GENE_DETECTION_7 <- 0.01

# =============================================================================
# 3. Helpers
# =============================================================================

theme_fig6_7 <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 1.5,
        hjust = 0.5,
        color = "#202A33"
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 0.5,
        hjust = 0.5,
        color = "#53636F"
      ),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = base_size + 3,
        color = "#202A33"
      ),
      plot.tag.position = c(0.01, 0.99),
      axis.title = ggplot2::element_text(
        face = "bold",
        color = "#202A33"
      ),
      axis.text = ggplot2::element_text(
        color = "#202A33"
      ),
      panel.grid.major = ggplot2::element_line(
        color = "#E8EDF0",
        linewidth = 0.28
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill = NA,
        color = "#CBD4DA",
        linewidth = 0.45
      ),
      strip.background = ggplot2::element_rect(
        fill = "#F1F5F7",
        color = "#CBD4DA",
        linewidth = 0.42
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        color = "#202A33"
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      )
    )
}


theme_traj_7 <- function(base_size = 8) {
  theme_fig6_7(base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}


format_p_7 <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}


zscore_safe_7 <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (sum(is.finite(x)) < 2 || !is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}


build_score_7 <- function(
    expression_matrix,
    genes,
    cells,
    min_detection = 0.01,
    minimum_genes = 3
) {
  genes <- unique(intersect(genes, rownames(expression_matrix)))
  cells <- unique(intersect(cells, colnames(expression_matrix)))

  if (length(genes) < minimum_genes) {
    return(
      list(
        score = setNames(rep(NA_real_, length(cells)), cells),
        genes = genes,
        z_matrix = matrix(
          numeric(0),
          nrow = 0,
          ncol = length(cells),
          dimnames = list(NULL, cells)
        )
      )
    )
  }

  mat <- as.matrix(
    expression_matrix[
      genes,
      cells,
      drop = FALSE
    ]
  )

  detection <- rowMeans(mat > 0, na.rm = TRUE)
  mat <- mat[
    is.finite(detection) & detection >= min_detection,
    ,
    drop = FALSE
  ]

  if (nrow(mat) < minimum_genes) {
    return(
      list(
        score = setNames(rep(NA_real_, length(cells)), cells),
        genes = rownames(mat),
        z_matrix = matrix(
          numeric(0),
          nrow = 0,
          ncol = length(cells),
          dimnames = list(NULL, cells)
        )
      )
    )
  }

  zmat <- matrix(
    0,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )

  for (i in seq_len(nrow(mat))) {
    zmat[i, ] <- zscore_safe_7(mat[i, ])
  }

  score <- colMeans(zmat, na.rm = TRUE)
  names(score) <- colnames(zmat)

  list(
    score = score,
    genes = rownames(mat),
    z_matrix = zmat
  )
}


save_panel_7 <- function(plot_object, filename, width, height, dpi = 600) {
  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggplot2::ggsave(
    filename = file.path(PANEL_DIR_7, paste0(filename, ".pdf")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = pdf_device,
    bg = "white",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    filename = file.path(PANEL_DIR_7, paste0(filename, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = ragg::agg_png,
    bg = "white",
    limitsize = FALSE
  )
}

# =============================================================================
# 4. Standardize current trajectory metadata
# =============================================================================

traj7 <- trajectory_df

if (!"Cell" %in% colnames(traj7)) {
  traj7$Cell <- rownames(traj7)
}

if (!"PseudotimeNorm" %in% colnames(traj7)) {
  if (!"Pseudotime" %in% colnames(traj7)) {
    stop("trajectory_df has neither PseudotimeNorm nor Pseudotime.")
  }
  pt_min_7 <- min(traj7$Pseudotime, na.rm = TRUE)
  pt_max_7 <- max(traj7$Pseudotime, na.rm = TRUE)
  traj7$PseudotimeNorm <- (
    traj7$Pseudotime - pt_min_7
  ) / (
    pt_max_7 - pt_min_7
  )
}

if (!"Sample" %in% colnames(traj7)) {
  sample_candidates_7 <- c("Fig6_Sample", "orig.ident")
  sample_hit_7 <- sample_candidates_7[sample_candidates_7 %in% colnames(traj7)][1]
  if (is.na(sample_hit_7)) stop("Cannot find sample metadata in trajectory_df.")
  traj7$Sample <- as.character(traj7[[sample_hit_7]])
}

if (!"Condition" %in% colnames(traj7)) {
  condition_candidates_7 <- c("Fig6_Condition", "tissue_type")
  condition_hit_7 <- condition_candidates_7[condition_candidates_7 %in% colnames(traj7)][1]
  if (is.na(condition_hit_7)) stop("Cannot find condition metadata in trajectory_df.")
  traj7$Condition <- as.character(traj7[[condition_hit_7]])
}

traj7$Condition <- factor(
  as.character(traj7$Condition),
  levels = c("CIT", "HR")
)

if (!"EndothelialState" %in% colnames(traj7)) {
  state_candidates_7 <- c(
    "Fig6_State",
    "Endothelial_Annotation",
    "Annotation"
  )
  state_hit_7 <- state_candidates_7[state_candidates_7 %in% colnames(traj7)][1]
  if (is.na(state_hit_7)) stop("Cannot find endothelial-state metadata in trajectory_df.")
  traj7$EndothelialState <- as.character(traj7[[state_hit_7]])
}

traj7$EndothelialState <- factor(
  as.character(traj7$EndothelialState),
  levels = CAPILLARY_STATES_7
)

coord1_7 <- c("DDRTree_1", "DDRTree1")[c("DDRTree_1", "DDRTree1") %in% colnames(traj7)][1]
coord2_7 <- c("DDRTree_2", "DDRTree2")[c("DDRTree_2", "DDRTree2") %in% colnames(traj7)][1]

if (is.na(coord1_7) || is.na(coord2_7)) {
  coords7 <- as.data.frame(
    t(
      monocle::reducedDimS(cds_fig6)
    )
  )
  coords7$Cell <- rownames(coords7)
  colnames(coords7)[1:2] <- c("DDRTree_1", "DDRTree_2")

  traj7 <- traj7 %>%
    dplyr::left_join(
      coords7[, c("Cell", "DDRTree_1", "DDRTree_2")],
      by = "Cell"
    )

  coord1_7 <- "DDRTree_1"
  coord2_7 <- "DDRTree_2"
}

if (!"SERPINE1" %in% colnames(traj7)) {
  if (!"SERPINE1" %in% rownames(rna_data)) stop("SERPINE1 is not present in rna_data.")
  traj7$SERPINE1 <- as.numeric(
    rna_data[
      "SERPINE1",
      traj7$Cell,
      drop = TRUE
    ]
  )
}

# =============================================================================
# 5. Gene programs focused on premature-senescence-like pathology
# =============================================================================

# Independent p53/p21-cell-cycle-arrest machinery.
# SERPINE1 is intentionally NOT included here.
P53_P21_ARREST_GENES_7 <- c(
  "CDKN1A",
  "GADD45A",
  "GADD45B",
  "BTG2",
  "DDB2",
  "MDM2",
  "SESN1",
  "SESN2",
  "RRM2B",
  "BBC3",
  "TP53I3",
  "GDF15"
)

OXIDATIVE_STRESS_GENES_7 <- c(
  "HMOX1",
  "NQO1",
  "SOD1",
  "SOD2",
  "GPX1",
  "GPX3",
  "TXNRD1",
  "PRDX1",
  "PRDX2",
  "GCLC",
  "GCLM",
  "ATF3",
  "DUSP1",
  "JUN",
  "FOS"
)

# Independent inflammatory/SASP score.
# SERPINE1 is excluded to avoid self-inclusion.
INFLAMMATORY_SASP_GENES_7 <- c(
  "ICAM1",
  "VCAM1",
  "SELE",
  "SELP",
  "IL6",
  "CXCL8",
  "CXCL1",
  "CXCL2",
  "CXCL3",
  "CCL2",
  "TNFAIP3",
  "NFKBIA",
  "PTGS2",
  "MMP9"
)

ECM_REMODELING_GENES_7 <- c(
  "FN1",
  "COL4A1",
  "COL4A2",
  "MMP14",
  "MMP2",
  "PLAU",
  "PLAUR",
  "ANGPT2",
  "ICAM1",
  "VCAM1"
)

BARRIER_HOMEOSTASIS_GENES_7 <- c(
  "CLDN5",
  "CDH5",
  "TJP1",
  "OCLN",
  "TEK",
  "KDR",
  "ROBO4",
  "ESAM"
)

trajectory_cells_7 <- intersect(
  traj7$Cell,
  colnames(rna_data)
)

p53_obj_7 <- build_score_7(
  expression_matrix = rna_data,
  genes = P53_P21_ARREST_GENES_7,
  cells = trajectory_cells_7,
  min_detection = MIN_GENE_DETECTION_7,
  minimum_genes = 3
)

oxidative_obj_7 <- build_score_7(
  expression_matrix = rna_data,
  genes = OXIDATIVE_STRESS_GENES_7,
  cells = trajectory_cells_7,
  min_detection = MIN_GENE_DETECTION_7,
  minimum_genes = 4
)

sasp_obj_7 <- build_score_7(
  expression_matrix = rna_data,
  genes = INFLAMMATORY_SASP_GENES_7,
  cells = trajectory_cells_7,
  min_detection = MIN_GENE_DETECTION_7,
  minimum_genes = 3
)

ecm_obj_7 <- build_score_7(
  expression_matrix = rna_data,
  genes = ECM_REMODELING_GENES_7,
  cells = trajectory_cells_7,
  min_detection = MIN_GENE_DETECTION_7,
  minimum_genes = 3
)

barrier_obj_7 <- build_score_7(
  expression_matrix = rna_data,
  genes = BARRIER_HOMEOSTASIS_GENES_7,
  cells = trajectory_cells_7,
  min_detection = MIN_GENE_DETECTION_7,
  minimum_genes = 3
)

traj7$p53_p21_arrest <- as.numeric(
  p53_obj_7$score[
    traj7$Cell
  ]
)

traj7$Oxidative_stress <- as.numeric(
  oxidative_obj_7$score[
    traj7$Cell
  ]
)

traj7$Inflammatory_SASP <- as.numeric(
  sasp_obj_7$score[
    traj7$Cell
  ]
)

traj7$ECM_remodeling <- as.numeric(
  ecm_obj_7$score[
    traj7$Cell
  ]
)

traj7$Barrier_homeostasis <- as.numeric(
  barrier_obj_7$score[
    traj7$Cell
  ]
)

# Higher value means stronger remodeling with relative loss of barrier program.
traj7$ECM_barrier_injury <- (
  zscore_safe_7(traj7$ECM_remodeling) -
    zscore_safe_7(traj7$Barrier_homeostasis)
) / 2

write.csv(
  data.frame(
    Program = c(
      rep("p53/p21-associated cell-cycle arrest", length(p53_obj_7$genes)),
      rep("Oxidative stress", length(oxidative_obj_7$genes)),
      rep("Inflammatory / SASP", length(sasp_obj_7$genes)),
      rep("ECM remodeling", length(ecm_obj_7$genes)),
      rep("Barrier homeostasis", length(barrier_obj_7$genes))
    ),
    Gene = c(
      p53_obj_7$genes,
      oxidative_obj_7$genes,
      sasp_obj_7$genes,
      ecm_obj_7$genes,
      barrier_obj_7$genes
    ),
    stringsAsFactors = FALSE
  ),
  file.path(TABLE_DIR_7, "Figure5_program_genes_used.csv"),
  row.names = FALSE
)

# =============================================================================
# 6. Panel A | state continuum
# =============================================================================

cds_state_7 <- cds_fig6
Biobase::pData(cds_state_7)$StateLabel7 <- as.character(
  traj7$EndothelialState[
    match(
      rownames(Biobase::pData(cds_state_7)),
      traj7$Cell
    )
  ]
)

state_centers_7 <- traj7 %>%
  dplyr::filter(!is.na(EndothelialState)) %>%
  dplyr::group_by(EndothelialState) %>%
  dplyr::summarise(
    X = stats::median(.data[[coord1_7]], na.rm = TRUE),
    Y = stats::median(.data[[coord2_7]], na.rm = TRUE),
    .groups = "drop"
  )

p7A <- monocle::plot_cell_trajectory(
  cds_state_7,
  color_by = "StateLabel7",
  show_tree = TRUE,
  show_backbone = TRUE,
  backbone_color = "#263640",
  show_branch_points = TRUE,
  show_state_number = FALSE,
  cell_size = 0.58,
  cell_link_size = 0.62
) +
  ggplot2::scale_color_manual(
    values = STATE_COLORS_7,
    drop = FALSE
  ) +
  ggrepel::geom_label_repel(
    data = state_centers_7,
    ggplot2::aes(
      x = X,
      y = Y,
      label = EndothelialState
    ),
    inherit.aes = FALSE,
    size = 2.25,
    fontface = "bold",
    fill = scales::alpha("white", 0.90),
    label.size = 0.18,
    box.padding = 0.28,
    point.padding = 0.10,
    min.segment.length = 0,
    segment.color = "#7B8790",
    seed = 1234,
    show.legend = FALSE
  ) +
  ggplot2::labs(
    tag = "A",
    title = "Capillary endothelial injury continuum",
    subtitle = "Homeostatic-to-stress-responsive endothelial state transition",
    x = "DDRTree 1",
    y = "DDRTree 2",
    color = NULL
  ) +
  theme_traj_7(8) +
  ggplot2::theme(
    legend.position = "right",
    legend.text = ggplot2::element_text(size = 6.2)
  )

# =============================================================================
# 7. Panel B | pseudotime topology with Early/Late labels
# =============================================================================

low_pt_cells_7 <- traj7 %>%
  dplyr::filter(
    PseudotimeNorm <= stats::quantile(PseudotimeNorm, 0.02, na.rm = TRUE)
  )

high_pt_cells_7 <- traj7 %>%
  dplyr::filter(
    PseudotimeNorm >= stats::quantile(PseudotimeNorm, 0.98, na.rm = TRUE)
  )

pt_labels_7 <- data.frame(
  Label = c("Early", "Late"),
  X = c(
    stats::median(low_pt_cells_7[[coord1_7]], na.rm = TRUE),
    stats::median(high_pt_cells_7[[coord1_7]], na.rm = TRUE)
  ),
  Y = c(
    stats::median(low_pt_cells_7[[coord2_7]], na.rm = TRUE),
    stats::median(high_pt_cells_7[[coord2_7]], na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

p7B <- monocle::plot_cell_trajectory(
  cds_fig6,
  color_by = "Pseudotime",
  show_tree = TRUE,
  show_backbone = TRUE,
  backbone_color = "#263640",
  show_branch_points = TRUE,
  show_state_number = TRUE,
  cell_size = 0.60,
  cell_link_size = 0.62
) +
  ggplot2::scale_color_gradientn(
    colors = PSEUDOTIME_COLORS_7,
    name = "Pseudotime"
  ) +
  ggrepel::geom_label_repel(
    data = pt_labels_7,
    ggplot2::aes(
      x = X,
      y = Y,
      label = Label
    ),
    inherit.aes = FALSE,
    size = 2.7,
    fontface = "bold",
    fill = scales::alpha("white", 0.92),
    label.size = 0.20,
    segment.color = "#616B73",
    seed = 1234,
    show.legend = FALSE
  ) +
  ggplot2::labs(
    tag = "B",
    title = "Pseudotime topology on DDRTree",
    subtitle = "Early-to-late transcriptional-state ordering; pseudotime is not physical clock time",
    x = "DDRTree 1",
    y = "DDRTree 2"
  ) +
  theme_traj_7(8)

# =============================================================================
# 8. Panel C | SERPINE1 distribution on DDRTree
# =============================================================================

cds_serpine_7 <- cds_fig6
Biobase::pData(cds_serpine_7)$SERPINE1_plot <- traj7$SERPINE1[
  match(
    rownames(Biobase::pData(cds_serpine_7)),
    traj7$Cell
  )
]

p7C <- monocle::plot_cell_trajectory(
  cds_serpine_7,
  color_by = "SERPINE1_plot",
  show_tree = TRUE,
  show_backbone = TRUE,
  backbone_color = "#37474F",
  show_branch_points = TRUE,
  show_state_number = FALSE,
  cell_size = 0.62,
  cell_link_size = 0.62
) +
  ggplot2::scale_color_gradientn(
    colors = c(
      "#DCE6EF",
      "#F7F7F7",
      "#F4A582",
      "#D6604D",
      "#B2182B"
    ),
    name = "SERPINE1\nexpression"
  ) +
  ggplot2::labs(
    tag = "C",
    title = expression(
      italic(SERPINE1) * " distribution on the injury trajectory"
    ),
    subtitle = "High expression localizes toward the late stress-responsive trajectory arm",
    x = "DDRTree 1",
    y = "DDRTree 2"
  ) +
  theme_traj_7(8)

# =============================================================================
# 9. Panel D | reperfusion pathology programs along pseudotime
# =============================================================================

traj7$ProgramBin <- cut(
  traj7$PseudotimeNorm,
  breaks = seq(
    0,
    1,
    length.out = N_PROGRAM_BINS_7 + 1
  ),
  include.lowest = TRUE,
  labels = FALSE
)

program_sample_bin_7 <- traj7 %>%
  dplyr::group_by(
    Sample,
    Condition,
    ProgramBin
  ) %>%
  dplyr::summarise(
    Cells = dplyr::n(),
    Pseudotime = mean(PseudotimeNorm, na.rm = TRUE),
    `p53/p21-associated cell-cycle arrest` = mean(p53_p21_arrest, na.rm = TRUE),
    `Oxidative stress` = mean(Oxidative_stress, na.rm = TRUE),
    `Inflammatory / SASP` = mean(Inflammatory_SASP, na.rm = TRUE),
    `ECM remodeling` = mean(ECM_remodeling, na.rm = TRUE),
    `Barrier homeostasis` = mean(Barrier_homeostasis, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    Cells >= MIN_CELLS_PER_BIN_7
  )

program_long_7 <- program_sample_bin_7 %>%
  tidyr::pivot_longer(
    cols = c(
      `p53/p21-associated cell-cycle arrest`,
      `Oxidative stress`,
      `Inflammatory / SASP`,
      `ECM remodeling`,
      `Barrier homeostasis`
    ),
    names_to = "Program",
    values_to = "Score"
  )

program_long_7$Program <- factor(
  program_long_7$Program,
  levels = names(PROGRAM_COLORS_7)
)

p7D <- ggplot2::ggplot(
  program_long_7,
  ggplot2::aes(
    x = Pseudotime,
    y = Score,
    color = Condition,
    fill = Condition
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#98A3AB",
    linewidth = 0.35
  ) +
  ggplot2::geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,
    linewidth = 0.95,
    alpha = 0.10
  ) +
  ggplot2::facet_wrap(
    ~ Program,
    nrow = 1,
    scales = "free_y"
  ) +
  ggplot2::scale_color_manual(
    values = CONDITION_COLORS_7
  ) +
  ggplot2::scale_fill_manual(
    values = CONDITION_COLORS_7
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  ggplot2::labs(
    tag = "D",
    title = "Reperfusion-associated premature-senescence programs along pseudotime",
    subtitle = "Independent arrest, oxidative-stress, SASP, ECM-remodeling and barrier-homeostasis programs",
    x = "Normalized pseudotime",
    y = "Program score",
    color = NULL,
    fill = NULL
  ) +
  theme_fig6_7(7.2) +
  ggplot2::theme(
    legend.position = "top",
    strip.text = ggplot2::element_text(size = 5.6)
  )

# =============================================================================
# 10. Canonical cellular-senescence genes
#     Used for data-driven late-up selection in Panel E.
# =============================================================================

if (
  exists("senescence_score_obj", envir = .GlobalEnv, inherits = FALSE) &&
    !is.null(senescence_score_obj$genes)
) {
  senescence_genes_7 <- unique(
    as.character(senescence_score_obj$genes)
  )
} else if (
  exists("senescence_genes", envir = .GlobalEnv, inherits = FALSE)
) {
  senescence_genes_7 <- unique(
    as.character(senescence_genes)
  )
} else {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(
      "No existing senescence gene set found and msigdbr is unavailable."
    )
  }

  msig_c2_7 <- msigdbr::msigdbr(
    species = "Homo sapiens",
    category = "C2"
  )

  sen_name_hits_7 <- unique(
    msig_c2_7$gs_name[
      grepl(
        "CELLULAR_SENESCENCE",
        msig_c2_7$gs_name,
        ignore.case = TRUE
      )
    ]
  )

  if (length(sen_name_hits_7) < 1) {
    stop("Cannot locate a cellular-senescence gene set in msigdbr C2.")
  }

  preferred_sen_name_7 <- if (
    "REACTOME_CELLULAR_SENESCENCE" %in% sen_name_hits_7
  ) {
    "REACTOME_CELLULAR_SENESCENCE"
  } else {
    sen_name_hits_7[1]
  }

  senescence_genes_7 <- unique(
    as.character(
      msig_c2_7$gene_symbol[
        msig_c2_7$gs_name == preferred_sen_name_7
      ]
    )
  )
}

senescence_genes_7 <- unique(
  intersect(
    senescence_genes_7,
    rownames(rna_data)
  )
)

if (!"SERPINE1" %in% senescence_genes_7 && "SERPINE1" %in% rownames(rna_data)) {
  senescence_genes_7 <- unique(
    c(
      "SERPINE1",
      senescence_genes_7
    )
  )
}

# =============================================================================
# 11. Panel E | LATE-UP senescence-related genes
#
# Selection is data-driven:
#   1. test canonical senescence genes versus pseudotime
#   2. compute late-minus-early expression delta
#   3. prioritize q<0.05 AND late delta > 0
#   4. if too few, use the largest positive late deltas
# =============================================================================

fd7 <- Biobase::fData(cds_fig6)
if (!"gene_short_name" %in% colnames(fd7)) {
  fd7$gene_short_name <- rownames(fd7)
}

sen_gene_ids_7 <- unique(
  rownames(fd7)[
    rownames(fd7) %in% senescence_genes_7 |
      as.character(fd7$gene_short_name) %in% senescence_genes_7
  ]
)

sen_gene_ids_7 <- intersect(
  sen_gene_ids_7,
  rownames(cds_fig6)
)

sen_dynamic_test_7 <- monocle::differentialGeneTest(
  cds_fig6[
    sen_gene_ids_7,
  ],
  fullModelFormulaStr = "~sm.ns(Pseudotime)",
  cores = 1
)

sen_dynamic_test_7 <- as.data.frame(
  sen_dynamic_test_7
)

sen_dynamic_test_7$GeneID <- rownames(
  sen_dynamic_test_7
)

sen_dynamic_test_7$Gene <- as.character(
  sen_dynamic_test_7$gene_short_name
)

sen_dynamic_test_7$Gene[
  is.na(sen_dynamic_test_7$Gene) |
    !nzchar(sen_dynamic_test_7$Gene)
] <- sen_dynamic_test_7$GeneID[
  is.na(sen_dynamic_test_7$Gene) |
    !nzchar(sen_dynamic_test_7$Gene)
]

sen_dynamic_test_7 <- sen_dynamic_test_7 %>%
  dplyr::distinct(
    Gene,
    .keep_all = TRUE
  )

sen_expr_genes_7 <- intersect(
  sen_dynamic_test_7$Gene,
  rownames(rna_data)
)

sen_expr_mat_7 <- as.matrix(
  rna_data[
    sen_expr_genes_7,
    traj7$Cell,
    drop = FALSE
  ]
)

sen_detect_7 <- rowMeans(
  sen_expr_mat_7 > 0,
  na.rm = TRUE
)

sen_expr_mat_7 <- sen_expr_mat_7[
  sen_detect_7 >= MIN_GENE_DETECTION_7,
  ,
  drop = FALSE
]

sen_expr_genes_7 <- rownames(
  sen_expr_mat_7
)

early_cells_7 <- traj7$Cell[
  traj7$PseudotimeNorm <= 0.20
]

late_cells_7 <- traj7$Cell[
  traj7$PseudotimeNorm >= 0.80
]

late_delta_7 <- rowMeans(
  sen_expr_mat_7[
    ,
    late_cells_7,
    drop = FALSE
  ],
  na.rm = TRUE
) - rowMeans(
  sen_expr_mat_7[
    ,
    early_cells_7,
    drop = FALSE
  ],
  na.rm = TRUE
)

late_gene_table_7 <- data.frame(
  Gene = names(late_delta_7),
  LateMinusEarly = as.numeric(late_delta_7),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    sen_dynamic_test_7 %>%
      dplyr::select(
        Gene,
        pval,
        qval
      ),
    by = "Gene"
  ) %>%
  dplyr::arrange(
    dplyr::desc(LateMinusEarly),
    qval,
    pval
  )

late_sig_genes_7 <- late_gene_table_7 %>%
  dplyr::filter(
    is.finite(qval),
    qval < 0.05,
    LateMinusEarly > 0
  ) %>%
  dplyr::pull(Gene)

if (length(late_sig_genes_7) >= 12) {
  late_heat_genes_7 <- head(
    late_sig_genes_7,
    N_LATE_HEAT_GENES_7
  )
  late_selection_rule_7 <- "q<0.05 and late > early"
} else {
  late_heat_genes_7 <- late_gene_table_7 %>%
    dplyr::filter(
      LateMinusEarly > 0
    ) %>%
    dplyr::slice_head(
      n = N_LATE_HEAT_GENES_7
    ) %>%
    dplyr::pull(Gene)
  late_selection_rule_7 <- "top positive late-minus-early senescence genes"
}

if (
  "SERPINE1" %in% late_gene_table_7$Gene &&
    late_gene_table_7$LateMinusEarly[
      late_gene_table_7$Gene == "SERPINE1"
    ][1] > 0
) {
  late_heat_genes_7 <- unique(
    c(
      "SERPINE1",
      late_heat_genes_7
    )
  )
}

late_heat_genes_7 <- head(
  late_heat_genes_7,
  N_LATE_HEAT_GENES_7
)

write.csv(
  late_gene_table_7,
  file.path(TABLE_DIR_7, "PanelE_late_up_senescence_gene_statistics.csv"),
  row.names = FALSE
)

traj7$HeatBin7 <- cut(
  traj7$PseudotimeNorm,
  breaks = seq(
    0,
    1,
    length.out = N_HEAT_BINS_7 + 1
  ),
  include.lowest = TRUE,
  labels = FALSE
)

heat_bins_7 <- sort(
  unique(
    traj7$HeatBin7
  )
)

late_heat_matrix_7 <- sapply(
  heat_bins_7,
  function(bin_i) {
    cells_i <- traj7$Cell[
      traj7$HeatBin7 == bin_i
    ]

    Matrix::rowMeans(
      rna_data[
        late_heat_genes_7,
        cells_i,
        drop = FALSE
      ]
    )
  }
)

rownames(late_heat_matrix_7) <- late_heat_genes_7

late_heat_z_7 <- matrix(
  0,
  nrow = nrow(late_heat_matrix_7),
  ncol = ncol(late_heat_matrix_7),
  dimnames = dimnames(late_heat_matrix_7)
)

for (i in seq_len(nrow(late_heat_matrix_7))) {
  late_heat_z_7[i, ] <- zscore_safe_7(
    late_heat_matrix_7[i, ]
  )
}

late_heat_order_7 <- rownames(late_heat_z_7)[
  order(
    late_gene_table_7$LateMinusEarly[
      match(
        rownames(late_heat_z_7),
        late_gene_table_7$Gene
      )
    ]
  )
]

late_heat_z_7 <- late_heat_z_7[
  late_heat_order_7,
  ,
  drop = FALSE
]

late_heat_long_7 <- as.data.frame(
  as.table(
    late_heat_z_7
  ),
  stringsAsFactors = FALSE
)

colnames(late_heat_long_7) <- c(
  "Gene",
  "Bin",
  "Zscore"
)

late_heat_long_7$BinIndex <- match(
  late_heat_long_7$Bin,
  colnames(late_heat_z_7)
)

if (all(is.na(late_heat_long_7$BinIndex))) {
  late_heat_long_7$BinIndex <- match(
    late_heat_long_7$Bin,
    unique(late_heat_long_7$Bin)
  )
}

late_heat_long_7$Pseudotime <- (
  late_heat_long_7$BinIndex - 0.5
) / ncol(late_heat_z_7)

late_heat_long_7$Gene <- factor(
  late_heat_long_7$Gene,
  levels = rev(
    rownames(late_heat_z_7)
  )
)

serpine_y_7 <- which(
  levels(late_heat_long_7$Gene) == "SERPINE1"
)

p7E <- ggplot2::ggplot(
  late_heat_long_7,
  ggplot2::aes(
    x = Pseudotime,
    y = Gene,
    fill = Zscore
  )
) +
  ggplot2::geom_tile(linewidth = 0) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Expression\n(z-score)"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1.10),
    breaks = c(0, 0.5, 1),
    labels = c("Early", "Mid", "Late"),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    tag = "E",
    title = "Late-activated senescence-related genes",
    subtitle = paste0(
      "Data-driven selection: ",
      late_selection_rule_7
    ),
    x = "Endothelial injury pseudotime",
    y = NULL
  ) +
  theme_fig6_7(6.7) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      face = "italic",
      size = 5.2
    ),
    axis.ticks.y = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    legend.position = "bottom"
  )

if (length(serpine_y_7) == 1) {
  p7E <- p7E +
    ggplot2::annotate(
      "segment",
      x = 1.005,
      xend = 1.045,
      y = serpine_y_7,
      yend = serpine_y_7,
      color = "#B2182B",
      linewidth = 0.8,
      arrow = grid::arrow(
        length = grid::unit(1.4, "mm")
      )
    ) +
    ggplot2::annotate(
      "text",
      x = 1.052,
      y = serpine_y_7,
      label = "SERPINE1",
      hjust = 0,
      color = "#B2182B",
      fontface = "bold.italic",
      size = 2.35
    )
}

# =============================================================================
# 12. Panel F | SERPINE1-associated dynamic module + GO
# =============================================================================

# Use the already-derived module if available.
if (
  exists("serpine_module_genes", envir = .GlobalEnv, inherits = FALSE) &&
    nrow(as.data.frame(serpine_module_genes)) > 0
) {
  module_gene_table_7 <- as.data.frame(
    serpine_module_genes
  )

  if ("Gene" %in% colnames(module_gene_table_7)) {
    module_genes_7 <- unique(
      as.character(module_gene_table_7$Gene)
    )
  } else {
    module_genes_7 <- rownames(
      module_gene_table_7
    )
  }
} else {
  if (!exists("dynamic_test", envir = .GlobalEnv, inherits = FALSE)) {
    stop(
      "Panel F requires serpine_module_genes or dynamic_test from the current session."
    )
  }

  dynamic7 <- as.data.frame(dynamic_test)

  if (!"Gene" %in% colnames(dynamic7)) {
    if ("gene_short_name" %in% colnames(dynamic7)) {
      dynamic7$Gene <- as.character(dynamic7$gene_short_name)
    } else {
      dynamic7$Gene <- rownames(dynamic7)
    }
  }

  dynamic_genes_7 <- dynamic7 %>%
    dplyr::filter(
      is.finite(qval),
      qval < 0.05,
      Gene %in% rownames(rna_data)
    ) %>%
    dplyr::arrange(qval, pval) %>%
    dplyr::slice_head(n = 250) %>%
    dplyr::pull(Gene)

  dynamic_genes_7 <- unique(
    c(
      "SERPINE1",
      dynamic_genes_7
    )
  )

  dynamic_expr_7 <- as.matrix(
    rna_data[
      dynamic_genes_7,
      traj7$Cell,
      drop = FALSE
    ]
  )

  dynamic_z_7 <- t(
    apply(
      dynamic_expr_7,
      1,
      zscore_safe_7
    )
  )

  rownames(dynamic_z_7) <- rownames(dynamic_expr_7)
  colnames(dynamic_z_7) <- colnames(dynamic_expr_7)

  serpine_profile_7 <- dynamic_z_7[
    "SERPINE1",
  ]

  serp_cor_7 <- apply(
    dynamic_z_7,
    1,
    function(x) {
      suppressWarnings(
        stats::cor(
          x,
          serpine_profile_7,
          method = "spearman",
          use = "pairwise.complete.obs"
        )
      )
    }
  )

  module_genes_7 <- names(
    sort(
      serp_cor_7,
      decreasing = TRUE
    )
  )

  module_genes_7 <- module_genes_7[
    is.finite(serp_cor_7[module_genes_7]) &
      serp_cor_7[module_genes_7] > 0.25
  ]

  module_genes_7 <- unique(
    c(
      "SERPINE1",
      module_genes_7
    )
  )
}

module_genes_7 <- unique(
  intersect(
    module_genes_7,
    rownames(rna_data)
  )
)

# Rank module genes by late-minus-early difference for display only.
module_expr_7 <- as.matrix(
  rna_data[
    module_genes_7,
    traj7$Cell,
    drop = FALSE
  ]
)

module_late_delta_7 <- rowMeans(
  module_expr_7[, late_cells_7, drop = FALSE],
  na.rm = TRUE
) - rowMeans(
  module_expr_7[, early_cells_7, drop = FALSE],
  na.rm = TRUE
)

module_display_genes_7 <- names(
  sort(
    module_late_delta_7,
    decreasing = TRUE
  )
)

module_display_genes_7 <- head(
  module_display_genes_7,
  N_MODULE_GENES_7
)

if (
  "SERPINE1" %in% module_genes_7 &&
    !"SERPINE1" %in% module_display_genes_7
) {
  module_display_genes_7 <- unique(
    c(
      "SERPINE1",
      module_display_genes_7
    )
  )
  module_display_genes_7 <- head(
    module_display_genes_7,
    N_MODULE_GENES_7
  )
}

module_heat_matrix_7 <- sapply(
  heat_bins_7,
  function(bin_i) {
    cells_i <- traj7$Cell[
      traj7$HeatBin7 == bin_i
    ]

    Matrix::rowMeans(
      rna_data[
        module_display_genes_7,
        cells_i,
        drop = FALSE
      ]
    )
  }
)

rownames(module_heat_matrix_7) <- module_display_genes_7

module_heat_z_7 <- matrix(
  0,
  nrow = nrow(module_heat_matrix_7),
  ncol = ncol(module_heat_matrix_7),
  dimnames = dimnames(module_heat_matrix_7)
)

for (i in seq_len(nrow(module_heat_matrix_7))) {
  module_heat_z_7[i, ] <- zscore_safe_7(
    module_heat_matrix_7[i, ]
  )
}

module_heat_order_7 <- rownames(module_heat_z_7)[
  order(
    module_late_delta_7[
      rownames(module_heat_z_7)
    ]
  )
]

module_heat_z_7 <- module_heat_z_7[
  module_heat_order_7,
  ,
  drop = FALSE
]

module_heat_long_7 <- as.data.frame(
  as.table(module_heat_z_7),
  stringsAsFactors = FALSE
)

colnames(module_heat_long_7) <- c(
  "Gene",
  "Bin",
  "Zscore"
)

module_heat_long_7$BinIndex <- match(
  module_heat_long_7$Bin,
  colnames(module_heat_z_7)
)

if (all(is.na(module_heat_long_7$BinIndex))) {
  module_heat_long_7$BinIndex <- match(
    module_heat_long_7$Bin,
    unique(module_heat_long_7$Bin)
  )
}

module_heat_long_7$Pseudotime <- (
  module_heat_long_7$BinIndex - 0.5
) / ncol(module_heat_z_7)

module_heat_long_7$Gene <- factor(
  module_heat_long_7$Gene,
  levels = rev(
    rownames(module_heat_z_7)
  )
)

serpine_module_y_7 <- which(
  levels(module_heat_long_7$Gene) == "SERPINE1"
)

# Try to retrieve SERPINE1 pseudotime P/Q from existing dynamic_test.
serpine_dynamic_p_7 <- NA_real_
serpine_dynamic_q_7 <- NA_real_

if (exists("dynamic_test", envir = .GlobalEnv, inherits = FALSE)) {
  dynamic7 <- as.data.frame(dynamic_test)

  if (!"Gene" %in% colnames(dynamic7)) {
    if ("gene_short_name" %in% colnames(dynamic7)) {
      dynamic7$Gene <- as.character(dynamic7$gene_short_name)
    } else {
      dynamic7$Gene <- rownames(dynamic7)
    }
  }

  serp_row_7 <- dynamic7[
    dynamic7$Gene == "SERPINE1",
    ,
    drop = FALSE
  ]

  if (nrow(serp_row_7) >= 1) {
    serpine_dynamic_p_7 <- serp_row_7$pval[1]
    serpine_dynamic_q_7 <- serp_row_7$qval[1]
  }
}

if (!is.finite(serpine_dynamic_p_7) || !is.finite(serpine_dynamic_q_7)) {
  serpine_gene_id_7 <- rownames(fd7)[
    rownames(fd7) == "SERPINE1" |
      as.character(fd7$gene_short_name) == "SERPINE1"
  ][1]

  if (!is.na(serpine_gene_id_7)) {
    serp_test_7 <- monocle::differentialGeneTest(
      cds_fig6[
        serpine_gene_id_7,
      ],
      fullModelFormulaStr = "~sm.ns(Pseudotime)",
      cores = 1
    )

    serpine_dynamic_p_7 <- serp_test_7$pval[1]
    serpine_dynamic_q_7 <- serp_test_7$qval[1]
  }
}

p7F_heat <- ggplot2::ggplot(
  module_heat_long_7,
  ggplot2::aes(
    x = Pseudotime,
    y = Gene,
    fill = Zscore
  )
) +
  ggplot2::geom_tile(linewidth = 0) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Expression\n(z-score)"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1.10),
    breaks = c(0, 0.5, 1),
    labels = c("Early", "Mid", "Late"),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    tag = "F",
    title = expression(
      italic(SERPINE1) * "-associated late injury module"
    ),
    subtitle = paste0(
      "SERPINE1 pseudotime: P=",
      format_p_7(serpine_dynamic_p_7),
      ", q=",
      format_p_7(serpine_dynamic_q_7)
    ),
    x = "Endothelial injury pseudotime",
    y = NULL
  ) +
  theme_fig6_7(6.4) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      face = "italic",
      size = 4.9
    ),
    axis.ticks.y = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    legend.position = "bottom"
  )

if (length(serpine_module_y_7) == 1) {
  p7F_heat <- p7F_heat +
    ggplot2::annotate(
      "segment",
      x = 1.005,
      xend = 1.045,
      y = serpine_module_y_7,
      yend = serpine_module_y_7,
      color = "#B2182B",
      linewidth = 0.8,
      arrow = grid::arrow(
        length = grid::unit(1.4, "mm")
      )
    ) +
    ggplot2::annotate(
      "text",
      x = 1.052,
      y = serpine_module_y_7,
      label = "SERPINE1",
      hjust = 0,
      color = "#B2182B",
      fontface = "bold.italic",
      size = 2.25
    )
}

if (
  exists("serpine_module_go", envir = .GlobalEnv, inherits = FALSE)
) {
  go7 <- as.data.frame(
    serpine_module_go
  )
} else {
  go7 <- data.frame()
}

if (nrow(go7) > 0) {
  if (!"Description" %in% colnames(go7)) {
    stop("serpine_module_go has no Description column.")
  }

  if (!"p.adjust" %in% colnames(go7)) {
    stop("serpine_module_go has no p.adjust column.")
  }

  go7 <- go7 %>%
    dplyr::filter(
      is.finite(p.adjust)
    ) %>%
    dplyr::arrange(
      p.adjust,
      pvalue
    ) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::mutate(
      FDRStrength = -log10(pmax(p.adjust, 1e-300)),
      InjuryRelated = grepl(
        paste(
          c(
            "senesc",
            "cell cycle",
            "inflamm",
            "cytokine",
            "chemokine",
            "oxidative",
            "apopt",
            "DNA damage",
            "stress",
            "wound",
            "endothelial"
          ),
          collapse = "|"
        ),
        Description,
        ignore.case = TRUE
      ),
      DescriptionWrap = vapply(
        Description,
        function(x) {
          paste(
            strwrap(x, width = 31),
            collapse = "\n"
          )
        },
        character(1)
      )
    )

  go7$DescriptionWrap <- factor(
    go7$DescriptionWrap,
    levels = rev(go7$DescriptionWrap)
  )

  n_sig_go_7 <- sum(
    go7$p.adjust < 0.05,
    na.rm = TRUE
  )

  p7F_go <- ggplot2::ggplot(
    go7,
    ggplot2::aes(
      x = FDRStrength,
      y = DescriptionWrap,
      fill = InjuryRelated
    )
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_vline(
      xintercept = -log10(0.05),
      linetype = "dashed",
      linewidth = 0.35,
      color = "#7D8890"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "FALSE" = "#8C96C6",
        "TRUE" = "#B2182B"
      ),
      labels = c(
        "FALSE" = "Other GO BP",
        "TRUE" = "Senescence / injury related"
      ),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Functional architecture of the SERPINE1 module",
      subtitle = paste0(
        n_sig_go_7,
        " displayed GO BP terms reached FDR < 0.05"
      ),
      x = expression(-log[10] * "(FDR)"),
      y = NULL
    ) +
    theme_fig6_7(6.1) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 5.0),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom"
    )
} else {
  p7F_go <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = "GO table not available in current session",
      size = 3
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void()
}

p7F <- patchwork::wrap_plots(
  p7F_heat,
  p7F_go,
  ncol = 2,
  widths = c(1.35, 1.00)
)

# =============================================================================

# Save only the retained trajectory panels A-F and their analysis data.
trajectory_panels <- list(A=p7A,B=p7B,C=p7C,D=p7D,E=p7E,F=p7F)
for(nm in names(trajectory_panels))
 save_publication_plot(trajectory_panels[[nm]],file.path(repo_output("Figure5/Trajectory"),paste0("Figure5",nm)),9,6)
saveRDS(list(panels=trajectory_panels),file.path(repo_output("Figure5/Trajectory"),"Figure5_A_F_panels.rds"))
trajectory_data <- list(trajectory_metadata=traj7,
 p53_p21_arrest_genes=p53_obj_7$genes,oxidative_stress_genes=oxidative_obj_7$genes,
 inflammatory_SASP_genes_SERPINE1_excluded=sasp_obj_7$genes,
 ECM_remodeling_genes=ecm_obj_7$genes,barrier_homeostasis_genes=barrier_obj_7$genes,
 late_up_senescence_gene_statistics=late_gene_table_7,
 SERPINE1_module_genes=module_genes_7)
saveRDS(trajectory_data,file.path(repo_output("Figure5/Trajectory"),"Figure5_TRAJECTORY_DATA_ONLY.rds"))
upper <- patchwork::wrap_plots(
 patchwork::wrap_plots(p7A,p7B,p7C,nrow=1),p7D,
 patchwork::wrap_plots(p7E,p7F,nrow=1),ncol=1,heights=c(1,.75,1.1))
save_publication_plot(upper,file.path(repo_output("Figure5/Trajectory"),"Figure5_A_F"),18,14)

}

# =============================================================================
# STAGE figure5_network_master
# =============================================================================
PIPELINE$figure5_network_master <- function() {
# scripts/07_build_network_master.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

# =============================================================================
# Figure 6 - Build ONE compact master RDS
#
# Purpose:
#   Save ONLY data objects needed for Figure 6.
#   Do NOT save Figure 4/5 plot objects.
#   Do NOT rerun Figure 4/5.
#
# Input:
#   sce.Rdata
#   Figure5_analysis_objects_PREMIUM.rds
#
# Output:
#   Figure6_MASTER_DATA_ONLY.rds
# =============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tibble)
})

SCE_RDATA <- repo_input("sce_rdata")
FIG5_RDS <- repo_input("endothelial_bundle")
OUTPUT_RDS <- repo_output("Shared/Figure6_MASTER_DATA_ONLY.rds")


# =============================================================================
# 1. Restore sce
# =============================================================================

if (!file.exists(SCE_RDATA)) {
  stop("Cannot find: ", SCE_RDATA)
}

loaded_names <- load(
  SCE_RDATA,
  envir = .GlobalEnv
)

if (!exists("sce", envir = .GlobalEnv, inherits = FALSE)) {
  stop(
    "sce.Rdata loaded, but object 'sce' was not found. Objects: ",
    paste(loaded_names, collapse = ", ")
  )
}

sce <- get(
  "sce",
  envir = .GlobalEnv
)

if (!inherits(sce, "Seurat")) {
  stop("sce is not a Seurat object.")
}


# =============================================================================
# 2. Restore Figure 5 data objects
# =============================================================================

if (!file.exists(FIG5_RDS)) {
  stop("Cannot find: ", FIG5_RDS)
}

fig5_obj <- readRDS(
  FIG5_RDS
)

if (!is.list(fig5_obj)) {
  stop("Figure5_analysis_objects_PREMIUM.rds is not a list.")
}

required_names <- c(
  "endo_clean",
  "endo_pathway",
  "key_state"
)

missing_names <- setdiff(
  required_names,
  names(fig5_obj)
)

if (length(missing_names) > 0) {
  stop(
    "Figure5 RDS is missing: ",
    paste(missing_names, collapse = ", ")
  )
}

endo_clean <- fig5_obj$endo_clean
endo_pathway <- fig5_obj$endo_pathway
KEY_STATE <- as.character(fig5_obj$key_state)

if (!inherits(endo_clean, "Seurat")) {
  stop("fig5_obj$endo_clean is not a Seurat object.")
}

if (!inherits(endo_pathway, "Seurat")) {
  stop("fig5_obj$endo_pathway is not a Seurat object.")
}


# =============================================================================
# 3. Detect metadata fields
# =============================================================================

detect_meta_field <- function(object, candidates, label) {
  hit <- candidates[
    candidates %in%
      colnames(object@meta.data)
  ]

  if (length(hit) == 0) {
    stop(
      "Cannot find ",
      label,
      " metadata column. Candidates: ",
      paste(candidates, collapse = ", ")
    )
  }

  hit[1]
}

SCE_SAMPLE_COL <- detect_meta_field(
  sce,
  c("sample_id", "orig.ident", "sample", "Sample"),
  "sce sample"
)

SCE_CONDITION_COL <- detect_meta_field(
  sce,
  c("tissue_type", "condition", "Condition", "group", "Group"),
  "sce condition"
)

ENDO_SAMPLE_COL <- detect_meta_field(
  endo_clean,
  c("orig.ident", "sample_id", "sample", "Sample"),
  "endo sample"
)

ENDO_CONDITION_COL <- detect_meta_field(
  endo_clean,
  c("tissue_type", "condition", "Condition", "group", "Group"),
  "endo condition"
)

ENDO_STATE_COL <- detect_meta_field(
  endo_clean,
  c(
    "Endothelial_Annotation",
    "endothelial_annotation",
    "Endothelial_State",
    "celltype",
    "Celltype"
  ),
  "endothelial state"
)


# =============================================================================
# 4. Core checks
# =============================================================================

if (!KEY_STATE %in% unique(as.character(endo_clean@meta.data[[ENDO_STATE_COL]]))) {
  stop("KEY_STATE not found in endo_clean: ", KEY_STATE)
}

common_endo_cells <- intersect(
  colnames(sce),
  colnames(endo_clean)
)

overlap_fraction <- length(common_endo_cells) / max(1, ncol(endo_clean))

if (overlap_fraction < 0.95) {
  warning(
    "Less than 95% of endo_clean barcodes are present in sce."
  )
}


# =============================================================================
# 5. Receiver mapping for Figure 6 CellChat
# =============================================================================

receiver_map <- endo_clean@meta.data %>%
  tibble::rownames_to_column("Cell") %>%
  dplyr::transmute(
    Cell = Cell,
    Sample = as.character(.data[[ENDO_SAMPLE_COL]]),
    Condition = as.character(.data[[ENDO_CONDITION_COL]]),
    Endothelial_State = as.character(.data[[ENDO_STATE_COL]]),
    Is_Key_State = Endothelial_State == KEY_STATE
  )


# =============================================================================
# 6. Keep ONLY useful Figure 5 analysis data
#
# No panels, no ggplot objects, no patchwork objects.
# Missing elements are stored as NULL.
# =============================================================================

fig5_data <- list(
  key_state_DE = fig5_obj[["key_state_DE"]],
  key_state_KEGG = fig5_obj[["key_state_KEGG"]],
  key_state_Hallmark_GSEA = fig5_obj[["key_state_Hallmark_GSEA"]],
  SERPINE1_high_low_DE = fig5_obj[["SERPINE1_high_low_DE"]],
  SERPINE1_high_low_DE_full = fig5_obj[["SERPINE1_high_low_DE_full"]],
  SERPINE1_high_low_Hallmark_GSEA = fig5_obj[["SERPINE1_high_low_Hallmark_GSEA"]],
  SERPINE1_high_enrichment_OR = fig5_obj[["SERPINE1_high_enrichment_OR"]],
  SERPINE1_high_enrichment_OR_full = fig5_obj[["SERPINE1_high_enrichment_OR_full"]]
)


# =============================================================================
# 7. Build compact master object
# =============================================================================

Figure6_MASTER <- list(

  core = list(
    sce = sce,
    endo_clean = endo_clean,
    endo_pathway = endo_pathway,
    key_state = KEY_STATE
  ),

  metadata = list(
    sce_sample_col = SCE_SAMPLE_COL,
    sce_condition_col = SCE_CONDITION_COL,
    endothelial_sample_col = ENDO_SAMPLE_COL,
    endothelial_condition_col = ENDO_CONDITION_COL,
    endothelial_state_col = ENDO_STATE_COL,
    common_endothelial_cells = common_endo_cells,
    endothelial_barcode_overlap_fraction = overlap_fraction,
    receiver_map = receiver_map
  ),

  figure5_data = fig5_data,

  figure6_plan = list(
    trajectory_object = "core$endo_clean",
    pathway_object = "core$endo_pathway",
    communication_object = "core$sce",
    receiver_state = KEY_STATE,
    trajectory_states = c(
      "Homeostatic capillary EC",
      "Transitional capillary EC",
      "Stress-responsive capillary EC"
    ),
    condition_levels = c("CIT", "HR"),
    statistical_unit = "biological sample"
  )
)


# =============================================================================
# 8. Save ONE RDS
# =============================================================================

saveRDS(
  Figure6_MASTER,
  OUTPUT_RDS,
  compress = FALSE
)

assign(
  "Figure6_MASTER",
  Figure6_MASTER,
  envir = .GlobalEnv
)


# =============================================================================
# 9. Summary
# =============================================================================

cat(
  "\n============================================================\n",
  "Figure6 data-only master RDS created\n",
  "============================================================\n",
  "Output: ", OUTPUT_RDS, "\n",
  "sce cells: ", ncol(sce), "\n",
  "endo_clean cells: ", ncol(endo_clean), "\n",
  "key state: ", KEY_STATE, "\n",
  "barcode overlap: ", sprintf("%.3f", overlap_fraction), "\n",
  "\nNo Figure 4/5 plot objects were saved.\n",
  sep = ""
)

cat(
  "\nRestore later with:\n",
  "x <- readRDS('", OUTPUT_RDS, "')\n",
  "sce <- x$core$sce\n",
  "endo_clean <- x$core$endo_clean\n",
  "endo_pathway <- x$core$endo_pathway\n",
  "KEY_STATE <- x$core$key_state\n",
  sep = ""
)




}

# =============================================================================
# STAGE figure5_cellchat_GJ
# =============================================================================
PIPELINE$figure5_cellchat_GJ <- function() {
# scripts/07_cellchat_networks.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.
# =============================================================================
# Reperfusion-associated signaling to SERPINE1-high stress-responsive capillary EC
#
# IMPORTANT:
#   1) Run in a FRESH R session using the MAIN library (igraph >= 2.0).
#   2) Do NOT load Monocle2 / the separate igraph 1.6.0 library in this session.
#   3) This script reads Figure5_MASTER_DATA_ONLY.rds directly.
#
# Main panels:
#   H  Global cell-cell communication changes: CIT vs HR
#   I  Sender-cell changes toward Stress-responsive capillary EC
#   J  Key ligand-receptor pairs targeting Stress-responsive capillary EC
#   K  HR-enhanced sender -> ligand-receptor -> receiver routes
#   L  Biological-sample validation of selected signaling axes
#
# Statistical design:
#   - CellChat is used for condition-level hypothesis generation.
#   - Biological samples are used for Panel L validation.
#   - CellChat population.size = FALSE to reduce abundance-driven effects.
#
# Output:
#   30_Figure5_CellChat_H_L/
# =============================================================================

options(stringsAsFactors = FALSE)
set.seed(1234)

# -----------------------------------------------------------------------------
# 0. Library path / environment
# -----------------------------------------------------------------------------
MAIN_LIB <- ""
if (dir.exists(MAIN_LIB)) {
  .libPaths(c(MAIN_LIB, setdiff(.libPaths(), MAIN_LIB)))
}

if ("igraph" %in% loadedNamespaces()) {
  if (utils::packageVersion("igraph") < "2.0.0") {
    stop(
      "igraph < 2.0 is already loaded. Restart R, put the configured project library first, ",
      "and rerun this CellChat script."
    )
  }
}

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(ragg)
})

cat(
  "\n============================================================\n",
  "Figure 5 G-J CellChat environment\n",
  "============================================================\n",
  "R: ", R.version.string, "\n",
  "CellChat: ", as.character(utils::packageVersion("CellChat")), "\n",
  "igraph: ", as.character(utils::packageVersion("igraph")), "\n",
  "Seurat: ", as.character(utils::packageVersion("Seurat")), "\n",
  sep = ""
)

if (utils::packageVersion("igraph") < "2.0.0") {
  stop("This CellChat session should use the main igraph >= 2.0 library.")
}

# -----------------------------------------------------------------------------
# 1. Paths / parameters
# -----------------------------------------------------------------------------
MASTER_RDS <- repo_input("network_master")

OUTDIR <- repo_output("Figure5/CellChat")
PLOT_DIR <- file.path(OUTDIR, "01_panels")
TABLE_DIR <- file.path(OUTDIR, "02_tables")
OBJECT_DIR <- file.path(OUTDIR, "03_objects")
FINAL_DIR <- file.path(OUTDIR, "04_final")

invisible(lapply(
  c(PLOT_DIR, TABLE_DIR, OBJECT_DIR, FINAL_DIR),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

if (!file.exists(MASTER_RDS)) {
  stop("Cannot find ", MASTER_RDS)
}

RECEIVER <- "Stress-responsive capillary EC"
MIN_CELLS_PER_GROUP <- 20
MIN_CELLS_PER_SAMPLE_GROUP <- 10
MAX_CELLS_PER_SAMPLE_GROUP <- 500
MAX_GLOBAL_EDGES_PER_CONDITION <- 35
TOP_SENDERS_I <- 8
TOP_LR_J <- 10
TOP_ROUTES_K <- 6
TOP_AXES_L <- 4
RANDOM_SEED <- 1234

CONDITION_COLORS <- c(
  "CIT" = "#2166AC",
  "HR" = "#D73027"
)

DELTA_COLORS <- c(
  "Decrease" = "#2166AC",
  "Increase" = "#D73027"
)

# Publication-style cell-population colors used in Panel H.
CELLCHAT_GROUP_COLORS <- c(
  "Stress-responsive capillary EC" = "#F2B701",
  "Other endothelial EC" = "#46A6A6",
  "Macrophage" = "#D95F5F",
  "Monocyte" = "#E68A5B",
  "Neutrophil" = "#8C78B8",
  "NK cell" = "#9A5AB4",
  "T cell" = "#59A14F",
  "B cell" = "#4E79A7",
  "Dendritic cell" = "#2F80A3",
  "Mast cell" = "#D77FB3",
  "AT1" = "#76B7B2",
  "AT2" = "#4C9ED9",
  "Fibroblast" = "#B07AA1",
  "Pericyte/SMC" = "#8A9A5B",
  "Other epithelial" = "#C49A6C",
  "Other" = "#9CA3AF"
)

# Condition-specific edge colors used in the global network.
GLOBAL_EDGE_COLORS <- c(
  "CIT" = "#4F81BD",
  "HR" = "#D95B43"
)

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
find_first_column <- function(meta, candidates, object_name = "metadata") {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) < 1) {
    stop(
      "Could not find a required column in ", object_name, ". Tried: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[1]
}

format_p <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

zscore_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (sum(is.finite(x)) < 2 || !is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

theme_fig6 <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 1.4,
        hjust = 0.5,
        color = "#202A33"
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 0.4,
        hjust = 0.5,
        color = "#596874"
      ),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = base_size + 3,
        color = "#202A33"
      ),
      plot.tag.position = c(0.01, 0.99),
      axis.title = ggplot2::element_text(face = "bold", color = "#202A33"),
      axis.text = ggplot2::element_text(color = "#202A33"),
      panel.grid.major = ggplot2::element_line(color = "#E8EDF0", linewidth = 0.28),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, color = "#CBD4DA", linewidth = 0.45),
      strip.background = ggplot2::element_rect(fill = "#F1F5F7", color = "#CBD4DA"),
      strip.text = ggplot2::element_text(face = "bold", color = "#202A33"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

save_panel <- function(p, filename, width, height, dpi = 600) {
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

  ggplot2::ggsave(
    filename = file.path(PLOT_DIR, paste0(filename, ".pdf")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = pdf_dev,
    bg = "white",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    filename = file.path(PLOT_DIR, paste0(filename, ".png")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = ragg::agg_png,
    bg = "white",
    limitsize = FALSE
  )
}

collapse_major_celltype <- function(x) {
  x0 <- as.character(x)
  xl <- tolower(x0)

  out <- rep("Other", length(xl))

  out[grepl("macroph", xl)] <- "Macrophage"
  out[grepl("monocyte|mono\\b", xl)] <- "Monocyte"
  out[grepl("neutroph", xl)] <- "Neutrophil"
  out[grepl("fibro", xl)] <- "Fibroblast"
  out[grepl("alveolar.*type.?2|alveolar.*ii|\\bat2\\b", xl)] <- "AT2"
  out[grepl("alveolar.*type.?1|alveolar.*i|\\bat1\\b", xl)] <- "AT1"
  out[grepl("dendritic|\\bdc\\b", xl)] <- "Dendritic cell"
  out[grepl("mast", xl)] <- "Mast cell"
  out[grepl("natural killer|\\bnk\\b", xl)] <- "NK cell"
  out[grepl("\\bb cell\\b|b-cell|b lymph|plasma", xl)] <- "B cell"
  out[grepl("\\bt cell\\b|t-cell|t lymph|cd4|cd8", xl)] <- "T cell"
  out[grepl("pericyte|smooth muscle|\\bsmc\\b", xl)] <- "Pericyte/SMC"
  out[grepl("epithelial|ciliated|club|basal", xl)] <- "Other epithelial"
  out[grepl("endothelial|\\bec\\b", xl)] <- "Other endothelial EC"

  out
}

wrap_group_label <- function(x) {
  out <- x
  out[out == "Stress-responsive capillary EC"] <- "Stress-responsive\ncapillary EC"
  out[out == "Other endothelial EC"] <- "Other\nendothelial EC"
  out[out == "Dendritic cell"] <- "Dendritic\ncell"
  out[out == "Other epithelial"] <- "Other\nepithelial"
  out[out == "Pericyte/SMC"] <- "Pericyte/\nSMC"
  out
}

sample_stratified_downsample <- function(meta, max_per_sample_group = 500, seed = 1234) {
  set.seed(seed)

  meta %>%
    dplyr::group_by(Condition, Sample, CellChat_Group) %>%
    dplyr::group_modify(
      ~ {
        n_take <- min(nrow(.x), max_per_sample_group)
        .x[sample(seq_len(nrow(.x)), n_take), , drop = FALSE]
      }
    ) %>%
    dplyr::ungroup()
}

build_cellchat <- function(data_input, meta_input, group_levels, db_use) {
  meta_input$CellChat_Group <- factor(
    as.character(meta_input$CellChat_Group),
    levels = group_levels
  )

  cellchat <- CellChat::createCellChat(object = data_input)
  cellchat <- CellChat::addMeta(cellchat, meta = meta_input)
  cellchat <- CellChat::setIdent(cellchat, ident.use = "CellChat_Group")

  cellchat@DB <- db_use
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)

  cellchat <- CellChat::computeCommunProb(
    cellchat,
    raw.use = TRUE,
    population.size = FALSE
  )

  cellchat <- CellChat::filterCommunication(
    cellchat,
    min.cells = MIN_CELLS_PER_GROUP
  )

  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  cellchat <- CellChat::netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )

  cellchat
}

prepare_comm_table <- function(cellchat, condition_name) {
  x <- CellChat::subsetCommunication(
    cellchat,
    slot.name = "net",
    thresh = 0.05
  )

  if (nrow(x) == 0) {
    return(data.frame())
  }

  required_cols <- c("source", "target", "ligand", "receptor", "prob", "pval")
  missing_cols <- setdiff(required_cols, colnames(x))
  if (length(missing_cols) > 0) {
    stop(
      "subsetCommunication() is missing expected columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!"interaction_name" %in% colnames(x)) {
    x$interaction_name <- paste(x$ligand, x$receptor, sep = "_")
  }

  if (!"pathway_name" %in% colnames(x)) {
    x$pathway_name <- NA_character_
  }

  x %>%
    dplyr::mutate(
      source = as.character(source),
      target = as.character(target),
      ligand = as.character(ligand),
      receptor = as.character(receptor),
      interaction_name = as.character(interaction_name),
      pathway_name = as.character(pathway_name),
      Condition = condition_name,
      PairID = paste(source, target, ligand, receptor, sep = "|")
    ) %>%
    dplyr::group_by(PairID, source, target, ligand, receptor) %>%
    dplyr::summarise(
      interaction_name = dplyr::first(interaction_name),
      pathway_name = paste(unique(stats::na.omit(pathway_name)), collapse = ";"),
      prob = max(prob, na.rm = TRUE),
      pval = min(pval, na.rm = TRUE),
      Condition = dplyr::first(Condition),
      .groups = "drop"
    )
}

get_net_edge_table <- function(cellchat, condition_name) {
  mat <- cellchat@net$weight
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("source", "target", "weight")

  df %>%
    dplyr::mutate(
      source = as.character(source),
      target = as.character(target),
      weight = as.numeric(weight),
      Condition = condition_name
    ) %>%
    dplyr::filter(is.finite(weight), weight > 0)
}

mean_complex_expression <- function(expr_mat, feature_name, cells) {
  if (length(cells) < 1 || is.na(feature_name) || !nzchar(feature_name)) {
    return(NA_real_)
  }

  components <- unique(unlist(strsplit(feature_name, "[_+&:]")))
  components <- components[nzchar(components)]
  components <- intersect(components, rownames(expr_mat))

  if (length(components) < 1) {
    return(NA_real_)
  }

  gene_means <- Matrix::rowMeans(
    expr_mat[components, cells, drop = FALSE]
  )

  mean(gene_means, na.rm = TRUE)
}

safe_sample_association <- function(df) {
  df <- df %>%
    dplyr::filter(
      is.finite(InteractionScore),
      is.finite(ReceiverSERPINE1),
      !is.na(Condition)
    )

  rho <- NA_real_
  rho_p <- NA_real_
  adjusted_beta <- NA_real_
  adjusted_p <- NA_real_

  if (
    nrow(df) >= 6 &&
      length(unique(df$InteractionScore)) >= 3 &&
      length(unique(df$ReceiverSERPINE1)) >= 3
  ) {
    ct <- suppressWarnings(
      stats::cor.test(
        df$InteractionScore,
        df$ReceiverSERPINE1,
        method = "spearman",
        exact = FALSE
      )
    )
    rho <- unname(ct$estimate)
    rho_p <- ct$p.value
  }

  if (
    nrow(df) >= 8 &&
      stats::sd(df$InteractionScore) > 0 &&
      stats::sd(df$ReceiverSERPINE1) > 0 &&
      length(unique(df$Condition)) == 2
  ) {
    df$Interaction_z <- zscore_safe(df$InteractionScore)
    df$Receiver_z <- zscore_safe(df$ReceiverSERPINE1)

    fit <- stats::lm(
      Receiver_z ~ Interaction_z + Condition,
      data = df
    )

    coef_tab <- summary(fit)$coefficients

    if ("Interaction_z" %in% rownames(coef_tab)) {
      adjusted_beta <- coef_tab["Interaction_z", "Estimate"]
      adjusted_p <- coef_tab["Interaction_z", "Pr(>|t|)"]
    }
  }

  data.frame(
    N = nrow(df),
    Spearman_rho = rho,
    Spearman_P = rho_p,
    Adjusted_Beta = adjusted_beta,
    Adjusted_P = adjusted_p,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 3. Load master objects and map endothelial receiver back to whole lung
# -----------------------------------------------------------------------------
master <- readRDS(MASTER_RDS)

sce <- master$core$sce
endo_clean <- master$core$endo_clean
key_state <- as.character(master$core$key_state)

if (!inherits(sce, "Seurat")) stop("master$core$sce is not a Seurat object.")
if (!inherits(endo_clean, "Seurat")) stop("master$core$endo_clean is not a Seurat object.")

if (!identical(key_state, RECEIVER)) {
  warning(
    "Master key_state is '", key_state,
    "' but this script receiver is '", RECEIVER, "'."
  )
}

Seurat::DefaultAssay(sce) <- "RNA"
Seurat::DefaultAssay(endo_clean) <- "RNA"

sce_meta <- sce@meta.data
endo_meta <- endo_clean@meta.data

sce_celltype_col <- find_first_column(
  sce_meta,
  c(
    "celltype", "CellType", "Cell_Type", "cell_type", "Celltype",
    "Annotation", "annotation", "celltype_final", "cell.type"
  ),
  "whole-lung sce metadata"
)

sce_condition_col <- find_first_column(
  sce_meta,
  c("tissue_type", "Condition", "condition", "group", "Group"),
  "whole-lung sce metadata"
)

sce_sample_col <- find_first_column(
  sce_meta,
  c("orig.ident", "Sample", "sample", "sample_id", "SampleID"),
  "whole-lung sce metadata"
)

endo_annotation_col <- find_first_column(
  endo_meta,
  c(
    "Endothelial_Annotation", "EndothelialAnnotation",
    "endo_annotation", "endo_state", "Annotation", "annotation"
  ),
  "endo_clean metadata"
)

endo_state_map <- as.character(endo_meta[[endo_annotation_col]])
names(endo_state_map) <- rownames(endo_meta)

meta_full <- data.frame(
  Cell = colnames(sce),
  Sample = as.character(sce_meta[[sce_sample_col]]),
  Condition = as.character(sce_meta[[sce_condition_col]]),
  OriginalCellType = as.character(sce_meta[[sce_celltype_col]]),
  stringsAsFactors = FALSE
)

meta_full$Condition <- factor(meta_full$Condition, levels = c("CIT", "HR"))
meta_full$CellChat_Group <- collapse_major_celltype(meta_full$OriginalCellType)

mapped_endo <- endo_state_map[meta_full$Cell]

meta_full$CellChat_Group[
  !is.na(mapped_endo) & mapped_endo == RECEIVER
] <- RECEIVER

meta_full$CellChat_Group[
  !is.na(mapped_endo) & mapped_endo != RECEIVER
] <- "Other endothelial EC"

# Any remaining broad endothelial cells are retained as Other endothelial EC.
meta_full$CellChat_Group[
  grepl("endothelial|\\bec\\b", tolower(meta_full$OriginalCellType)) &
    meta_full$CellChat_Group != RECEIVER
] <- "Other endothelial EC"

meta_full <- meta_full %>%
  dplyr::filter(
    Condition %in% c("CIT", "HR"),
    !is.na(Sample),
    nzchar(Sample),
    !is.na(CellChat_Group),
    CellChat_Group != "Other"
  )

# Keep groups adequately represented in BOTH conditions.
group_counts <- meta_full %>%
  dplyr::count(Condition, CellChat_Group, name = "Cells") %>%
  tidyr::complete(
    Condition = factor(c("CIT", "HR"), levels = c("CIT", "HR")),
    CellChat_Group,
    fill = list(Cells = 0)
  )

common_groups <- group_counts %>%
  dplyr::group_by(CellChat_Group) %>%
  dplyr::summarise(
    MinCells = min(Cells),
    .groups = "drop"
  ) %>%
  dplyr::filter(MinCells >= MIN_CELLS_PER_GROUP) %>%
  dplyr::pull(CellChat_Group)

if (!RECEIVER %in% common_groups) {
  stop(
    "Receiver '", RECEIVER,
    "' does not have >= ", MIN_CELLS_PER_GROUP,
    " cells in both CIT and HR."
  )
}

preferred_group_order <- c(
  RECEIVER,
  "Other endothelial EC",
  "Macrophage",
  "Monocyte",
  "Fibroblast",
  "Neutrophil",
  "AT2",
  "AT1",
  "T cell",
  "B cell",
  "NK cell",
  "Dendritic cell",
  "Mast cell",
  "Pericyte/SMC",
  "Other epithelial"
)

group_levels <- c(
  intersect(preferred_group_order, common_groups),
  setdiff(sort(common_groups), preferred_group_order)
)

meta_chat_full <- meta_full %>%
  dplyr::filter(CellChat_Group %in% group_levels)

meta_chat_full$CellChat_Group <- factor(
  meta_chat_full$CellChat_Group,
  levels = group_levels
)

write.csv(
  group_counts,
  file.path(TABLE_DIR, "CellChat_group_counts_full.csv"),
  row.names = FALSE
)

write.csv(
  meta_chat_full,
  file.path(TABLE_DIR, "CellChat_cell_metadata_full.csv"),
  row.names = FALSE
)

cat("\nCellChat groups retained:\n")
print(group_levels)

cat("\nFull cell counts by condition/group:\n")
print(
  meta_chat_full %>%
    dplyr::count(Condition, CellChat_Group)
)

# -----------------------------------------------------------------------------
# 4. Balanced cell selection for CellChat
# -----------------------------------------------------------------------------
meta_chat_ds <- sample_stratified_downsample(
  meta_chat_full,
  max_per_sample_group = MAX_CELLS_PER_SAMPLE_GROUP,
  seed = RANDOM_SEED
)

meta_chat_ds$CellChat_Group <- factor(
  as.character(meta_chat_ds$CellChat_Group),
  levels = group_levels
)

write.csv(
  meta_chat_ds,
  file.path(TABLE_DIR, "CellChat_cell_metadata_downsampled.csv"),
  row.names = FALSE
)

cat("\nDownsampled cell counts by condition/group:\n")
print(
  meta_chat_ds %>%
    dplyr::count(Condition, CellChat_Group)
)

expr_sce <- SeuratObject::GetAssayData(
  sce,
  assay = "RNA",
  slot = "data"
)

# -----------------------------------------------------------------------------
# 5. Run CellChat separately for CIT and HR
# -----------------------------------------------------------------------------
CellChatDB.use <- CellChat::CellChatDB.human

cellchat_list <- list()

for (cond in c("CIT", "HR")) {
  message("\nRunning CellChat for ", cond, "...")

  meta_i <- meta_chat_ds %>%
    dplyr::filter(Condition == cond)

  cells_i <- meta_i$Cell
  meta_i2 <- meta_i %>%
    dplyr::select(Cell, CellChat_Group) %>%
    tibble::column_to_rownames("Cell")

  data_i <- expr_sce[, cells_i, drop = FALSE]

  cellchat_i <- build_cellchat(
    data_input = data_i,
    meta_input = meta_i2,
    group_levels = group_levels,
    db_use = CellChatDB.use
  )

  cellchat_list[[cond]] <- cellchat_i

  saveRDS(
    cellchat_i,
    file.path(OBJECT_DIR, paste0("CellChat_", cond, ".rds")),
    compress = FALSE
  )
}

cellchat_CIT <- cellchat_list[["CIT"]]
cellchat_HR <- cellchat_list[["HR"]]

cellchat_merged <- CellChat::mergeCellChat(
  object.list = cellchat_list,
  add.names = names(cellchat_list)
)

saveRDS(
  cellchat_merged,
  file.path(OBJECT_DIR, "CellChat_CIT_HR_merged.rds"),
  compress = FALSE
)

# -----------------------------------------------------------------------------
# 6. Communication tables and differential receiver-focused analysis
# -----------------------------------------------------------------------------
comm_CIT <- prepare_comm_table(cellchat_CIT, "CIT")
comm_HR <- prepare_comm_table(cellchat_HR, "HR")

write.csv(
  comm_CIT,
  file.path(TABLE_DIR, "CellChat_CIT_significant_interactions.csv"),
  row.names = FALSE
)

write.csv(
  comm_HR,
  file.path(TABLE_DIR, "CellChat_HR_significant_interactions.csv"),
  row.names = FALSE
)

receiver_CIT <- comm_CIT %>%
  dplyr::filter(target == RECEIVER)

receiver_HR <- comm_HR %>%
  dplyr::filter(target == RECEIVER)

pair_diff <- dplyr::full_join(
  receiver_CIT %>%
    dplyr::select(
      PairID, source, target, ligand, receptor,
      interaction_name_CIT = interaction_name,
      pathway_CIT = pathway_name,
      prob_CIT = prob,
      pval_CIT = pval
    ),
  receiver_HR %>%
    dplyr::select(
      PairID, source, target, ligand, receptor,
      interaction_name_HR = interaction_name,
      pathway_HR = pathway_name,
      prob_HR = prob,
      pval_HR = pval
    ),
  by = c("PairID", "source", "target", "ligand", "receptor")
) %>%
  dplyr::mutate(
    prob_CIT = dplyr::coalesce(prob_CIT, 0),
    prob_HR = dplyr::coalesce(prob_HR, 0),
    pval_CIT = dplyr::coalesce(pval_CIT, 1),
    pval_HR = dplyr::coalesce(pval_HR, 1),
    interaction_name = dplyr::coalesce(interaction_name_HR, interaction_name_CIT),
    pathway_name = dplyr::coalesce(pathway_HR, pathway_CIT),
    DeltaProb = prob_HR - prob_CIT,
    Log2FC = log2((prob_HR + 1e-6) / (prob_CIT + 1e-6)),
    IncreasedInHR = DeltaProb > 0,
    LR_Label = paste0(ligand, " -> ", receptor),
    AxisLabel = paste0(source, " | ", ligand, " -> ", receptor)
  ) %>%
  dplyr::arrange(dplyr::desc(DeltaProb))

# All ligand-receptor labels use plain ASCII to avoid Windows/GBK mojibake.
pair_diff$LR_Label <- paste0(pair_diff$ligand, " -> ", pair_diff$receptor)
pair_diff$AxisLabel <- paste0(pair_diff$source, " | ", pair_diff$ligand, " -> ", pair_diff$receptor)

write.csv(
  pair_diff,
  file.path(TABLE_DIR, "Receiver_StressResponsiveEC_LR_differential.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 7. Panel G | global communication network CIT vs HR
#
# Visual refinement:
#   - node fill = cell population
#   - CIT edges = cool blue
#   - HR edges = warm red-orange
#   - Stress-responsive capillary EC = gold + extra halo
# -----------------------------------------------------------------------------

edge_CIT <- get_net_edge_table(cellchat_CIT, "CIT")
edge_HR <- get_net_edge_table(cellchat_HR, "HR")
edge_all <- dplyr::bind_rows(edge_CIT, edge_HR)

edge_plot <- edge_all %>%
  dplyr::group_by(Condition) %>%
  dplyr::slice_max(
    order_by = weight,
    n = MAX_GLOBAL_EDGES_PER_CONDITION,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

n_groups <- length(group_levels)

theta <- seq(
  0,
  2 * pi,
  length.out = n_groups + 1
)[seq_len(n_groups)] +
  pi / 2

node_pos <- data.frame(
  CellChat_Group = group_levels,
  x = cos(theta),
  y = sin(theta),
  stringsAsFactors = FALSE
)

node_pos$Label <- wrap_group_label(
  node_pos$CellChat_Group
)

node_pos$IsReceiver <- node_pos$CellChat_Group == RECEIVER

edge_plot <- edge_plot %>%
  dplyr::left_join(
    node_pos %>%
      dplyr::select(
        CellChat_Group,
        x,
        y
      ),
    by = c(
      "source" = "CellChat_Group"
    )
  ) %>%
  dplyr::rename(
    x_source = x,
    y_source = y
  ) %>%
  dplyr::left_join(
    node_pos %>%
      dplyr::select(
        CellChat_Group,
        x,
        y
      ),
    by = c(
      "target" = "CellChat_Group"
    )
  ) %>%
  dplyr::rename(
    x_target = x,
    y_target = y
  )

node_strength <- edge_all %>%
  dplyr::select(
    Condition,
    source,
    target,
    weight
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      source,
      target
    ),
    names_to = "Direction",
    values_to = "CellChat_Group"
  ) %>%
  dplyr::group_by(
    Condition,
    CellChat_Group
  ) %>%
  dplyr::summarise(
    TotalStrength = sum(
      weight,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

node_plot <- tidyr::crossing(
  Condition = c(
    "CIT",
    "HR"
  ),
  node_pos
) %>%
  dplyr::left_join(
    node_strength,
    by = c(
      "Condition",
      "CellChat_Group"
    )
  ) %>%
  dplyr::mutate(
    TotalStrength = dplyr::coalesce(
      TotalStrength,
      0
    )
  )

# Add safe fallback colors for any unexpected CellChat group.
missing_group_colors <- setdiff(
  unique(
    as.character(
      node_plot$CellChat_Group
    )
  ),
  names(
    CELLCHAT_GROUP_COLORS
  )
)

if (length(missing_group_colors) > 0) {

  extra_group_colors <- grDevices::hcl.colors(
    n = length(missing_group_colors),
    palette = "Dark 3"
  )

  names(extra_group_colors) <- missing_group_colors

  CELLCHAT_GROUP_COLORS <- c(
    CELLCHAT_GROUP_COLORS,
    extra_group_colors
  )
}

network_summary <- data.frame(
  Condition = c(
    "CIT",
    "HR"
  ),
  InteractionCount = c(
    sum(cellchat_CIT@net$count),
    sum(cellchat_HR@net$count)
  ),
  TotalStrength = c(
    sum(cellchat_CIT@net$weight),
    sum(cellchat_HR@net$weight)
  )
)

network_summary$FacetLabel <- paste0(
  network_summary$Condition,
  "\nInteractions = ",
  network_summary$InteractionCount,
  " | Strength = ",
  signif(
    network_summary$TotalStrength,
    3
  )
)

facet_map_H <- setNames(
  network_summary$FacetLabel,
  network_summary$Condition
)

pH <- ggplot2::ggplot() +

  ggplot2::geom_curve(
    data = edge_plot %>%
      dplyr::filter(
        source != target
      ),
    ggplot2::aes(
      x = x_source,
      y = y_source,
      xend = x_target,
      yend = y_target,
      linewidth = weight,
      alpha = weight,
      color = Condition
    ),
    curvature = 0.18,
    arrow = grid::arrow(
      length = grid::unit(
        1.05,
        "mm"
      ),
      type = "closed"
    ),
    lineend = "round"
  ) +

  ggplot2::geom_point(
    data = node_plot,
    ggplot2::aes(
      x = x,
      y = y,
      size = TotalStrength,
      fill = CellChat_Group
    ),
    shape = 21,
    color = "#263238",
    stroke = 0.55
  ) +

  # Additional halo around the receiver.
  ggplot2::geom_point(
    data = node_plot %>%
      dplyr::filter(
        CellChat_Group == RECEIVER
      ),
    ggplot2::aes(
      x = x,
      y = y
    ),
    shape = 21,
    fill = NA,
    color = "#8C6500",
    stroke = 1.15,
    size = 8.0,
    inherit.aes = FALSE
  ) +

  ggplot2::geom_text(
    data = node_plot,
    ggplot2::aes(
      x = 1.20 * x,
      y = 1.20 * y,
      label = Label,
      fontface = ifelse(
        CellChat_Group == RECEIVER,
        "bold",
        "plain"
      )
    ),
    size = 2.35,
    color = "#202A33",
    lineheight = 0.88
  ) +

  ggplot2::facet_wrap(
    ~ Condition,
    nrow = 1,
    labeller = ggplot2::labeller(
      Condition = facet_map_H
    )
  ) +

  ggplot2::scale_color_manual(
    values = GLOBAL_EDGE_COLORS,
    guide = "none"
  ) +

  ggplot2::scale_fill_manual(
    values = CELLCHAT_GROUP_COLORS,
    guide = "none"
  ) +

  ggplot2::scale_size_continuous(
    range = c(
      3.2,
      7.1
    ),
    guide = "none"
  ) +

  ggplot2::scale_linewidth_continuous(
    range = c(
      0.20,
      1.65
    ),
    guide = "none"
  ) +

  ggplot2::scale_alpha_continuous(
    range = c(
      0.14,
      0.72
    ),
    guide = "none"
  ) +

  ggplot2::coord_equal(
    xlim = c(
      -1.45,
      1.45
    ),
    ylim = c(
      -1.42,
      1.42
    ),
    clip = "off"
  ) +

  ggplot2::labs(
    tag = "H",
    title = "Global cell-cell communication changes",
    subtitle = "Node colors denote cell populations; blue and red-orange edges denote CIT and HR networks",
    x = NULL,
    y = NULL
  ) +

  ggplot2::theme_void(
    base_size = 8
  ) +

  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5,
      color = "#202A33"
    ),
    plot.subtitle = ggplot2::element_text(
      hjust = 0.5,
      color = "#596874",
      size = 6.8
    ),
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 7.3,
      color = "#202A33"
    ),
    panel.spacing = grid::unit(
      5,
      "mm"
    ),
    plot.margin = ggplot2::margin(
      8,
      12,
      8,
      12
    )
  )

# -----------------------------------------------------------------------------
# 8. Panel H | sender cells to receiver
# -----------------------------------------------------------------------------

sender_CIT <- receiver_CIT %>%
  dplyr::group_by(
    source
  ) %>%
  dplyr::summarise(
    Strength_CIT = sum(
      prob,
      na.rm = TRUE
    ),
    Interactions_CIT = dplyr::n(),
    .groups = "drop"
  )

sender_HR <- receiver_HR %>%
  dplyr::group_by(
    source
  ) %>%
  dplyr::summarise(
    Strength_HR = sum(
      prob,
      na.rm = TRUE
    ),
    Interactions_HR = dplyr::n(),
    .groups = "drop"
  )

sender_diff <- dplyr::full_join(
  sender_CIT,
  sender_HR,
  by = "source"
) %>%
  dplyr::mutate(
    Strength_CIT = dplyr::coalesce(
      Strength_CIT,
      0
    ),
    Strength_HR = dplyr::coalesce(
      Strength_HR,
      0
    ),
    Interactions_CIT = dplyr::coalesce(
      Interactions_CIT,
      0L
    ),
    Interactions_HR = dplyr::coalesce(
      Interactions_HR,
      0L
    ),
    DeltaStrength = Strength_HR - Strength_CIT,
    Direction = ifelse(
      DeltaStrength >= 0,
      "Increase",
      "Decrease"
    )
  ) %>%
  dplyr::filter(
    source != RECEIVER
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      abs(
        DeltaStrength
      )
    )
  )

write.csv(
  sender_diff,
  file.path(
    TABLE_DIR,
    "Receiver_sender_strength_differential.csv"
  ),
  row.names = FALSE
)

sender_plot <- sender_diff %>%
  dplyr::slice_head(
    n = TOP_SENDERS_I
  )

sender_plot$source <- factor(
  sender_plot$source,
  levels = rev(
    sender_plot$source
  )
)

delta_min_I <- min(
  c(
    0,
    sender_plot$DeltaStrength
  ),
  na.rm = TRUE
)

delta_max_I <- max(
  c(
    0,
    sender_plot$DeltaStrength
  ),
  na.rm = TRUE
)

delta_span_I <- delta_max_I -
  delta_min_I

if (
  !is.finite(delta_span_I) ||
    delta_span_I <= 0
) {
  delta_span_I <- 1
}

sender_plot <- sender_plot %>%
  dplyr::mutate(
    Label = paste0(
      "HR ",
      sprintf(
        "%.2f",
        Strength_HR
      ),
      " / CIT ",
      sprintf(
        "%.2f",
        Strength_CIT
      )
    ),
    LabelX = ifelse(
      DeltaStrength >= 0,
      DeltaStrength +
        0.035 *
          delta_span_I,
      DeltaStrength -
        0.035 *
          delta_span_I
    ),
    LabelHjust = ifelse(
      DeltaStrength >= 0,
      0,
      1
    )
  )

x_limits_I <- c(
  delta_min_I -
    0.22 *
      delta_span_I,
  delta_max_I +
    0.42 *
      delta_span_I
)

pI <- ggplot2::ggplot(
  sender_plot,
  ggplot2::aes(
    x = DeltaStrength,
    y = source,
    fill = Direction
  )
) +

  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "#8C969D",
    linewidth = 0.4
  ) +

  ggplot2::geom_col(
    width = 0.72
  ) +

  ggplot2::geom_text(
    ggplot2::aes(
      x = LabelX,
      label = Label,
      hjust = LabelHjust
    ),
    size = 2.05,
    color = "#202A33"
  ) +

  ggplot2::scale_fill_manual(
    values = DELTA_COLORS,
    guide = "none"
  ) +

  ggplot2::scale_x_continuous(
    limits = x_limits_I,
    expand = ggplot2::expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  ggplot2::labs(
    tag = "I",
    title = paste0(
      "Incoming signals to ",
      RECEIVER
    ),
    subtitle = "Change in aggregate sender strength: HR minus CIT",
    x = "Delta incoming communication strength",
    y = NULL
  ) +

  theme_fig6(
    7
  ) +

  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(
      6,
      8,
      6,
      8
    )
  )

# -----------------------------------------------------------------------------
# 9. Panel I | key ligand-receptor pairs
# -----------------------------------------------------------------------------

top_lr <- pair_diff %>%
  dplyr::filter(
    source != RECEIVER,
    prob_HR > 0,
    pval_HR < 0.05,
    DeltaProb > 0
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      DeltaProb
    ),
    dplyr::desc(
      prob_HR
    )
  ) %>%
  dplyr::slice_head(
    n = TOP_LR_J
  ) %>%
  dplyr::mutate(
    LR_Label = paste0(
      ligand,
      " -> ",
      receptor
    ),
    AxisLabel = paste0(
      source,
      " | ",
      ligand,
      " -> ",
      receptor
    ),
    AxisLabel_J = paste0(
      source,
      "\n",
      ligand,
      " -> ",
      receptor
    )
  )

if (nrow(top_lr) < 1) {
  stop(
    "No HR-enhanced significant ligand-receptor pairs were found for receiver ",
    RECEIVER,
    "."
  )
}

write.csv(
  top_lr,
  file.path(
    TABLE_DIR,
    "Top_HR_enhanced_receiver_LR_pairs_ENGLISH.csv"
  ),
  row.names = FALSE
)

lr_long <- dplyr::bind_rows(
  top_lr %>%
    dplyr::transmute(
      PairID,
      AxisLabel_J,
      Condition = "CIT",
      prob = prob_CIT,
      pval = pval_CIT
    ),
  top_lr %>%
    dplyr::transmute(
      PairID,
      AxisLabel_J,
      Condition = "HR",
      prob = prob_HR,
      pval = pval_HR
    )
) %>%
  dplyr::mutate(
    Condition = factor(
      Condition,
      levels = c(
        "CIT",
        "HR"
      )
    ),
    NegLog10P = -log10(
      pmax(
        pval,
        1e-300
      )
    )
  )

axis_order_J <- rev(
  top_lr$AxisLabel_J
)

lr_long$AxisLabel_J <- factor(
  lr_long$AxisLabel_J,
  levels = axis_order_J
)

pJ <- ggplot2::ggplot(
  lr_long,
  ggplot2::aes(
    x = Condition,
    y = AxisLabel_J,
    size = prob,
    fill = Condition,
    alpha = NegLog10P
  )
) +

  ggplot2::geom_point(
    shape = 21,
    color = "#37474F",
    stroke = 0.35
  ) +

  ggplot2::scale_fill_manual(
    values = CONDITION_COLORS
  ) +

  ggplot2::scale_size_continuous(
    range = c(
      1.5,
      6.5
    ),
    name = "Communication\nprobability"
  ) +

  ggplot2::scale_alpha_continuous(
    range = c(
      0.25,
      1
    ),
    name = "-log10(P)"
  ) +

  ggplot2::labs(
    tag = "J",
    title = "Key ligand-receptor pairs targeting the injury EC state",
    subtitle = "Top HR-enhanced incoming interactions",
    x = NULL,
    y = NULL
  ) +

  theme_fig6(
    6.6
  ) +

  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      size = 5.4,
      lineheight = 0.92
    ),
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "right",
    plot.margin = ggplot2::margin(
      6,
      8,
      6,
      8
    )
  )

# -----------------------------------------------------------------------------
# 10. Panel J | top HR-enhanced sender-LR-receiver routes
# -----------------------------------------------------------------------------

top_routes <- top_lr %>%
  dplyr::slice_head(
    n = TOP_ROUTES_K
  ) %>%
  dplyr::mutate(
    RouteOrder = rev(
      seq_len(
        dplyr::n()
      )
    ),
    RouteLabel = paste0(
      ligand,
      " -> ",
      receptor
    )
  )

write.csv(
  top_routes,
  file.path(
    TABLE_DIR,
    "PanelK_top_HR_enhanced_routes_ENGLISH.csv"
  ),
  row.names = FALSE
)

route_edges <- dplyr::bind_rows(
  top_routes %>%
    dplyr::transmute(
      x = 0.12,
      xend = 0.86,
      y = RouteOrder,
      yend = RouteOrder,
      DeltaProb,
      Segment = "Sender-to-LR"
    ),
  top_routes %>%
    dplyr::transmute(
      x = 1.14,
      xend = 1.88,
      y = RouteOrder,
      yend = RouteOrder,
      DeltaProb,
      Segment = "LR-to-receiver"
    )
)

pK <- ggplot2::ggplot() +

  ggplot2::geom_segment(
    data = route_edges,
    ggplot2::aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend,
      linewidth = DeltaProb,
      color = DeltaProb
    ),
    arrow = grid::arrow(
      length = grid::unit(
        1.3,
        "mm"
      ),
      type = "closed"
    )
  ) +

  ggplot2::geom_label(
    data = top_routes,
    ggplot2::aes(
      x = 0,
      y = RouteOrder,
      label = source
    ),
    fill = "#E8EEF2",
    color = "#263238",
    label.size = 0.22,
    size = 2.35
  ) +

  ggplot2::geom_label(
    data = top_routes,
    ggplot2::aes(
      x = 1,
      y = RouteOrder,
      label = RouteLabel
    ),
    fill = "#FFF3E0",
    color = "#7A3E00",
    label.size = 0.22,
    size = 2.20
  ) +

  ggplot2::geom_label(
    data = data.frame(
      x = 2,
      y = mean(
        top_routes$RouteOrder
      ),
      label = RECEIVER
    ),
    ggplot2::aes(
      x = x,
      y = y,
      label = label
    ),
    fill = "#FDE2E2",
    color = "#8B1E1E",
    fontface = "bold",
    label.size = 0.30,
    size = 2.45
  ) +

  ggplot2::scale_color_gradient(
    low = "#F4A582",
    high = "#B2182B",
    name = "HR-CIT\nDelta probability"
  ) +

  ggplot2::scale_linewidth_continuous(
    range = c(
      0.5,
      2.0
    ),
    guide = "none"
  ) +

  ggplot2::annotate(
    "text",
    x = 0,
    y = max(top_routes$RouteOrder) + 0.75,
    label = "Sender",
    fontface = "bold",
    size = 3
  ) +

  ggplot2::annotate(
    "text",
    x = 1,
    y = max(top_routes$RouteOrder) + 0.75,
    label = "Ligand -> receptor",
    fontface = "bold",
    size = 3
  ) +

  ggplot2::annotate(
    "text",
    x = 2,
    y = max(top_routes$RouteOrder) + 0.75,
    label = "Receiver",
    fontface = "bold",
    size = 3
  ) +

  ggplot2::coord_cartesian(
    xlim = c(
      -0.28,
      2.42
    ),
    ylim = c(
      0.3,
      max(top_routes$RouteOrder) + 1.0
    ),
    clip = "off"
  ) +

  ggplot2::labs(
    tag = "K",
    title = "HR-enhanced signaling routes converging on the SERPINE1-high injury EC state",
    subtitle = "Top incoming ligand-receptor axes ranked by communication-probability increase",
    x = NULL,
    y = NULL
  ) +

  ggplot2::theme_void(
    base_size = 8
  ) +

  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5,
      color = "#202A33"
    ),
    plot.subtitle = ggplot2::element_text(
      hjust = 0.5,
      color = "#596874"
    ),
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "right",
    plot.margin = ggplot2::margin(
      8,
      22,
      8,
      18
    )
  )

# -----------------------------------------------------------------------------

network_panels <- list(G=pH,H=pI,I=pJ,J=pK)
for (nm in names(network_panels)) save_publication_plot(network_panels[[nm]],file.path(repo_output("Figure5/CellChat"),paste0("Figure5",nm)),9,7)
saveRDS(network_panels,file.path(repo_output("Figure5/CellChat"),"Figure5_G_J_panels.rds"))

}

# =============================================================================
# STAGE figure5_paired_proxy_K
# =============================================================================
PIPELINE$figure5_paired_proxy_K <- function() {
# Figure 5K: expression-derived ligand-receptor proxies, not CellChat probabilities.
require_packages(c("SeuratObject","Matrix","ggplot2"))
sce<-load_sce();md<-sce@meta.data
for(k in c("orig.ident","celltype"))if(!k%in%names(md))stop("Missing atlas metadata: ",k)
e<-readRDS(repo_input("endothelial_bundle"))
if(is.null(e$endo_clean))stop("Expected endo_clean in the final endothelial bundle.")
endo<-e$endo_clean;em<-endo@meta.data
state_col<-CONFIG$endothelial_state_column
if(!state_col%in%names(em))stop("Configure the actual endothelial-state column.")
receiver<-rownames(em)[as.character(em[[state_col]])==CONFIG$stress_state]
expr<-tryCatch(SeuratObject::GetAssayData(sce,assay="RNA",layer="data"),error=function(e)SeuratObject::GetAssayData(sce,assay="RNA",slot="data"))
axes<-data.frame(Sender=c(CONFIG$nk_label,CONFIG$nk_label,CONFIG$monocyte_label,CONFIG$macrophage_label),Ligand=c("GZMA","ITGB2","ITGB2","ITGB2"),Receptor=c("F2RL3","ICAM1","ICAM1","ICAM1"))
map<-read.csv(file.path(REPO_ROOT,"metadata/single_cell_sample_map.csv"))
rows<-list()
for(i in seq_len(nrow(axes)))for(j in seq_len(nrow(map))) {
 a<-axes[i,];sample<-map$Sample[j]
 send<-rownames(md)[md$orig.ident==sample & md$celltype==a$Sender]
 recv<-intersect(receiver,rownames(md)[md$orig.ident==sample])
 if(!all(c(a$Ligand,a$Receptor)%in%rownames(expr)))stop("Ligand or receptor absent from RNA assay.")
 rows[[length(rows)+1L]]<-data.frame(Axis=paste(a$Sender,paste(a$Ligand,a$Receptor,sep="-")),Sample=sample,Donor=map$Donor[j],Condition=map$Condition[j],Sender_cells=length(send),Receiver_cells=length(recv),Score=if(length(send)>=10L&&length(recv)>=10L)sqrt(mean(expr[a$Ligand,send])*mean(expr[a$Receptor,recv])) else NA_real_)
}
d<-do.call(rbind,rows)
# Enumerate all sign allocations of the nonzero paired signed ranks. This
# preserves the exact paired signed-rank method for L; this exact paired test is retained for the final proxy panel.
exact_signed_rank<-function(delta) {
 delta<-delta[is.finite(delta)&delta!=0];n<-length(delta)
 if(!n)return(1)
 ranks<-rank(abs(delta),ties.method="average");observed<-abs(sum(sign(delta)*ranks))
 signs<-as.matrix(expand.grid(rep(list(c(-1,1)),n)))
 mean(abs(drop(signs%*%ranks))>=observed-1e-10)
}
st<-do.call(rbind,lapply(unique(d$Axis),function(a) {
 z<-d[d$Axis==a,];w<-reshape(z[c("Donor","Condition","Score")],idvar="Donor",timevar="Condition",direction="wide")
 w<-w[is.finite(w$Score.CIT)&is.finite(w$Score.HR),]
 data.frame(Axis=a,N_pairs=nrow(w),P=if(nrow(w))exact_signed_rank(w$Score.HR-w$Score.CIT) else NA_real_)
}))
st$BH_FDR<-p.adjust(st$P,method="BH",n=4L)
d$Condition<-factor(d$Condition,levels=c("CIT","HR"))
st$Label<-paste0("P = ",format.pval(st$P,digits=3),"\nBH-FDR = ",format.pval(st$BH_FDR,digits=3))
p<-ggplot2::ggplot(d,ggplot2::aes(Condition,Score))+ggplot2::geom_line(ggplot2::aes(group=Donor),colour="grey65")+ggplot2::geom_point(ggplot2::aes(colour=Condition),size=2.5)+ggplot2::facet_wrap(~Axis,scales="free_y",nrow=1)+ggplot2::geom_text(data=st,ggplot2::aes(x=1.5,y=Inf,label=Label),inherit.aes=FALSE,vjust=1.2,size=3)+ggplot2::scale_colour_manual(values=c(CIT="#468CCB",HR="#CC4F4E"))+ggplot2::theme_classic()+ggplot2::theme(legend.position="none")+ggplot2::labs(x=NULL,y="Expression-derived interaction score")
out<-repo_output("Figure5/Paired_K")
write.csv(d,file.path(out,"Figure5K_sample_scores.csv"),row.names=FALSE);write.csv(st,file.path(out,"Figure5K_paired_statistics.csv"),row.names=FALSE)
save_publication_plot(p,file.path(out,"Figure5K"),15,4.5);saveRDS(p,file.path(out,"Figure5K_panel.rds"))

}

# =============================================================================
# STAGE figure6_perturbation
# =============================================================================
PIPELINE$figure6_perturbation <- function() {
# scripts/08_virtual_perturbation_compute.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

# =============================================================================
# SERPINE1 virtual knockout | OPTIMIZED FINAL
#
# PRIMARY:
#   HR capillary endothelial continuum
#   - Homeostatic capillary EC
#   - Transitional capillary EC
#   - Stress-responsive capillary EC
#
# VALIDATION:
#   HR Stress-responsive capillary EC
#
# Key optimizations
#   1. Fixed data source: Figure6_MASTER_DATA_ONLY.rds
#   2. ALL eligible cells are retained in the input pool
#      (no external downsampling before scTenifoldKnk)
#   3. nc_nCells = 1000 and nc_nNet = 20
#   4. Primary and validation use INDEPENDENT feature universes
#   5. Each context uses 2500 HVGs + mechanistic anchor genes + SERPINE1
#   6. Ranked GSEA uses the full scTenifoldKnk perturbation Z profile
#      rather than ORA on only a few FDR-significant genes
#   7. Primary-vs-validation comparison is performed only AFTER both GRNs
#      are independently constructed
#   8. Epithelial-marker QC is reported because SFTPC appeared in the
#      previous primary virtual-KO result; no cells/genes are automatically
#      removed by this QC.
#
# IMPORTANT INTERPRETATION
#   diffRegulation$distance / Z / FC describe NETWORK PERTURBATION.
#   They are NOT expression fold changes and must not be described as
#   transcriptional up-regulation or down-regulation.
#
# Recommended:
#   Run in a fresh R session using the main R library, not the Monocle2 library.
# =============================================================================


# =============================================================================
# 0. Packages
# =============================================================================

required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "scTenifoldNet",
  "scTenifoldKnk",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "patchwork",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "ragg",
  "scales"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(scTenifoldNet)
  library(scTenifoldKnk)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ragg)
  library(scales)
})

set.seed(123)


# =============================================================================
# 1. Settings
# =============================================================================

MASTER_RDS <- repo_input("network_master")

TARGET_GENE <- "SERPINE1"
CONDITION_KEEP <- "HR"

CAPILLARY_STATES <- c(
  "Homeostatic capillary EC",
  "Transitional capillary EC",
  "Stress-responsive capillary EC"
)

# Feature selection
N_HVGS_PRIMARY <- 2500
N_HVGS_VALIDATION <- 2500

# scTenifoldKnk network ensemble
NC_NNET <- 20
NC_NCELLS_TARGET <- 1000
TD_K <- 3

# endo_clean has already undergone QC, so internal scTenifoldKnk QC is kept
# deliberately permissive to avoid silently removing eligible cells.
QC_MT_THRESHOLD <- 1
QC_MIN_LSIZE <- 1

# Gene detection thresholds before network construction.
# These are context-specific; no common feature universe is forced.
MIN_DETECTION_FRACTION <- 0.01
MIN_DETECTION_CELLS <- 10

# Statistics
PADJ_CUTOFF <- 0.05
TOP_LABEL_N <- 15
TOP_GENE_N <- 20
TOP_GO_N <- 12

# Ranked GSEA
GSEA_MIN_GS_SIZE_TARGETED <- 4
GSEA_MIN_GS_SIZE_GO <- 10
GSEA_MAX_GS_SIZE <- 500

# Optional extra seed stability analysis.
# FALSE is recommended for the first optimized run because 20-network ensembles
# are already computationally intensive.
RUN_EXTRA_SEED_STABILITY <- FALSE
EXTRA_SEEDS <- c(101, 202, 303)

# Output
OUTDIR <- repo_output("Figure6/Compute")

DIR_QC <- file.path(OUTDIR, "01_QC")
DIR_PRIMARY <- file.path(OUTDIR, "02_PRIMARY_capillary")
DIR_VALIDATION <- file.path(OUTDIR, "03_VALIDATION_stressEC")
DIR_COMPARE <- file.path(OUTDIR, "04_COMPARISON")
DIR_OBJECT <- file.path(OUTDIR, "05_OBJECTS")
DIR_FINAL <- file.path(OUTDIR, "06_FINAL_FIGURE")

invisible(
  lapply(
    c(
      OUTDIR,
      DIR_QC,
      DIR_PRIMARY,
      DIR_VALIDATION,
      DIR_COMPARE,
      DIR_OBJECT,
      DIR_FINAL
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)


# =============================================================================
# 2. Helpers
# =============================================================================

normalize_meta_value <- function(x) {
  trimws(as.character(x))
}


find_meta_col_by_values <- function(meta, values, preferred = character(), label = "metadata") {

  values <- trimws(as.character(values))

  contains_all <- function(nm) {
    vals <- unique(normalize_meta_value(meta[[nm]]))
    all(values %in% vals)
  }

  preferred_present <- preferred[
    preferred %in% colnames(meta)
  ]

  if (length(preferred_present) > 0) {
    preferred_hits <- preferred_present[
      vapply(
        preferred_present,
        contains_all,
        logical(1)
      )
    ]

    if (length(preferred_hits) > 0) {
      return(preferred_hits[1])
    }
  }

  all_hits <- colnames(meta)[
    vapply(
      colnames(meta),
      contains_all,
      logical(1)
    )
  ]

  if (length(all_hits) > 0) {
    return(all_hits[1])
  }

  stop(
    "Could not find a metadata column containing requested ",
    label,
    " values: ",
    paste(values, collapse = ", ")
  )
}


pick_meta_col <- function(meta, candidates, label) {

  hit <- candidates[
    candidates %in% colnames(meta)
  ][1]

  if (is.na(hit)) {
    stop(
      "Could not identify ",
      label,
      " metadata column. Tried: ",
      paste(candidates, collapse = ", ")
    )
  }

  hit
}


get_raw_counts <- function(obj) {

  DefaultAssay(obj) <- "RNA"

  if ("Assay5" %in% class(obj[["RNA"]])) {
    obj <- SeuratObject::JoinLayers(
      obj,
      assay = "RNA"
    )
  }

  mat <- tryCatch(
    {
      Seurat::GetAssayData(
        obj,
        assay = "RNA",
        layer = "counts"
      )
    },
    error = function(e) {
      Seurat::GetAssayData(
        obj,
        assay = "RNA",
        slot = "counts"
      )
    }
  )

  list(
    object = obj,
    counts = mat
  )
}


theme_pub <- function(base_size = 9) {

  ggplot2::theme_bw(
    base_size = base_size
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "sans",
        color = "#202A33"
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = base_size + 1.4
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = base_size - 0.4,
        color = "#66737D"
      ),
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.text = ggplot2::element_text(
        color = "#263540"
      ),
      strip.background = ggplot2::element_rect(
        fill = "#F2F5F7",
        color = "#CBD3D8"
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "#EDF1F4",
        linewidth = 0.3
      ),
      panel.border = ggplot2::element_rect(
        color = "#596771",
        fill = NA,
        linewidth = 0.6
      )
    )
}


save_plot <- function(p, filename, width, height, folder, dpi = 600) {

  pdf_device <- if (
    capabilities("cairo")
  ) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggplot2::ggsave(
    filename = file.path(
      folder,
      paste0(filename, ".pdf")
    ),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = pdf_device,
    bg = "white",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    filename = file.path(
      folder,
      paste0(filename, ".png")
    ),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = ragg::agg_png,
    bg = "white",
    limitsize = FALSE
  )
}


get_variable_genes <- function(obj, nfeatures) {

  obj <- Seurat::NormalizeData(
    obj,
    verbose = FALSE
  )

  obj <- Seurat::FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = FALSE
  )

  Seurat::VariableFeatures(obj)
}


prepare_context_matrix <- function(
    obj,
    counts,
    target_gene,
    anchored_genes,
    n_hvgs,
    min_detection_fraction,
    min_detection_cells,
    context_name
) {

  hvgs <- get_variable_genes(
    obj,
    n_hvgs
  )

  hvgs <- hvgs[
    !grepl(
      "^(MT-|RPL|RPS)",
      hvgs,
      ignore.case = TRUE
    )
  ]

  genes_candidate <- unique(
    c(
      target_gene,
      hvgs,
      anchored_genes
    )
  )

  genes_candidate <- intersect(
    genes_candidate,
    rownames(counts)
  )

  mat0 <- counts[
    genes_candidate,
    ,
    drop = FALSE
  ]

  min_detect <- max(
    min_detection_cells,
    ceiling(
      min_detection_fraction *
        ncol(mat0)
    )
  )

  detected <- Matrix::rowSums(
    mat0 > 0
  )

  keep <- detected >= min_detect

  keep[
    rownames(mat0) == target_gene
  ] <- TRUE

  mat <- mat0[
    keep,
    ,
    drop = FALSE
  ]

  mat <- mat[
    ,
    Matrix::colSums(mat) > 0,
    drop = FALSE
  ]

  feature_qc <- data.frame(
    Context = context_name,
    CandidateGenes = nrow(mat0),
    RetainedGenes = nrow(mat),
    Cells = ncol(mat),
    MinDetectionCells = min_detect,
    stringsAsFactors = FALSE
  )

  list(
    matrix = mat,
    hvgs = hvgs,
    feature_qc = feature_qc
  )
}


prepare_knk_input <- function(mat) {

  knk_version <- utils::packageVersion(
    "scTenifoldKnk"
  )

  dense <- as.matrix(
    mat
  )

  if (knk_version < "1.1.0") {
    out <- as.data.frame(
      dense,
      check.names = FALSE
    )
  } else {
    out <- dense
  }

  out
}


run_scTenifoldKnk_compat <- function(
    input_counts,
    target_gene,
    nc_nNet,
    nc_nCells,
    td_K,
    qc_mtThreshold,
    qc_minLSize,
    seed
) {

  fun <- scTenifoldKnk::scTenifoldKnk

  fml <- names(
    formals(fun)
  )

  args <- list(
    countMatrix = input_counts,
    gKO = target_gene
  )

  if ("nc_nNet" %in% fml) {
    args$nc_nNet <- nc_nNet
  }

  if ("nc_nCells" %in% fml) {
    args$nc_nCells <- nc_nCells
  }

  if ("td_K" %in% fml) {
    args$td_K <- td_K
  }

  if ("qc_mtThreshold" %in% fml) {
    args$qc_mtThreshold <- qc_mtThreshold
  }

  if ("qc_minLSize" %in% fml) {
    args$qc_minLSize <- qc_minLSize
  }

  if ("qc_minLibSize" %in% fml) {
    args$qc_minLibSize <- qc_minLSize
  }

  if ("verbose" %in% fml) {
    args$verbose <- TRUE
  }

  set.seed(seed)

  do.call(
    fun,
    args
  )
}


process_diff_result <- function(result, target_gene, padj_cutoff) {

  if (is.null(result$diffRegulation)) {
    stop(
      "result$diffRegulation was not found."
    )
  }

  df <- result$diffRegulation

  needed <- c(
    "gene",
    "distance",
    "Z",
    "FC",
    "p.value",
    "p.adj"
  )

  missing_cols <- setdiff(
    needed,
    colnames(df)
  )

  if (length(missing_cols) > 0) {
    stop(
      "diffRegulation is missing columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }

  df <- df[
    df$gene != target_gene,
    ,
    drop = FALSE
  ]

  df$distance <- as.numeric(df$distance)
  df$Z <- as.numeric(df$Z)
  df$FC <- as.numeric(df$FC)
  df$p.value <- as.numeric(df$p.value)
  df$p.adj <- as.numeric(df$p.adj)

  df <- df[
    is.finite(df$distance) &
      is.finite(df$Z) &
      is.finite(df$FC) &
      !is.na(df$p.adj),
    ,
    drop = FALSE
  ]

  df$p.adj.plot <- pmax(
    df$p.adj,
    1e-300
  )

  df$negLog10Padj <- -log10(
    df$p.adj.plot
  )

  df$Significance <- ifelse(
    df$p.adj < padj_cutoff,
    "FDR < 0.05",
    "Not significant"
  )

  df <- df[
    order(
      df$p.adj,
      -df$distance
    ),
    ,
    drop = FALSE
  ]

  df
}


make_rank_vector <- function(df) {

  rank_df <- df %>%
    dplyr::select(
      gene,
      Z
    ) %>%
    dplyr::filter(
      is.finite(Z),
      !is.na(gene),
      gene != ""
    ) %>%
    dplyr::group_by(
      gene
    ) %>%
    dplyr::summarise(
      Z = max(
        Z,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  ranks <- rank_df$Z

  names(ranks) <- rank_df$gene

  sort(
    ranks,
    decreasing = TRUE
  )
}


run_targeted_gsea <- function(rank_vector, program_sets) {

  term2gene <- dplyr::bind_rows(
    lapply(
      names(program_sets),
      function(term_i) {
        data.frame(
          term = term_i,
          gene = program_sets[[term_i]],
          stringsAsFactors = FALSE
        )
      }
    )
  )

  term2gene <- term2gene %>%
    dplyr::filter(
      gene %in% names(rank_vector)
    )

  if (nrow(term2gene) == 0) {
    return(
      data.frame()
    )
  }

  gsea_fun <- clusterProfiler::GSEA

  gsea_formals <- names(
    formals(gsea_fun)
  )

  args <- list(
    geneList = rank_vector,
    TERM2GENE = term2gene,
    minGSSize = GSEA_MIN_GS_SIZE_TARGETED,
    maxGSSize = GSEA_MAX_GS_SIZE,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )

  if ("eps" %in% gsea_formals) {
    args$eps <- 0
  }

  if ("seed" %in% gsea_formals) {
    args$seed <- TRUE
  }

  x <- suppressMessages(
    do.call(
      gsea_fun,
      args
    )
  )

  as.data.frame(x)
}


run_go_gsea <- function(rank_vector) {

  symbol_rank_df <- data.frame(
    SYMBOL = names(rank_vector),
    RankScore = as.numeric(rank_vector),
    stringsAsFactors = FALSE
  )

  map_df <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(symbol_rank_df$SYMBOL),
      keytype = "SYMBOL",
      columns = c(
        "SYMBOL",
        "ENTREZID"
      )
    )
  )

  entrez_rank_df <- dplyr::inner_join(
    symbol_rank_df,
    map_df,
    by = "SYMBOL"
  ) %>%
    dplyr::filter(
      !is.na(ENTREZID),
      ENTREZID != ""
    ) %>%
    dplyr::group_by(
      ENTREZID
    ) %>%
    dplyr::summarise(
      RankScore = max(
        RankScore,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  entrez_rank <- entrez_rank_df$RankScore

  names(entrez_rank) <- entrez_rank_df$ENTREZID

  entrez_rank <- sort(
    entrez_rank,
    decreasing = TRUE
  )

  gsego_fun <- clusterProfiler::gseGO

  gsego_formals <- names(
    formals(gsego_fun)
  )

  args <- list(
    geneList = entrez_rank,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = GSEA_MIN_GS_SIZE_GO,
    maxGSSize = GSEA_MAX_GS_SIZE,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )

  if ("eps" %in% gsego_formals) {
    args$eps <- 0
  }

  if ("seed" %in% gsego_formals) {
    args$seed <- TRUE
  }

  x <- suppressMessages(
    do.call(
      gsego_fun,
      args
    )
  )

  as.data.frame(x)
}


make_perturbation_plot <- function(df, title_text, subtitle_text) {

  label_df <- df %>%
    dplyr::filter(
      p.adj < PADJ_CUTOFF
    ) %>%
    dplyr::arrange(
      dplyr::desc(distance)
    ) %>%
    dplyr::slice_head(
      n = TOP_LABEL_N
    )

  if (nrow(label_df) < TOP_LABEL_N) {
    additional <- df %>%
      dplyr::filter(
        !(gene %in% label_df$gene)
      ) %>%
      dplyr::arrange(
        dplyr::desc(distance)
      ) %>%
      dplyr::slice_head(
        n = TOP_LABEL_N -
          nrow(label_df)
      )

    label_df <- dplyr::bind_rows(
      label_df,
      additional
    )
  }

  ycap <- min(
    120,
    max(
      10,
      as.numeric(
        stats::quantile(
          df$negLog10Padj,
          probs = 0.995,
          na.rm = TRUE
        )
      )
    )
  )

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Z,
      y = pmin(
        negLog10Padj,
        ycap
      )
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = Significance
      ),
      alpha = 0.75,
      size = 1.6
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "FDR < 0.05" = "#D73027",
        "Not significant" = "#C7CED3"
      )
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(
        PADJ_CUTOFF
      ),
      linetype = "dashed",
      color = "#6E7A83",
      linewidth = 0.45
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#A6B0B7",
      linewidth = 0.4
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(
        label = gene
      ),
      size = 3.0,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.20,
      segment.color = "#7B858D",
      min.segment.length = 0
    ) +
    ggplot2::labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Perturbation Z-score",
      y = expression(
        -log[10](
          "adjusted P"
        )
      ),
      color = NULL
    ) +
    theme_pub(
      9
    ) +
    ggplot2::theme(
      legend.position = "top"
    )
}


make_top_gene_plot <- function(df, title_text) {

  plot_df <- df %>%
    dplyr::arrange(
      dplyr::desc(distance)
    ) %>%
    dplyr::slice_head(
      n = TOP_GENE_N
    )

  plot_df$gene <- factor(
    plot_df$gene,
    levels = rev(
      plot_df$gene
    )
  )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = Z,
      y = gene,
      color = Significance
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = min(
          Z,
          0
        ),
        xend = Z,
        yend = gene
      ),
      linewidth = 0.7,
      color = "#C3CBD1"
    ) +
    ggplot2::geom_point(
      size = 3
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#87939B",
      linewidth = 0.45
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "FDR < 0.05" = "#D73027",
        "Not significant" = "#78909C"
      )
    ) +
    ggplot2::labs(
      title = title_text,
      subtitle = "Ranked by virtual-KO network perturbation",
      x = "Perturbation Z-score",
      y = NULL,
      color = NULL
    ) +
    theme_pub(
      9
    ) +
    ggplot2::theme(
      legend.position = "top"
    )
}


# =============================================================================
# 3. Load Figure6 master object
# =============================================================================

if (!file.exists(MASTER_RDS)) {
  stop(
    "Cannot find ",
    MASTER_RDS,
    " in the current working directory."
  )
}

master <- readRDS(
  MASTER_RDS
)

if (is.null(master$core$endo_clean)) {
  stop(
    "master$core$endo_clean was not found."
  )
}

if (is.null(master$core$key_state)) {
  stop(
    "master$core$key_state was not found."
  )
}

endo_clean <- master$core$endo_clean

KEY_STATE <- as.character(
  master$core$key_state
)

if (!inherits(endo_clean, "Seurat")) {
  stop(
    "master$core$endo_clean is not a Seurat object."
  )
}

if (length(KEY_STATE) != 1) {
  stop(
    "master$core$key_state should contain exactly one label."
  )
}

CAPILLARY_STATES[3] <- KEY_STATE

DefaultAssay(
  endo_clean
) <- "RNA"

meta <- endo_clean@meta.data

cat(
  "\n============================================================\n",
  "SERPINE1 optimized virtual KO\n",
  "============================================================\n",
  "Data source: ", MASTER_RDS, "\n",
  "endo_clean cells: ", ncol(endo_clean), "\n",
  "KEY_STATE: ", KEY_STATE, "\n",
  sep = ""
)


# =============================================================================
# 4. Detect metadata columns by actual values
# =============================================================================

annotation_col <- find_meta_col_by_values(
  meta,
  values = CAPILLARY_STATES,
  preferred = c(
    "Annotation",
    "annotation",
    "EndothelialSubtype",
    "endothelial_subtype",
    "EndoSubtype",
    "endo_subtype",
    "Fig5_State",
    "final_state",
    "subtype",
    "Subtype",
    "endo_state",
    "Endo_state",
    "celltype",
    "Celltype"
  ),
  label = "capillary-state"
)

condition_col <- find_meta_col_by_values(
  meta,
  values = CONDITION_KEEP,
  preferred = c(
    "Condition",
    "condition",
    "Group",
    "group"
  ),
  label = "condition"
)

sample_col <- pick_meta_col(
  meta,
  c(
    "Sample",
    "sample",
    "orig.ident",
    "Orig.ident"
  ),
  "sample"
)

cat(
  "\nMetadata columns:\n",
  "  annotation_col = ", annotation_col, "\n",
  "  condition_col  = ", condition_col, "\n",
  "  sample_col     = ", sample_col, "\n",
  sep = ""
)


# =============================================================================
# 5. Select ALL eligible cells
# =============================================================================

if (!(TARGET_GENE %in% rownames(endo_clean))) {
  stop(
    TARGET_GENE,
    " was not found in the RNA assay."
  )
}

condition_vec <- normalize_meta_value(
  meta[[condition_col]]
)

state_vec <- normalize_meta_value(
  meta[[annotation_col]]
)

primary_cells <- rownames(meta)[
  condition_vec == CONDITION_KEEP &
    state_vec %in% CAPILLARY_STATES
]

validation_cells <- rownames(meta)[
  condition_vec == CONDITION_KEEP &
    state_vec == KEY_STATE
]

obj_primary <- subset(
  endo_clean,
  cells = primary_cells
)

obj_validation <- subset(
  endo_clean,
  cells = validation_cells
)

cat(
  "\nALL eligible input cells retained:\n",
  "  Primary HR capillary continuum: ", ncol(obj_primary), "\n",
  "  Validation HR ", KEY_STATE, ": ", ncol(obj_validation), "\n",
  sep = ""
)

if (ncol(obj_primary) < NC_NCELLS_TARGET) {
  warning(
    "Primary context has fewer than ",
    NC_NCELLS_TARGET,
    " cells."
  )
}

if (ncol(obj_validation) < NC_NCELLS_TARGET) {
  warning(
    "Validation context has fewer than ",
    NC_NCELLS_TARGET,
    " cells."
  )
}


# =============================================================================
# 6. Sample/state QC
# =============================================================================

qc_primary_state <- obj_primary@meta.data %>%
  tibble::rownames_to_column(
    "Cell"
  ) %>%
  dplyr::transmute(
    Sample = normalize_meta_value(
      .data[[sample_col]]
    ),
    State = normalize_meta_value(
      .data[[annotation_col]]
    )
  ) %>%
  dplyr::count(
    Sample,
    State,
    name = "Cells"
  )

qc_validation_state <- obj_validation@meta.data %>%
  tibble::rownames_to_column(
    "Cell"
  ) %>%
  dplyr::transmute(
    Sample = normalize_meta_value(
      .data[[sample_col]]
    ),
    State = normalize_meta_value(
      .data[[annotation_col]]
    )
  ) %>%
  dplyr::count(
    Sample,
    State,
    name = "Cells"
  )

write.csv(
  qc_primary_state,
  file.path(
    DIR_QC,
    "PRIMARY_allEligible_sample_state_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  qc_validation_state,
  file.path(
    DIR_QC,
    "VALIDATION_allEligible_sample_counts.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 7. Extract raw counts
# =============================================================================

tmp_primary <- get_raw_counts(
  obj_primary
)

obj_primary <- tmp_primary$object
counts_primary <- tmp_primary$counts

tmp_validation <- get_raw_counts(
  obj_validation
)

obj_validation <- tmp_validation$object
counts_validation <- tmp_validation$counts

rm(
  tmp_primary,
  tmp_validation
)


# =============================================================================
# 8. SERPINE1 and epithelial-marker QC
# =============================================================================

EPITHELIAL_MARKERS <- c(
  "SFTPC",
  "SFTPA1",
  "SFTPA2",
  "EPCAM",
  "SCGB1A1",
  "SCGB3A1"
)

ENDO_MARKERS <- c(
  "PECAM1",
  "VWF",
  "CLDN5",
  "CDH5",
  "EMCN",
  "KDR"
)

qc_marker_genes <- unique(
  c(
    TARGET_GENE,
    EPITHELIAL_MARKERS,
    ENDO_MARKERS
  )
)

marker_qc_one <- function(counts, obj, context_name) {

  genes_present <- intersect(
    qc_marker_genes,
    rownames(counts)
  )

  meta_i <- obj@meta.data[
    colnames(counts),
    ,
    drop = FALSE
  ]

  rows <- lapply(
    genes_present,
    function(g) {

      x <- as.numeric(
        counts[g, , drop = TRUE]
      )

      data.frame(
        Context = context_name,
        Gene = g,
        PositiveCells = sum(
          x > 0
        ),
        PositiveFraction = mean(
          x > 0
        ),
        MeanRawCount = mean(
          x
        ),
        MedianRawCount = stats::median(
          x
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  dplyr::bind_rows(
    rows
  )
}

marker_qc <- dplyr::bind_rows(
  marker_qc_one(
    counts_primary,
    obj_primary,
    "Primary: HR capillary continuum"
  ),
  marker_qc_one(
    counts_validation,
    obj_validation,
    paste0(
      "Validation: HR ",
      KEY_STATE
    )
  )
)

write.csv(
  marker_qc,
  file.path(
    DIR_QC,
    "SERPINE1_and_epithelial_marker_QC.csv"
  ),
  row.names = FALSE
)

cat(
  "\nMarker QC:\n"
)

print(
  marker_qc
)


# =============================================================================
# 9. Mechanistic programs
#
# SERPINE1 is excluded from the program sets to avoid trivial self-enrichment.
# =============================================================================

program_sets <- list(

  `p53-p21 / cell-cycle arrest` = c(
    "CDKN1A",
    "CDKN2B",
    "GADD45A",
    "GADD45B",
    "BTG2",
    "DDB2",
    "MDM2",
    "SESN1",
    "SESN2",
    "GDF15",
    "ATM",
    "CHEK1",
    "CHEK2"
  ),

  `Inflammatory / SASP` = c(
    "ICAM1",
    "VCAM1",
    "SELE",
    "SELP",
    "IL6",
    "CXCL8",
    "CXCL1",
    "CXCL2",
    "CXCL3",
    "CCL2",
    "TNFAIP3",
    "NFKBIA",
    "IRF1",
    "IL1RL1"
  ),

  `Oxidative stress` = c(
    "HMOX1",
    "NQO1",
    "SOD2",
    "GPX3",
    "TXNRD1",
    "PRDX1",
    "PRDX2",
    "DUSP1",
    "JUN",
    "FOS",
    "ATF3",
    "GCLC",
    "GCLM",
    "HIF1A"
  ),

  `ECM remodeling` = c(
    "FN1",
    "COL4A1",
    "COL4A2",
    "COL1A1",
    "COL1A2",
    "MMP2",
    "MMP14",
    "TIMP1",
    "TIMP3",
    "THBS1",
    "PLAU",
    "PLAUR",
    "ANGPT2",
    "SPARCL1"
  ),

  `Barrier homeostasis` = c(
    "CLDN5",
    "CDH5",
    "TJP1",
    "OCLN",
    "TEK",
    "KDR",
    "ROBO4",
    "ESAM",
    "PECAM1",
    "PTPRB",
    "EMCN",
    "VWF"
  ),

  `Fibrinolysis / coagulation` = c(
    "PLAT",
    "PLAU",
    "PLAUR",
    "F3",
    "THBD",
    "PROCR",
    "VWF",
    "F2R",
    "TFPI",
    "THBS1"
  )
)

program_sets <- lapply(
  program_sets,
  function(x) {
    unique(
      setdiff(
        x,
        TARGET_GENE
      )
    )
  }
)

anchored_genes <- unique(
  unlist(
    program_sets,
    use.names = FALSE
  )
)


# =============================================================================
# 10. Independent feature universes
# =============================================================================

prep_primary <- prepare_context_matrix(
  obj = obj_primary,
  counts = counts_primary,
  target_gene = TARGET_GENE,
  anchored_genes = anchored_genes,
  n_hvgs = N_HVGS_PRIMARY,
  min_detection_fraction = MIN_DETECTION_FRACTION,
  min_detection_cells = MIN_DETECTION_CELLS,
  context_name = "Primary: HR capillary continuum"
)

mat_primary <- prep_primary$matrix

prep_validation <- prepare_context_matrix(
  obj = obj_validation,
  counts = counts_validation,
  target_gene = TARGET_GENE,
  anchored_genes = anchored_genes,
  n_hvgs = N_HVGS_VALIDATION,
  min_detection_fraction = MIN_DETECTION_FRACTION,
  min_detection_cells = MIN_DETECTION_CELLS,
  context_name = paste0(
    "Validation: HR ",
    KEY_STATE
  )
)

mat_validation <- prep_validation$matrix

feature_qc <- dplyr::bind_rows(
  prep_primary$feature_qc,
  prep_validation$feature_qc
)

write.csv(
  feature_qc,
  file.path(
    DIR_QC,
    "independent_feature_universe_QC.csv"
  ),
  row.names = FALSE
)

writeLines(
  rownames(mat_primary),
  con = file.path(
    DIR_PRIMARY,
    "PRIMARY_feature_universe.txt"
  )
)

writeLines(
  rownames(mat_validation),
  con = file.path(
    DIR_VALIDATION,
    "VALIDATION_feature_universe.txt"
  )
)

cat(
  "\nIndependent feature universes:\n"
)

print(
  feature_qc
)


# =============================================================================
# 11. Effective network cell settings and theoretical pool coverage
# =============================================================================

nc_nCells_primary <- min(
  NC_NCELLS_TARGET,
  ncol(mat_primary)
)

nc_nCells_validation <- min(
  NC_NCELLS_TARGET,
  ncol(mat_validation)
)

coverage_primary <- 1 -
  (
    1 -
      nc_nCells_primary /
        ncol(mat_primary)
  ) ^ NC_NNET

coverage_validation <- 1 -
  (
    1 -
      nc_nCells_validation /
        ncol(mat_validation)
  ) ^ NC_NNET

network_qc <- data.frame(
  Context = c(
    "Primary: HR capillary continuum",
    paste0(
      "Validation: HR ",
      KEY_STATE
    )
  ),
  InputPoolCells = c(
    ncol(mat_primary),
    ncol(mat_validation)
  ),
  CellsPerNetwork = c(
    nc_nCells_primary,
    nc_nCells_validation
  ),
  NumberOfNetworks = c(
    NC_NNET,
    NC_NNET
  ),
  ApproxProbabilityCellAppearsAtLeastOnce = c(
    coverage_primary,
    coverage_validation
  ),
  stringsAsFactors = FALSE
)

write.csv(
  network_qc,
  file.path(
    DIR_QC,
    "network_cell_coverage_QC.csv"
  ),
  row.names = FALSE
)

cat(
  "\nNetwork ensemble QC:\n"
)

print(
  network_qc
)


# =============================================================================
# 12. Save sparse input matrices
# =============================================================================

saveRDS(
  mat_primary,
  file.path(
    DIR_OBJECT,
    "PRIMARY_rawCounts_independentFeatures_sparse.rds"
  ),
  compress = FALSE
)

saveRDS(
  mat_validation,
  file.path(
    DIR_OBJECT,
    "VALIDATION_rawCounts_independentFeatures_sparse.rds"
  ),
  compress = FALSE
)


# =============================================================================
# 13. PRIMARY scTenifoldKnk
# =============================================================================

cat(
  "\n============================================================\n",
  "PRIMARY SERPINE1 virtual KO\n",
  "HR capillary endothelial continuum\n",
  "============================================================\n"
)

input_primary <- prepare_knk_input(
  mat_primary
)

result_primary <- run_scTenifoldKnk_compat(
  input_counts = input_primary,
  target_gene = TARGET_GENE,
  nc_nNet = NC_NNET,
  nc_nCells = nc_nCells_primary,
  td_K = TD_K,
  qc_mtThreshold = QC_MT_THRESHOLD,
  qc_minLSize = QC_MIN_LSIZE,
  seed = 123
)

saveRDS(
  result_primary,
  file.path(
    DIR_OBJECT,
    "PRIMARY_SERPINE1_scTenifoldKnk_result.rds"
  ),
  compress = FALSE
)

df_primary <- process_diff_result(
  result_primary,
  TARGET_GENE,
  PADJ_CUTOFF
)

sig_primary <- df_primary[
  df_primary$p.adj < PADJ_CUTOFF,
  ,
  drop = FALSE
]

write.csv(
  df_primary,
  file.path(
    DIR_PRIMARY,
    "PRIMARY_diffRegulation_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  sig_primary,
  file.path(
    DIR_PRIMARY,
    "PRIMARY_diffRegulation_FDR005.csv"
  ),
  row.names = FALSE
)

rank_primary <- make_rank_vector(
  df_primary
)

targeted_gsea_primary <- run_targeted_gsea(
  rank_primary,
  program_sets
)

go_gsea_primary <- run_go_gsea(
  rank_primary
)

write.csv(
  targeted_gsea_primary,
  file.path(
    DIR_PRIMARY,
    "PRIMARY_targeted_program_GSEA.csv"
  ),
  row.names = FALSE
)

write.csv(
  go_gsea_primary,
  file.path(
    DIR_PRIMARY,
    "PRIMARY_GO_BP_GSEA.csv"
  ),
  row.names = FALSE
)

# Free the dense input before running validation.
rm(
  input_primary
)

gc()


# =============================================================================
# 14. VALIDATION scTenifoldKnk
# =============================================================================

cat(
  "\n============================================================\n",
  "VALIDATION SERPINE1 virtual KO\n",
  "HR ", KEY_STATE, "\n",
  "============================================================\n",
  sep = ""
)

input_validation <- prepare_knk_input(
  mat_validation
)

result_validation <- run_scTenifoldKnk_compat(
  input_counts = input_validation,
  target_gene = TARGET_GENE,
  nc_nNet = NC_NNET,
  nc_nCells = nc_nCells_validation,
  td_K = TD_K,
  qc_mtThreshold = QC_MT_THRESHOLD,
  qc_minLSize = QC_MIN_LSIZE,
  seed = 123
)

saveRDS(
  result_validation,
  file.path(
    DIR_OBJECT,
    "VALIDATION_SERPINE1_scTenifoldKnk_result.rds"
  ),
  compress = FALSE
)

df_validation <- process_diff_result(
  result_validation,
  TARGET_GENE,
  PADJ_CUTOFF
)

sig_validation <- df_validation[
  df_validation$p.adj < PADJ_CUTOFF,
  ,
  drop = FALSE
]

write.csv(
  df_validation,
  file.path(
    DIR_VALIDATION,
    "VALIDATION_diffRegulation_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  sig_validation,
  file.path(
    DIR_VALIDATION,
    "VALIDATION_diffRegulation_FDR005.csv"
  ),
  row.names = FALSE
)

rank_validation <- make_rank_vector(
  df_validation
)

targeted_gsea_validation <- run_targeted_gsea(
  rank_validation,
  program_sets
)

go_gsea_validation <- run_go_gsea(
  rank_validation
)

write.csv(
  targeted_gsea_validation,
  file.path(
    DIR_VALIDATION,
    "VALIDATION_targeted_program_GSEA.csv"
  ),
  row.names = FALSE
)

write.csv(
  go_gsea_validation,
  file.path(
    DIR_VALIDATION,
    "VALIDATION_GO_BP_GSEA.csv"
  ),
  row.names = FALSE
)

rm(
  input_validation
)

gc()


# =============================================================================
# 15. Compare results AFTER independent network construction
# =============================================================================

compare_df <- dplyr::inner_join(
  df_primary %>%
    dplyr::select(
      gene,
      distance_primary = distance,
      Z_primary = Z,
      FC_primary = FC,
      padj_primary = p.adj
    ),
  df_validation %>%
    dplyr::select(
      gene,
      distance_validation = distance,
      Z_validation = Z,
      FC_validation = FC,
      padj_validation = p.adj
    ),
  by = "gene"
) %>%
  dplyr::mutate(
    SigPrimary = padj_primary < PADJ_CUTOFF,
    SigValidation = padj_validation < PADJ_CUTOFF,
    ConcordanceClass = dplyr::case_when(
      SigPrimary & SigValidation ~ "Significant in both",
      SigPrimary ~ "Primary only",
      SigValidation ~ "Validation only",
      TRUE ~ "Not significant"
    ),
    log10Distance_primary = log10(
      pmax(
        distance_primary,
        1e-12
      )
    ),
    log10Distance_validation = log10(
      pmax(
        distance_validation,
        1e-12
      )
    )
  )

rho_distance <- suppressWarnings(
  stats::cor.test(
    compare_df$distance_primary,
    compare_df$distance_validation,
    method = "spearman",
    exact = FALSE
  )
)

shared_sig_genes <- intersect(
  sig_primary$gene,
  sig_validation$gene
)

union_sig_genes <- union(
  sig_primary$gene,
  sig_validation$gene
)

jaccard_sig <- if (
  length(union_sig_genes) > 0
) {
  length(shared_sig_genes) /
    length(union_sig_genes)
} else {
  NA_real_
}

universe_compare <- unique(
  compare_df$gene
)

primary_sig_common <- intersect(
  sig_primary$gene,
  universe_compare
)

validation_sig_common <- intersect(
  sig_validation$gene,
  universe_compare
)

overlap_p <- if (
  length(shared_sig_genes) > 0
) {
  stats::phyper(
    q = length(shared_sig_genes) - 1,
    m = length(primary_sig_common),
    n = length(universe_compare) -
      length(primary_sig_common),
    k = length(validation_sig_common),
    lower.tail = FALSE
  )
} else {
  1
}

comparison_summary <- data.frame(
  Metric = c(
    "Primary feature universe",
    "Validation feature universe",
    "Common genes available for comparison",
    "Primary FDR-significant perturbed genes",
    "Validation FDR-significant perturbed genes",
    "Shared FDR-significant genes",
    "Significant-gene Jaccard",
    "Spearman rho of perturbation distance",
    "Spearman P of perturbation distance",
    "Hypergeometric P for significant-gene overlap"
  ),
  Value = c(
    nrow(mat_primary),
    nrow(mat_validation),
    nrow(compare_df),
    nrow(sig_primary),
    nrow(sig_validation),
    length(shared_sig_genes),
    jaccard_sig,
    unname(
      rho_distance$estimate
    ),
    rho_distance$p.value,
    overlap_p
  ),
  stringsAsFactors = FALSE
)

write.csv(
  compare_df,
  file.path(
    DIR_COMPARE,
    "PRIMARY_vs_VALIDATION_gene_level.csv"
  ),
  row.names = FALSE
)

write.csv(
  comparison_summary,
  file.path(
    DIR_COMPARE,
    "PRIMARY_vs_VALIDATION_summary.csv"
  ),
  row.names = FALSE
)

writeLines(
  shared_sig_genes,
  con = file.path(
    DIR_COMPARE,
    "shared_FDR005_genes.txt"
  )
)


# =============================================================================
# 16. Combine targeted GSEA tables
# =============================================================================

gsea_add_context <- function(df, context) {

  if (nrow(df) == 0) {
    return(
      data.frame()
    )
  }

  df$Context <- context

  df
}

targeted_gsea_compare <- dplyr::bind_rows(
  gsea_add_context(
    targeted_gsea_primary,
    "Primary: HR capillary continuum"
  ),
  gsea_add_context(
    targeted_gsea_validation,
    paste0(
      "Validation: HR ",
      KEY_STATE
    )
  )
)

write.csv(
  targeted_gsea_compare,
  file.path(
    DIR_COMPARE,
    "targeted_program_GSEA_PRIMARY_vs_VALIDATION.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 17. Optional extra seed stability
# =============================================================================

if (RUN_EXTRA_SEED_STABILITY) {

  stability_dir <- file.path(
    DIR_COMPARE,
    "seed_stability"
  )

  dir.create(
    stability_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  run_seed_context <- function(mat, context_name, n_cells, seeds) {

    rows <- list()

    for (seed_i in seeds) {

      cat(
        "\nSeed stability: ",
        context_name,
        " | seed ",
        seed_i,
        "\n",
        sep = ""
      )

      input_i <- prepare_knk_input(
        mat
      )

      result_i <- run_scTenifoldKnk_compat(
        input_counts = input_i,
        target_gene = TARGET_GENE,
        nc_nNet = NC_NNET,
        nc_nCells = n_cells,
        td_K = TD_K,
        qc_mtThreshold = QC_MT_THRESHOLD,
        qc_minLSize = QC_MIN_LSIZE,
        seed = seed_i
      )

      df_i <- process_diff_result(
        result_i,
        TARGET_GENE,
        PADJ_CUTOFF
      )

      rows[[as.character(seed_i)]] <- df_i %>%
        dplyr::transmute(
          Context = context_name,
          Seed = seed_i,
          gene,
          distance,
          Z,
          FC,
          p.adj,
          Significant = p.adj <
            PADJ_CUTOFF
        )

      rm(
        input_i,
        result_i,
        df_i
      )

      gc()
    }

    dplyr::bind_rows(
      rows
    )
  }

  stability_primary <- run_seed_context(
    mat_primary,
    "Primary",
    nc_nCells_primary,
    EXTRA_SEEDS
  )

  stability_validation <- run_seed_context(
    mat_validation,
    "Validation",
    nc_nCells_validation,
    EXTRA_SEEDS
  )

  stability_all <- dplyr::bind_rows(
    stability_primary,
    stability_validation
  )

  stability_summary <- stability_all %>%
    dplyr::group_by(
      Context,
      gene
    ) %>%
    dplyr::summarise(
      Runs = dplyr::n(),
      SignificantRuns = sum(
        Significant
      ),
      SignificantFraction = mean(
        Significant
      ),
      MedianZ = stats::median(
        Z,
        na.rm = TRUE
      ),
      MedianDistance = stats::median(
        distance,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  write.csv(
    stability_all,
    file.path(
      stability_dir,
      "seed_stability_all_runs.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    stability_summary,
    file.path(
      stability_dir,
      "seed_stability_summary.csv"
    ),
    row.names = FALSE
  )
}


# =============================================================================
# 18. Plot A | input design and network coverage
# =============================================================================

input_plot_df <- network_qc %>%
  dplyr::mutate(
    Context = factor(
      Context,
      levels = rev(
        Context
      )
    ),
    CoverageLabel = scales::percent(
      ApproxProbabilityCellAppearsAtLeastOnce,
      accuracy = 0.1
    ),
    Label = paste0(
      InputPoolCells,
      " cells; ",
      CellsPerNetwork,
      "/network x ",
      NumberOfNetworks,
      "; approx. pool coverage ",
      CoverageLabel
    )
  )

pA <- ggplot2::ggplot(
  input_plot_df,
  ggplot2::aes(
    x = InputPoolCells,
    y = Context
  )
) +
  ggplot2::geom_col(
    width = 0.62,
    fill = "#607D8B"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = Label
    ),
    hjust = -0.04,
    size = 3.2,
    fontface = "bold"
  ) +
  ggplot2::expand_limits(
    x = max(
      input_plot_df$InputPoolCells
    ) *
      1.60
  ) +
  ggplot2::labs(
    title = "Optimized virtual-knockout design",
    subtitle = "All eligible cells retained in each input pool; larger 20-network ensembles",
    x = "Eligible cells in input pool",
    y = NULL
  ) +
  theme_pub(
    9
  ) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank()
  )


# =============================================================================
# 19. Plot B/C | primary and validation perturbation landscapes
# =============================================================================

pB <- make_perturbation_plot(
  df_primary,
  title_text = "Primary: HR capillary continuum",
  subtitle_text = paste0(
    nrow(sig_primary),
    " genes at FDR < 0.05; full profile used for ranked GSEA"
  )
)

pC <- make_perturbation_plot(
  df_validation,
  title_text = paste0(
    "Validation: HR ",
    KEY_STATE
  ),
  subtitle_text = paste0(
    nrow(sig_validation),
    " genes at FDR < 0.05; full profile used for ranked GSEA"
  )
)


# =============================================================================
# 20. Plot D | top primary perturbed genes
# =============================================================================

pD <- make_top_gene_plot(
  df_primary,
  "Top SERPINE1-dependent network perturbations"
)


# =============================================================================
# 21. Plot E | independent-GRN concordance
# =============================================================================

label_compare <- compare_df %>%
  dplyr::filter(
    ConcordanceClass == "Significant in both"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      distance_primary +
        distance_validation
    )
  ) %>%
  dplyr::slice_head(
    n = 12
  )

if (nrow(label_compare) < 8) {

  extra_labels <- compare_df %>%
    dplyr::filter(
      !(gene %in% label_compare$gene)
    ) %>%
    dplyr::arrange(
      dplyr::desc(
        distance_primary +
          distance_validation
      )
    ) %>%
    dplyr::slice_head(
      n = 8 -
        nrow(label_compare)
    )

  label_compare <- dplyr::bind_rows(
    label_compare,
    extra_labels
  )
}

pE <- ggplot2::ggplot(
  compare_df,
  ggplot2::aes(
    x = log10Distance_primary,
    y = log10Distance_validation,
    color = ConcordanceClass
  )
) +
  ggplot2::geom_point(
    size = 1.7,
    alpha = 0.75
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Significant in both" = "#B2182B",
      "Primary only" = "#E07B72",
      "Validation only" = "#4C78A8",
      "Not significant" = "#C8CED3"
    )
  ) +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "#4A5560",
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  ggrepel::geom_text_repel(
    data = label_compare,
    ggplot2::aes(
      label = gene
    ),
    size = 3.0,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.20,
    segment.color = "#7B858D"
  ) +
  ggplot2::labs(
    title = "State-specific validation of SERPINE1 network dependency",
    subtitle = paste0(
      "Independent GRNs; Spearman rho = ",
      sprintf(
        "%.2f",
        unname(
          rho_distance$estimate
        )
      ),
      "; shared FDR-significant genes = ",
      length(
        shared_sig_genes
      )
    ),
    x = "Primary log10 perturbation distance",
    y = "Validation log10 perturbation distance",
    color = NULL
  ) +
  theme_pub(
    9
  ) +
  ggplot2::theme(
    legend.position = "top"
  )


# =============================================================================
# 22. Plot F | targeted ranked GSEA
# =============================================================================

if (nrow(targeted_gsea_compare) > 0) {

  targeted_plot_df <- targeted_gsea_compare %>%
    dplyr::mutate(
      ID = factor(
        ID,
        levels = rev(
          names(
            program_sets
          )
        )
      ),
      Context = factor(
        Context,
        levels = c(
          "Primary: HR capillary continuum",
          paste0(
            "Validation: HR ",
            KEY_STATE
          )
        )
      ),
      FDRLabel = ifelse(
        p.adjust < PADJ_CUTOFF,
        "FDR < 0.05",
        "FDR >= 0.05"
      )
    )

  pF <- ggplot2::ggplot(
    targeted_plot_df,
    ggplot2::aes(
      x = NES,
      y = ID
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#8D989F",
      linewidth = 0.45
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = NES,
        yend = ID
      ),
      color = "#BDC5CB",
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = -log10(
          pmax(
            p.adjust,
            1e-300
          )
        ),
        color = FDRLabel
      )
    ) +
    ggplot2::facet_wrap(
      ~ Context,
      ncol = 1
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "FDR < 0.05" = "#B2182B",
        "FDR >= 0.05" = "#7893AD"
      )
    ) +
    ggplot2::scale_size_continuous(
      range = c(
        2.5,
        6.2
      ),
      name = "-log10(FDR)"
    ) +
    ggplot2::labs(
      title = "Ranked pathway perturbation after SERPINE1 virtual KO",
      subtitle = "Positive NES indicates concentration toward the most perturbed genes",
      x = "Normalized enrichment score (NES)",
      y = NULL,
      color = NULL
    ) +
    theme_pub(
      8.5
    ) +
    ggplot2::theme(
      legend.position = "right"
    )

} else {

  pF <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "No targeted GSEA results available",
      size = 4
    ) +
    ggplot2::xlim(
      0,
      1
    ) +
    ggplot2::ylim(
      0,
      1
    ) +
    ggplot2::theme_void()
}


# =============================================================================
# 23. Plot G | Primary GO BP ranked GSEA
# =============================================================================

if (nrow(go_gsea_primary) > 0) {

  go_plot_df <- go_gsea_primary %>%
    dplyr::filter(
      is.finite(NES),
      NES > 0
    ) %>%
    dplyr::arrange(
      p.adjust,
      dplyr::desc(
        NES
      )
    ) %>%
    dplyr::slice_head(
      n = TOP_GO_N
    ) %>%
    dplyr::mutate(
      Description = factor(
        Description,
        levels = rev(
          Description
        )
      ),
      FDRLabel = ifelse(
        p.adjust < PADJ_CUTOFF,
        "FDR < 0.05",
        "FDR >= 0.05"
      )
    )

  if (nrow(go_plot_df) > 0) {

    pG <- ggplot2::ggplot(
      go_plot_df,
      ggplot2::aes(
        x = NES,
        y = Description,
        fill = FDRLabel
      )
    ) +
      ggplot2::geom_col(
        width = 0.72
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "FDR < 0.05" = "#B2182B",
          "FDR >= 0.05" = "#8998C0"
        )
      ) +
      ggplot2::labs(
        title = "Primary capillary functional perturbation architecture",
        subtitle = "Top positively enriched GO biological processes from the full perturbation ranking",
        x = "Normalized enrichment score (NES)",
        y = NULL,
        fill = NULL
      ) +
      theme_pub(
        8
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.y = ggplot2::element_text(
          size = 6.8
        )
      )

  } else {

    pG <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = "No positive GO BP GSEA terms available",
        size = 4
      ) +
      ggplot2::xlim(
        0,
        1
      ) +
      ggplot2::ylim(
        0,
        1
      ) +
      ggplot2::theme_void()
  }

} else {

  pG <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "No GO BP GSEA results available",
      size = 4
    ) +
    ggplot2::xlim(
      0,
      1
    ) +
    ggplot2::ylim(
      0,
      1
    ) +
    ggplot2::theme_void()
}


# =============================================================================
# 24. Plot H | epithelial-marker QC
# =============================================================================

marker_plot_df <- marker_qc %>%
  dplyr::filter(
    Gene %in% c(
      TARGET_GENE,
      EPITHELIAL_MARKERS,
      ENDO_MARKERS
    )
  ) %>%
  dplyr::mutate(
    GeneClass = dplyr::case_when(
      Gene == TARGET_GENE ~ "SERPINE1",
      Gene %in% EPITHELIAL_MARKERS ~ "Epithelial marker",
      Gene %in% ENDO_MARKERS ~ "Endothelial marker",
      TRUE ~ "Other"
    )
  )

pH <- ggplot2::ggplot(
  marker_plot_df,
  ggplot2::aes(
    x = PositiveFraction,
    y = Gene,
    color = GeneClass
  )
) +
  ggplot2::geom_point(
    size = 3.0
  ) +
  ggplot2::facet_wrap(
    ~ Context,
    ncol = 1
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "SERPINE1" = "#B2182B",
      "Epithelial marker" = "#D98C3B",
      "Endothelial marker" = "#3C78A8"
    )
  ) +
  ggplot2::labs(
    title = "Lineage-marker QC of virtual-KO input cells",
    subtitle = "Epithelial-marker detection is reported, not automatically filtered",
    x = "Raw-count detection fraction",
    y = NULL,
    color = NULL
  ) +
  theme_pub(
    8.5
  ) +
  ggplot2::theme(
    legend.position = "bottom"
  )


# =============================================================================
# 25. Save plots
# =============================================================================

save_plot(
  pA,
  "A_optimized_design",
  9.0,
  4.0,
  DIR_FINAL
)

save_plot(
  pB,
  "B_PRIMARY_perturbation",
  6.2,
  5.1,
  DIR_FINAL
)

save_plot(
  pC,
  "C_VALIDATION_perturbation",
  6.2,
  5.1,
  DIR_FINAL
)

save_plot(
  pD,
  "D_PRIMARY_top_perturbed",
  5.8,
  5.1,
  DIR_FINAL
)

save_plot(
  pE,
  "E_independentGRN_concordance",
  6.2,
  5.1,
  DIR_FINAL
)

save_plot(
  pF,
  "F_targeted_ranked_GSEA",
  7.2,
  6.0,
  DIR_FINAL
)

save_plot(
  pG,
  "G_PRIMARY_GO_BP_ranked_GSEA",
  7.2,
  5.6,
  DIR_FINAL
)

save_plot(
  pH,
  "H_input_lineage_marker_QC",
  7.0,
  5.5,
  DIR_QC
)


# =============================================================================
# 26. Final figure A-G
# =============================================================================

row1 <- patchwork::wrap_plots(
  pA,
  ncol = 1
)

row2 <- patchwork::wrap_plots(
  pB,
  pC,
  ncol = 2,
  widths = c(
    1,
    1
  )
)

row3 <- patchwork::wrap_plots(
  pD,
  pE,
  ncol = 2,
  widths = c(
    0.95,
    1.05
  )
)

row4 <- patchwork::wrap_plots(
  pF,
  pG,
  ncol = 2,
  widths = c(
    1.05,
    0.95
  )
)

Figure_SERPINE1_KO <- patchwork::wrap_plots(
  row1,
  row2,
  row3,
  row4,
  ncol = 1,
  heights = c(
    0.60,
    1.0,
    1.0,
    1.08
  )
) +
  patchwork::plot_annotation(
    title = expression(
      "Virtual knockout maps " *
        italic(
          SERPINE1
        ) *
        "-dependent regulatory perturbation across the reperfusion-injured capillary endothelial continuum"
    ),
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 14,
        color = "#202A33",
        margin = ggplot2::margin(
          b = 8
        )
      )
    )
  )

pdf_device <- if (
  capabilities("cairo")
) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

ggplot2::ggsave(
  filename = file.path(
    DIR_FINAL,
    "SERPINE1_virtualKO_OPTIMIZED_FINAL.pdf"
  ),
  plot = Figure_SERPINE1_KO,
  width = 13.5,
  height = 18.2,
  units = "in",
  device = pdf_device,
  bg = "white",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    DIR_FINAL,
    "SERPINE1_virtualKO_OPTIMIZED_FINAL.png"
  ),
  plot = Figure_SERPINE1_KO,
  width = 13.5,
  height = 18.2,
  units = "in",
  dpi = 500,
  device = ragg::agg_png,
  bg = "white",
  limitsize = FALSE
)


# =============================================================================
# 27. Save lightweight reusable bundle
#
# Large scTenifoldKnk result objects are stored separately and referenced by path
# to keep this bundle manageable.
# =============================================================================

bundle <- list(

  design = list(
    master_rds = MASTER_RDS,
    target_gene = TARGET_GENE,
    condition = CONDITION_KEEP,
    primary_states = CAPILLARY_STATES,
    validation_state = KEY_STATE,
    nc_nNet = NC_NNET,
    nc_nCells_primary = nc_nCells_primary,
    nc_nCells_validation = nc_nCells_validation,
    primary_pool_cells = ncol(mat_primary),
    validation_pool_cells = ncol(mat_validation)
  ),

  metadata = list(
    annotation_col = annotation_col,
    condition_col = condition_col,
    sample_col = sample_col
  ),

  qc = list(
    primary_sample_state = qc_primary_state,
    validation_sample_state = qc_validation_state,
    feature_qc = feature_qc,
    network_qc = network_qc,
    marker_qc = marker_qc
  ),

  primary = list(
    result_path = file.path(
      DIR_OBJECT,
      "PRIMARY_SERPINE1_scTenifoldKnk_result.rds"
    ),
    diff_all = df_primary,
    diff_sig = sig_primary,
    targeted_GSEA = targeted_gsea_primary,
    GO_BP_GSEA = go_gsea_primary
  ),

  validation = list(
    result_path = file.path(
      DIR_OBJECT,
      "VALIDATION_SERPINE1_scTenifoldKnk_result.rds"
    ),
    diff_all = df_validation,
    diff_sig = sig_validation,
    targeted_GSEA = targeted_gsea_validation,
    GO_BP_GSEA = go_gsea_validation
  ),

  comparison = list(
    gene_level = compare_df,
    summary = comparison_summary,
    shared_sig_genes = shared_sig_genes,
    targeted_GSEA = targeted_gsea_compare
  ),

  panels = list(
    A = pA,
    B = pB,
    C = pC,
    D = pD,
    E = pE,
    F = pF,
    G = pG,
    H_QC = pH
  ),

  assembled = Figure_SERPINE1_KO
)

saveRDS(
  bundle,
  file.path(
    DIR_OBJECT,
    "SERPINE1_virtualKO_OPTIMIZED_ANALYSIS_BUNDLE.rds"
  ),
  compress = FALSE
)


# =============================================================================
# 28. Summary
# =============================================================================

summary_lines <- c(
  paste0(
    "Data source: ",
    MASTER_RDS
  ),
  paste0(
    "Target gene: ",
    TARGET_GENE
  ),
  paste0(
    "Condition: ",
    CONDITION_KEEP
  ),
  paste0(
    "Primary states: ",
    paste(
      CAPILLARY_STATES,
      collapse = " | "
    )
  ),
  paste0(
    "Validation state: ",
    KEY_STATE
  ),
  paste0(
    "Primary input pool cells: ",
    ncol(
      mat_primary
    )
  ),
  paste0(
    "Validation input pool cells: ",
    ncol(
      mat_validation
    )
  ),
  paste0(
    "Primary features: ",
    nrow(
      mat_primary
    )
  ),
  paste0(
    "Validation features: ",
    nrow(
      mat_validation
    )
  ),
  paste0(
    "Cells per network: primary=",
    nc_nCells_primary,
    ", validation=",
    nc_nCells_validation
  ),
  paste0(
    "Number of networks: ",
    NC_NNET
  ),
  paste0(
    "Approx primary pool coverage: ",
    signif(
      coverage_primary,
      5
    )
  ),
  paste0(
    "Approx validation pool coverage: ",
    signif(
      coverage_validation,
      5
    )
  ),
  paste0(
    "Primary FDR<0.05 genes: ",
    nrow(
      sig_primary
    )
  ),
  paste0(
    "Validation FDR<0.05 genes: ",
    nrow(
      sig_validation
    )
  ),
  paste0(
    "Shared FDR<0.05 genes: ",
    length(
      shared_sig_genes
    )
  ),
  paste0(
    "Distance Spearman rho: ",
    signif(
      unname(
        rho_distance$estimate
      ),
      4
    )
  ),
  paste0(
    "Distance Spearman P: ",
    signif(
      rho_distance$p.value,
      4
    )
  ),
  paste0(
    "Significant overlap P: ",
    signif(
      overlap_p,
      4
    )
  ),
  paste0(
    "scTenifoldKnk version: ",
    as.character(
      utils::packageVersion(
        "scTenifoldKnk"
      )
    )
  )
)

writeLines(
  summary_lines,
  con = file.path(
    OUTDIR,
    "SERPINE1_virtualKO_OPTIMIZED_summary.txt"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    OUTDIR,
    "sessionInfo_SERPINE1_virtualKO_OPTIMIZED.txt"
  )
)


# =============================================================================
# 29. Console report
# =============================================================================

cat(
  "\n============================================================\n",
  "SERPINE1 optimized virtual KO completed\n",
  "============================================================\n",
  "PRIMARY\n",
  "  Eligible input pool: ", ncol(mat_primary), " cells\n",
  "  Features: ", nrow(mat_primary), "\n",
  "  Cells per network: ", nc_nCells_primary, "\n",
  "  Networks: ", NC_NNET, "\n",
  "  Approx pool coverage: ", scales::percent(coverage_primary, accuracy = 0.1), "\n",
  "  FDR-significant genes: ", nrow(sig_primary), "\n",
  "\nVALIDATION\n",
  "  Eligible input pool: ", ncol(mat_validation), " cells\n",
  "  Features: ", nrow(mat_validation), "\n",
  "  Cells per network: ", nc_nCells_validation, "\n",
  "  Networks: ", NC_NNET, "\n",
  "  Approx pool coverage: ", scales::percent(coverage_validation, accuracy = 0.1), "\n",
  "  FDR-significant genes: ", nrow(sig_validation), "\n",
  "\nCONCORDANCE\n",
  "  Common genes compared: ", nrow(compare_df), "\n",
  "  Shared FDR-significant genes: ", length(shared_sig_genes), "\n",
  "  Perturbation-distance Spearman rho: ", sprintf("%.3f", unname(rho_distance$estimate)), "\n",
  "  P: ", format.pval(rho_distance$p.value, digits = 3), "\n",
  "  Overlap P: ", format.pval(overlap_p, digits = 3), "\n",
  "\nFinal figure:\n",
  file.path(
    DIR_FINAL,
    "SERPINE1_virtualKO_OPTIMIZED_FINAL.pdf"
  ),
  "\n",
  "\nBundle:\n",
  file.path(
    DIR_OBJECT,
    "SERPINE1_virtualKO_OPTIMIZED_ANALYSIS_BUNDLE.rds"
  ),
  "\n",
  sep = ""
)

cat(
  "\nTop PRIMARY perturbed genes:\n"
)

print(
  df_primary %>%
    dplyr::select(
      gene,
      distance,
      Z,
      FC,
      p.adj
    ) %>%
    dplyr::slice_head(
      n = 15
    )
)

cat(
  "\nTop VALIDATION perturbed genes:\n"
)

print(
  df_validation %>%
    dplyr::select(
      gene,
      distance,
      Z,
      FC,
      p.adj
    ) %>%
    dplyr::slice_head(
      n = 15
    )
)

cat(
  "\nTargeted ranked GSEA comparison:\n"
)

if (nrow(targeted_gsea_compare) > 0) {
  print(
    targeted_gsea_compare %>%
      dplyr::select(
        Context,
        ID,
        setSize,
        enrichmentScore,
        NES,
        pvalue,
        p.adjust,
        core_enrichment
      )
  )
} else {
  cat(
    "No targeted GSEA results available.\n"
  )
}

cat(
  "\nIMPORTANT INTERPRETATION:\n",
  "distance/Z/FC are regulatory-network perturbation metrics.\n",
  "They are not transcriptional expression directions.\n",
  sep = ""
)



}

# =============================================================================
# STAGE figure6_perturbation_plots
# =============================================================================
PIPELINE$figure6_perturbation_plots <- function() {
# scripts/08_virtual_perturbation_plots.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

# =============================================================================
# Manuscript Figure 6 | SERPINE1 virtual-KO reviewer correction
# Reuse the saved *_ANALYSIS_BUNDLE.rds; DO NOT rerun scTenifoldKnk.
# A-G: nominal-vs-FDR correction; explain unsigned perturbation-ranked GSEA;
#      donors represented in both nested HR analysis pools.
# Compatible with R 4.3.x, ggplot2 >=3.4, patchwork >=1.2.
# =============================================================================
options(stringsAsFactors = FALSE)
required_pkgs <- c("ggplot2", "patchwork", "dplyr", "tidyr", "ggrepel")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)) stop("Please install packages: ", paste(missing_pkgs, collapse=", "))
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr)
  library(tidyr); library(ggrepel)
})

# --- 1. Paths. Set this to the exact file shown in your 05_OBJECTS directory.
ROOT <- repo_output("Figure6")
BUNDLE_FILE <- repo_input("perturbation_bundle")
if (!file.exists(BUNDLE_FILE)) {
  message("Bundle not found at default path. Choose your *_ANALYSIS_BUNDLE.rds file.")
  BUNDLE_FILE <- file.choose()
}
if (!file.exists(BUNDLE_FILE)) stop("RDS not found: ", BUNDLE_FILE)
OUT <- repo_output("Figure6/Final")
PANEL_DIR <- file.path(OUT, "Panels")
TABLE_DIR <- file.path(OUT, "Tables")
SUP_DIR <- file.path(OUT, "Supplementary")
for (dd in c(OUT,PANEL_DIR,TABLE_DIR,SUP_DIR)) dir.create(dd,recursive=TRUE,showWarnings=FALSE)
FONT <- "Arial"  # If Arial is unavailable on your machine, switch to "sans".
COL_PRIMARY <- "#2166AC"
COL_VALIDATION <- "#D73027"
COL_NONSIG <- "#C8D0D6"
FDR_CUT <- .05
PLOT_DPI <- 600

# --- 2. Bundle checks; no data recomputation.
x <- readRDS(BUNDLE_FILE)
req <- list(x$primary$diff_all,x$validation$diff_all,x$comparison$gene_level,
            x$comparison$targeted_GSEA,x$qc$network_qc)
if (any(vapply(req,is.null,logical(1)))) stop("Bundle missing primary/validation/compare/GSEA/network QC data")
check_fields <- function(df, wanted, label) {
  miss <- setdiff(wanted, names(df))
  if (length(miss)) stop(label, ": missing columns ", paste(miss,collapse=", "))
}
for (nm in c("primary","validation")) {
  check_fields(x[[nm]]$diff_all,c("gene","Z","p.value","p.adj","distance"),nm)
}
check_fields(x$comparison$gene_level,
             c("gene","log10Distance_primary","log10Distance_validation",
               "padj_primary","padj_validation"),"comparison$gene_level")
check_fields(x$comparison$targeted_GSEA,c("ID","Context","NES","p.adjust"),"targeted_GSEA")
check_fields(x$qc$network_qc,c("Context","InputPoolCells","CellsPerNetwork","NumberOfNetworks"),"network_qc")

# --- 3. Utility functions.
write_table <- function(df,filename) utils::write.csv(df,file.path(TABLE_DIR,filename),row.names=FALSE,na="")
save_plot <- function(p, path_base, width, height) {
  pdfdev <- if(capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(paste0(path_base,".pdf"),p,width=width,height=height,
                  units="in",device=pdfdev,bg="white",limitsize=FALSE)
  pngdev <- if(requireNamespace("ragg",quietly=TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(paste0(path_base,".png"),p,width=width,height=height,
                  units="in",dpi=PLOT_DPI,device=pngdev,bg="white",limitsize=FALSE)
}
fmt <- function(v) ifelse(is.na(v),"NA",formatC(v,format="fg",digits=3))
thm <- function(base_size=9) ggplot2::theme_bw(base_size=base_size,base_family=FONT) +
 theme(text=element_text(color="black",face="bold",family=FONT),
       axis.text=element_text(color="black",face="bold"),
       axis.title=element_text(color="black",face="bold"),
       strip.background=element_rect(fill="white",color="black",linewidth=.55),
       strip.text=element_text(color="black",face="bold"),
       legend.title=element_text(color="black",face="bold"),
       legend.text=element_text(color="black"),
       panel.border=element_rect(color="black",fill=NA,linewidth=.55),
       panel.grid.minor=element_blank(),
       panel.grid.major=element_line(color="#EFF2F4",linewidth=.25),
       plot.title=element_blank(),plot.subtitle=element_blank(),
       plot.tag=element_text(color="black",face="bold",size=12),
       plot.margin=margin(5,7,5,7))

prep <- function(df, context) {
  df <- as.data.frame(df)
  df$gene <- as.character(df$gene)
  df$Z <- as.numeric(df$Z)
  df$p.value <- as.numeric(df$p.value)
  df$p.adj <- as.numeric(df$p.adj)
  df$distance <- as.numeric(df$distance)
  df$Qlabel <- ifelse(is.finite(df$p.adj)&df$p.adj<FDR_CUT,
                     "BH-FDR < 0.05","Not FDR-significant")
  df$Yq <- -base::log10(pmax(df$p.adj,1e-300))
  df$Context <- context
  df
}
a <- prep(x$primary$diff_all,"Primary capillary continuum")
b <- prep(x$validation$diff_all,"Stress-responsive EC subset")
write_table(a,"01_Primary_ALL_genes_with_p_and_FDR.csv")
write_table(b,"02_Subset_ALL_genes_with_p_and_FDR.csv")
counts <- dplyr::bind_rows(
  data.frame(Context="Primary capillary continuum", Tested=sum(is.finite(a$p.adj)),
             NominalP_lt_0.05=sum(a$p.value<.05,na.rm=TRUE),
             BH_FDR_lt_0.05=sum(a$p.adj<.05,na.rm=TRUE)),
  data.frame(Context="Stress-responsive EC subset", Tested=sum(is.finite(b$p.adj)),
             NominalP_lt_0.05=sum(b$p.value<.05,na.rm=TRUE),
             BH_FDR_lt_0.05=sum(b$p.adj<.05,na.rm=TRUE)))
write_table(counts,"00_Nominal_vs_FDR_counts.csv")
print(counts)

# --- 4. Panel A: sampling setup, not an independent-donor validation.
nq <- as.data.frame(x$qc$network_qc)
nq$Group <- ifelse(grepl("^Primary",nq$Context),"Primary capillary continuum","Stress-responsive EC subset")
nq$Group <- factor(nq$Group, levels=c("Stress-responsive EC subset","Primary capillary continuum"))
if ("ApproxProbabilityCellAppearsAtLeastOnce" %in% colnames(nq)) {
  coverage_text <- paste0(" | sampling inclusion ",sprintf("%.1f",100*nq$ApproxProbabilityCellAppearsAtLeastOnce),"%")
} else coverage_text <- ""
nq$BarText <- paste0(format(nq$InputPoolCells,big.mark=",")," cells | ",
                     nq$CellsPerNetwork,"/network x ",nq$NumberOfNetworks,coverage_text)
pA <- ggplot(nq,aes(InputPoolCells,Group))+
 geom_col(fill="#E97B70",width=.58)+
 geom_text(aes(label=BarText),hjust=-.04,size=3.2,fontface="bold",color="black")+
 scale_x_continuous(limits=c(0,max(nq$InputPoolCells)*1.68),expand=expansion(mult=c(0,0)))+
 labs(x="Eligible HR cells in input pool",y=NULL)+thm(9)+
 theme(panel.grid.major.y=element_blank(),axis.text.y=element_text(size=9))

# --- 5. B/C: BH-FDR is the ONLY significance classification.
max_y <- max(c(a$Yq,b$Yq)[is.finite(c(a$Yq,b$Yq))],na.rm=TRUE)
ycap <- max(2, min(10, ceiling(max_y)))
xlim <- max(3, ceiling(max(abs(c(a$Z,b$Z)),na.rm=TRUE)*10)/10)
plot_land <- function(df, nlab=10) {
  ldf <- df |>
    dplyr::filter(is.finite(Z),is.finite(p.adj),p.adj<FDR_CUT) |>
    dplyr::arrange(p.adj,dplyr::desc(Z)) |>
    dplyr::slice_head(n=nlab)
  ggplot(df,aes(x=Z,y=pmin(Yq,ycap)))+
    geom_point(aes(color=Qlabel),size=1.5,alpha=.85,na.rm=TRUE)+
    geom_hline(yintercept=-log10(.05),linetype="dashed",linewidth=.5,color="#687989")+
    geom_vline(xintercept=0,linetype="dashed",linewidth=.4,color="#99A5AF")+
    ggrepel::geom_text_repel(data=ldf,aes(label=gene),size=3,
      family=FONT,fontface="bold",max.overlaps=Inf,box.padding=.25,point.padding=.15,
      color="black",min.segment.length=0)+
    scale_color_manual(values=c("BH-FDR < 0.05"=COL_VALIDATION,
                              "Not FDR-significant"=COL_NONSIG),drop=FALSE)+
    coord_cartesian(xlim=c(-xlim,xlim),ylim=c(0,ycap),clip="off")+
    labs(x="Network perturbation Z-score",y=expression(-log[10]("BH-FDR")),color=NULL)+
    thm(9)+theme(legend.position="top",legend.text=element_text(size=7.2))
}
pB <- plot_land(a,12)
pC <- plot_land(b,10)

# --- 6. D: 15 most affected primary genes; FDR-only highlighting.
top <- a |>
 dplyr::filter(is.finite(Z)) |>
 dplyr::arrange(dplyr::desc(Z)) |>
 dplyr::slice_head(n=15)
top$gene <- factor(top$gene,levels=rev(top$gene))
pD <- ggplot(top,aes(x=Z,y=gene))+
 geom_segment(aes(x=0,xend=Z,yend=gene),color="#BEC8D0",linewidth=.7)+
 geom_point(aes(color=Qlabel),size=2.7)+
 scale_color_manual(values=c("BH-FDR < 0.05"=COL_VALIDATION,
                             "Not FDR-significant"=COL_NONSIG),drop=FALSE)+
 labs(x="Network perturbation Z-score",y=NULL,color=NULL)+
 thm(9)+theme(legend.position="top",panel.grid.major.y=element_blank(),axis.text.y=element_text(size=8))
write_table(as.data.frame(top),"03_Panel_D_Top_15.csv")

# --- 7. E: nested subset concordance, NOT external validation.
e <- as.data.frame(x$comparison$gene_level)
e$PrimaryFDR <- as.numeric(e$padj_primary)
e$SubsetFDR <- as.numeric(e$padj_validation)
e$Class <- ifelse(e$PrimaryFDR<FDR_CUT & e$SubsetFDR<FDR_CUT,"Both",
             ifelse(e$PrimaryFDR<FDR_CUT,"Primary only",
             ifelse(e$SubsetFDR<FDR_CUT,"Subset only","Neither")))
e$Class[is.na(e$Class)] <- "Not comparable"
colors_e <- c("Both"="#A9222D","Primary only"="#E18A79",
              "Subset only"="#4875A8","Neither"="#CED6DB","Not comparable"="#ECEFF1")
keep_e <- is.finite(e$log10Distance_primary)&is.finite(e$log10Distance_validation)
rho_e <- if(sum(keep_e)>2) suppressWarnings(cor(e$log10Distance_primary[keep_e],
        e$log10Distance_validation[keep_e],method="spearman")) else NA_real_
predeclared <- c("ICAM1","IRF1","IL1RL1","PTPRB","EDNRB","VWF","SOD2","HES1")
e_labels <- e[e$gene %in% predeclared & keep_e,,drop=FALSE]
pE <- ggplot(e,aes(log10Distance_primary,log10Distance_validation))+
 geom_point(aes(color=Class),size=1.5,alpha=.8,na.rm=TRUE)+
 geom_smooth(method="lm",formula=y~x,se=FALSE,linetype="dashed",color="#75818B",linewidth=.65)+
 ggrepel::geom_text_repel(data=e_labels,aes(label=gene),family=FONT,fontface="bold",
        color="black",size=2.7,max.overlaps=Inf)+
 annotate("label",x=-Inf,y=Inf,hjust=-.05,vjust=1.1,
          label=paste0("Descriptive Spearman rho = ",fmt(rho_e)),
          size=2.9,fontface="bold",family=FONT,color="black")+
 scale_color_manual(values=colors_e)+
 labs(x="Primary log10 perturbation distance",y="Subset log10 perturbation distance",color=NULL)+
 thm(9)+theme(legend.position="top",legend.text=element_text(size=7))
write_table(e,"04_Panel_E_Gene_Overlap_FDR.csv")

# --- 8. F: NES is *relative position in unsigned perturbation ranking*.
g <- as.data.frame(x$comparison$targeted_GSEA)
g$NES <- as.numeric(g$NES)
g$p.adjust <- as.numeric(g$p.adjust)
g$Panel <- ifelse(grepl("^Primary",g$Context),"Primary","Subset")
g$Panel <- factor(g$Panel,levels=c("Primary","Subset"))
g$Qlabel <- ifelse(is.finite(g$p.adjust)&g$p.adjust<FDR_CUT,
                   "BH-FDR < 0.05","Not FDR-significant")
g$FDRscore <- pmin(-log10(pmax(g$p.adjust,1e-300)),8)
g$ID <- factor(g$ID,levels=rev(unique(g$ID)))
pF <- ggplot(g,aes(NES,ID))+
 geom_vline(xintercept=0,linetype="dashed",linewidth=.4,color="#98A5AE")+
 geom_segment(aes(x=0,xend=NES,yend=ID),color="#BEC9D0",linewidth=.65)+
 geom_point(aes(size=FDRscore,color=Qlabel),alpha=.95)+
 facet_wrap(~Panel,nrow=1)+
 scale_color_manual(values=c("BH-FDR < 0.05"=COL_VALIDATION,
                             "Not FDR-significant"=COL_PRIMARY),drop=FALSE)+
 scale_size_continuous(range=c(2.2,5.6),name="-log10(BH-FDR)")+
 labs(x="NES (higher vs lower network-perturbation rank)",y=NULL,color=NULL)+
 thm(8.5)+theme(legend.position="bottom",strip.text=element_text(size=8),
                axis.text.y=element_text(size=7))
write_table(g,"05_Panel_F_Precomputed_Program_GSEA.csv")

# --- 9. G: keep preselected genes, but FDR-significance only.
map <- data.frame(
 Module=c(rep("Inflammatory / adhesion",4),rep("Endothelial / barrier",4),
          rep("Vascular remodeling",4),rep("Stress response",4)),
 gene=c("ICAM1","IRF1","IL1RL1","CCL2","PTPRB","VWF","EMCN","KLF2",
        "EDNRB","HPGD","FCN3","THBD","SOD2","HES1","IER3","MT2A"),
 stringsAsFactors=FALSE)
m1 <- a[,c("gene","Z","p.value","p.adj")];m1$Panel <- "Primary"
m2 <- b[,c("gene","Z","p.value","p.adj")];m2$Panel <- "Subset"
mg <- dplyr::bind_rows(m1,m2)
gg <- tidyr::crossing(gene=map$gene,Panel=c("Primary","Subset")) |>
 dplyr::left_join(map,by="gene") |>
 dplyr::left_join(mg,by=c("gene","Panel")) |>
 dplyr::mutate(Tested=is.finite(Z)&is.finite(p.adj),
        Qlabel=ifelse(!Tested,"Not tested",ifelse(p.adj<FDR_CUT,"BH-FDR < 0.05","Not FDR-significant")),
        FDRscore=ifelse(Tested,pmin(-log10(pmax(p.adj,1e-300)),8),NA_real_))
gg$Panel <- factor(gg$Panel,levels=c("Primary","Subset"))
gg$Module <- factor(gg$Module,levels=c("Inflammatory / adhesion","Endothelial / barrier",
                                  "Vascular remodeling","Stress response"))
gg$gene <- factor(gg$gene,levels=rev(map$gene))
zt <- max(3,min(5,ceiling(max(abs(gg$Z),na.rm=TRUE)*10)/10))
pG <- ggplot() +
 geom_point(data=gg[!gg$Tested,,drop=FALSE],aes(Panel,gene),shape=4,size=2,
            color="#B9C2C8")+
 geom_point(data=gg[gg$Tested,,drop=FALSE],aes(Panel,gene,size=FDRscore,fill=Z,color=Qlabel),
            shape=21,stroke=.85)+
 facet_grid(Module~.,scales="free_y",space="free_y",switch="y")+
 scale_fill_gradient2(low="#427DAF",mid="#F5F6F7",high="#C73934",midpoint=0,
                      limits=c(-zt,zt),oob=scales::squish,name="Perturbation Z")+
 scale_color_manual(values=c("BH-FDR < 0.05"="#151515","Not FDR-significant"="#B5BDC5"),
                    name=NULL)+
 scale_size_continuous(range=c(1.6,5.3),name="-log10(BH-FDR)")+
 labs(x=NULL,y=NULL)+thm(8.5)+
 theme(legend.position="bottom",strip.placement="outside",axis.text.y=element_text(size=7.5),
       legend.text=element_text(size=7),strip.text.y.left=element_text(angle=0,size=7),
       panel.grid.major.y=element_blank())
write_table(gg,"06_Panel_G_Predeclared_Genes_FDR.csv")

# --- 10. Supplementary: actual donors represented; no pseudo-replicated inference.
combine_qc <- function(df, tag) {
 if(is.null(df)) return(NULL)
 if(!all(c("Sample","Cells") %in% names(df))) return(NULL)
 dd <- as.data.frame(df)
 dd$Context <- tag
 dd
}
qc_samples <- dplyr::bind_rows(
 combine_qc(x$qc$primary_sample_state,"Primary capillary continuum"),
 combine_qc(x$qc$validation_sample_state,"Stress-responsive EC subset"))
if(nrow(qc_samples)>0) {
  if(!"State" %in% names(qc_samples)) qc_samples$State <- "All"
  write_table(qc_samples,"07_Donor_Contribution_by_Context_and_State.csv")
  qcp <- qc_samples |>
    dplyr::group_by(Context,Sample) |>
    dplyr::summarise(Cells=sum(Cells),.groups="drop")
  qcp$Sample <- factor(qcp$Sample,levels=paste0("HR",1:6))
  pS <- ggplot(qcp,aes(Sample,Cells,fill=Context))+
    geom_col(width=.72,show.legend=FALSE)+
    facet_wrap(~Context,scales="free_y",ncol=1)+
    scale_fill_manual(values=c("Primary capillary continuum"=COL_PRIMARY,
                               "Stress-responsive EC subset"=COL_VALIDATION))+
    labs(x="Biological donor (HR samples only)",y="Eligible cells")+thm(9)
  save_plot(pS,file.path(SUP_DIR,"S1_Donor_contribution_to_virtualKO"),8.4,6)
  write_table(qcp,"08_Donor_Summary_and_Context_Counts.csv")
}
if(!is.null(x$qc$marker_qc)) write_table(as.data.frame(x$qc$marker_qc),"09_Epithelial_Marker_QC.csv")
if(!is.null(x$qc$feature_qc)) write_table(as.data.frame(x$qc$feature_qc),"10_Feature_Universe_QC.csv")
if(!is.null(x$design)) {
  utils::capture.output(str(x$design),file=file.path(TABLE_DIR,"11_Analysis_Design.txt"))
}

# --- 11. Export A-G and final layout.
panels <- list(A=pA,B=pB,C=pC,D=pD,E=pE,F=pF,G=pG)
sizes <- list(A=c(13.6,2.7),B=c(6.6,5),C=c(6.6,5),D=c(6.4,5),
              E=c(6.8,5),F=c(6.5,5.9),G=c(6.9,5.9))
for(nm in names(panels)) save_plot(panels[[nm]],
    file.path(PANEL_DIR,paste0("Figure6",nm,"_REVISED")),sizes[[nm]][1],sizes[[nm]][2])
allfig <- (pA / (pB+pC) / (pD+pE) / (pF+pG))+
  patchwork::plot_layout(heights=c(.52,1.02,1.02,1.20)) +
  patchwork::plot_annotation(tag_levels="A",theme=theme(
    plot.tag=element_text(color="black",face="bold",family=FONT,size=13)))
save_plot(allfig,file.path(OUT,"Figure6_REVISED_A_to_G_FDR_CORRECTED"),14.6,19.3)
saveRDS(list(panels=panels,figure=allfig,input_bundle=BUNDLE_FILE,
             FDR_summary=counts,donor_sample_counts=qc_samples),
        file.path(OUT,"Figure6_REVISED_PLOTTING_OBJECTS.rds"))

# --- 12. Honest audit note: do not claim completed biological sensitivity testing.
notes <- c(
 "Manuscript Figure 6 (local Figure7): virtual SERPINE1 KO",
 paste("Input analysis bundle:",BUNDLE_FILE),
 "No scTenifoldKnk / network reconstruction was run by this script.",
 "B/C/D/G significance is based on saved gene-level BH-FDR (p.adj < 0.05).",
 "F GSEA NES orders genes by their NETWORK PERTURBATION Z, not signed expression change.",
 "Positive/negative NES must NOT be called pathway activation/inhibition after KO.",
 "Primary and stress-EC subset share cells: subset sensitivity, NOT independent validation.",
 "One 20-network ensemble does NOT establish random-seed robustness.",
 "Pooled HR networks do NOT constitute six independent biological networks.",
 "To claim donor robustness, rerun scTenifoldKnk under predeclared leave-one-donor-out/seed settings.",
 "These data do not establish causal ICAM1 dependence or true expression response.",
 "If SFTPC was excluded from analysis, verify the feature universe and marker QC, not only plotted labels."
)
writeLines(notes,con=file.path(OUT,"READ_ME_FIRST_METHODOLOGY_NOTES.txt"),useBytes=TRUE)
cat("\nFigure6 review-corrected A-G saved in:\n",OUT,"\n",sep="")



}

# =============================================================================
# STAGE figure7_metabolic
# =============================================================================
PIPELINE$figure7_metabolic <- function() {
# scripts/09_metabolic_trajectories.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.


# =============================================================================
# FIG8 Stage 1 | Single-cell metabolic discovery
# No Monocle2 / DDRTree / VGAM / igraph required
#
# Inputs:
#   Figure6_MASTER_DATA_ONLY.rds
#   Figure6_PREMATURE_SENESCENCE_DATA_ONLY.rds
#
# Logic:
#   A. Unbiased screen of broad KEGG metabolic programs along Fig6 pseudotime
#   B. Heatmap of top dynamic metabolic programs
#   C. Relative transition ordering vs SERPINE1 / inflammatory-SASP / ICAM1
#
# Statistical unit:
#   Cells are first aggregated within Sample x pseudotime-bin.
#   GAM includes Condition plus sample random intercept.
#
# Interpretation:
#   "Earlier" means earlier along inferred pseudotime, not real-time causality.
# =============================================================================
options(stringsAsFactors = FALSE)
set.seed(1234)

required_pkgs <- c(
  "Seurat", "SeuratObject", "Matrix",
  "dplyr", "tidyr", "tibble",
  "mgcv", "msigdbr",
  "ggplot2", "patchwork", "scales", "ragg"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(mgcv)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ragg)
})

MASTER_RDS <- repo_input("network_master")
FIG6_DATA_RDS <- repo_input("trajectory_data")

OUTDIR <- repo_output("Figure7/Metabolic")
PANEL_DIR <- file.path(OUTDIR, "01_panels")
TABLE_DIR <- file.path(OUTDIR, "02_tables")
OBJECT_DIR <- file.path(OUTDIR, "03_objects")
FINAL_DIR <- file.path(OUTDIR, "04_final")

for (d in c(PANEL_DIR, TABLE_DIR, OBJECT_DIR, FINAL_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(MASTER_RDS)) stop("Cannot find ", MASTER_RDS)
if (!file.exists(FIG6_DATA_RDS)) stop("Cannot find ", FIG6_DATA_RDS)

master <- readRDS(MASTER_RDS)
fig6 <- readRDS(FIG6_DATA_RDS)

if (is.null(master$core$endo_clean)) stop("master$core$endo_clean is missing.")
if (is.null(fig6$trajectory_metadata)) stop("fig6$trajectory_metadata is missing.")

endo_clean <- master$core$endo_clean
traj <- as.data.frame(fig6$trajectory_metadata, stringsAsFactors = FALSE)

Seurat::DefaultAssay(endo_clean) <- "RNA"
expr <- SeuratObject::GetAssayData(endo_clean, assay = "RNA", slot = "data")

required_traj_cols <- c(
  "Cell", "Sample", "Condition", "EndothelialState",
  "Pseudotime", "PseudotimeNorm", "PT_Bin", "PT_Bin_Mid",
  "SERPINE1", "Inflammatory_SASP"
)

missing_traj_cols <- setdiff(required_traj_cols, colnames(traj))
if (length(missing_traj_cols) > 0) {
  stop("Missing trajectory columns: ", paste(missing_traj_cols, collapse = ", "))
}

CAPILLARY_STATES <- c(
  "Homeostatic capillary EC",
  "Transitional capillary EC",
  "Stress-responsive capillary EC"
)

traj <- traj %>%
  dplyr::filter(
    Cell %in% colnames(expr),
    EndothelialState %in% CAPILLARY_STATES,
    is.finite(PseudotimeNorm),
    is.finite(PT_Bin_Mid)
  ) %>%
  dplyr::mutate(
    Sample = factor(Sample),
    Condition = factor(Condition, levels = c("CIT", "HR")),
    EndothelialState = factor(EndothelialState, levels = CAPILLARY_STATES)
  )

trajectory_cells <- traj$Cell

cat(
  "\nMatched trajectory cells: ", length(trajectory_cells),
  "\nSamples: ", length(unique(traj$Sample)),
  "\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# Broad metabolic universe
# Lipid metabolism is included but not preselected as the answer.
# -----------------------------------------------------------------------------

METABOLIC_UNIVERSE <- tibble::tribble(
  ~Pathway, ~Category, ~Pattern,
  "Glycolysis / gluconeogenesis", "Carbohydrate", "GLYCOLYSIS_GLUCONEOGENESIS$",
  "Citrate cycle (TCA cycle)", "Carbohydrate", "CITRATE_CYCLE_TCA_CYCLE$",
  "Pentose phosphate pathway", "Carbohydrate", "PENTOSE_PHOSPHATE_PATHWAY$",
  "Fructose and mannose metabolism", "Carbohydrate", "FRUCTOSE_AND_MANNOSE_METABOLISM$",
  "Galactose metabolism", "Carbohydrate", "GALACTOSE_METABOLISM$",
  "Pyruvate metabolism", "Carbohydrate", "PYRUVATE_METABOLISM$",
  "Propanoate metabolism", "Carbohydrate", "PROPANOATE_METABOLISM$",
  "Butanoate metabolism", "Carbohydrate", "BUTANOATE_METABOLISM$",
  "Amino sugar and nucleotide sugar metabolism", "Carbohydrate", "AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM$",
  "Oxidative phosphorylation", "Energy", "OXIDATIVE_PHOSPHORYLATION$",
  "Fatty acid degradation", "Lipid", "FATTY_ACID_DEGRADATION$",
  "Fatty acid elongation", "Lipid", "FATTY_ACID_ELONGATION$",
  "Biosynthesis of unsaturated fatty acids", "Lipid", "BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS$",
  "Glycerolipid metabolism", "Lipid", "GLYCEROLIPID_METABOLISM$",
  "Glycerophospholipid metabolism", "Lipid", "GLYCEROPHOSPHOLIPID_METABOLISM$",
  "Ether lipid metabolism", "Lipid", "ETHER_LIPID_METABOLISM$",
  "Sphingolipid metabolism", "Lipid", "SPHINGOLIPID_METABOLISM$",
  "Arachidonic acid metabolism", "Lipid", "ARACHIDONIC_ACID_METABOLISM$",
  "Linoleic acid metabolism", "Lipid", "LINOLEIC_ACID_METABOLISM$",
  "Alpha-linolenic acid metabolism", "Lipid", "ALPHA_LINOLENIC_ACID_METABOLISM$",
  "Steroid biosynthesis", "Lipid", "STEROID_BIOSYNTHESIS$",
  "Steroid hormone biosynthesis", "Lipid", "STEROID_HORMONE_BIOSYNTHESIS$",
  "Primary bile acid biosynthesis", "Lipid", "PRIMARY_BILE_ACID_BIOSYNTHESIS$",
  "Glycosphingolipid biosynthesis - lacto/neolacto", "Lipid", "GLYCOSPHINGOLIPID_BIOSYNTHESIS_LACTO_AND_NEOLACTO_SERIES$",
  "Glycosphingolipid biosynthesis - ganglio", "Lipid", "GLYCOSPHINGOLIPID_BIOSYNTHESIS_GANGLIO_SERIES$",
  "Glycosphingolipid biosynthesis - globo/isoglobo", "Lipid", "GLYCOSPHINGOLIPID_BIOSYNTHESIS_GLOBO_AND_ISOGLOBO_SERIES$",
  "Alanine, aspartate and glutamate metabolism", "Amino acid", "ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM$",
  "Arginine and proline metabolism", "Amino acid", "ARGININE_AND_PROLINE_METABOLISM$",
  "Cysteine and methionine metabolism", "Amino acid", "CYSTEINE_AND_METHIONINE_METABOLISM$",
  "Glycine, serine and threonine metabolism", "Amino acid", "GLYCINE_SERINE_AND_THREONINE_METABOLISM$",
  "Valine, leucine and isoleucine degradation", "Amino acid", "VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION$",
  "Phenylalanine metabolism", "Amino acid", "PHENYLALANINE_METABOLISM$",
  "Tyrosine metabolism", "Amino acid", "TYROSINE_METABOLISM$",
  "Tryptophan metabolism", "Amino acid", "TRYPTOPHAN_METABOLISM$",
  "Histidine metabolism", "Amino acid", "HISTIDINE_METABOLISM$",
  "Lysine degradation", "Amino acid", "LYSINE_DEGRADATION$",
  "Purine metabolism", "Nucleotide", "PURINE_METABOLISM$",
  "Pyrimidine metabolism", "Nucleotide", "PYRIMIDINE_METABOLISM$",
  "Glutathione metabolism", "Cofactor / redox", "GLUTATHIONE_METABOLISM$",
  "Nicotinate and nicotinamide metabolism", "Cofactor / redox", "NICOTINATE_AND_NICOTINAMIDE_METABOLISM$",
  "Riboflavin metabolism", "Cofactor / redox", "RIBOFLAVIN_METABOLISM$",
  "Thiamine metabolism", "Cofactor / redox", "THIAMINE_METABOLISM$",
  "Porphyrin metabolism", "Cofactor / redox", "PORPHYRIN_METABOLISM$",
  "One-carbon pool by folate", "Cofactor / redox", "ONE_CARBON_POOL_BY_FOLATE$",
  "Inositol phosphate metabolism", "Other metabolism", "INOSITOL_PHOSPHATE_METABOLISM$",
  "N-glycan biosynthesis", "Glycan", "N_GLYCAN_BIOSYNTHESIS$",
  "Other glycan degradation", "Glycan", "OTHER_GLYCAN_DEGRADATION$"
)

METABOLIC_UNIVERSE$PathwayID <- sprintf("M%03d", seq_len(nrow(METABOLIC_UNIVERSE)))

# -----------------------------------------------------------------------------
# msigdbr 7.5.1 compatible retrieval
# -----------------------------------------------------------------------------

msig <- tryCatch(
  msigdbr::msigdbr(
    species = "Homo sapiens",
    category = "C2",
    subcategory = "CP:KEGG"
  ),
  error = function(e) {
    message("CP:KEGG call failed; falling back to all C2 gene sets.")
    msigdbr::msigdbr(species = "Homo sapiens", category = "C2")
  }
)

if (nrow(msig) == 0) stop("msigdbr returned no gene sets.")

msig$gs_name <- as.character(msig$gs_name)
msig$gene_symbol <- as.character(msig$gene_symbol)

MIN_GENE_DETECTION <- 0.01
MIN_MATCHED_GENES <- 6

gene_sets <- list()
gene_set_qc <- list()

for (i in seq_len(nrow(METABOLIC_UNIVERSE))) {
  pathway_id <- METABOLIC_UNIVERSE$PathwayID[i]
  matched_sets <- unique(
    msig$gs_name[
      grepl(METABOLIC_UNIVERSE$Pattern[i], msig$gs_name, ignore.case = TRUE)
    ]
  )

  genes_i <- unique(msig$gene_symbol[msig$gs_name %in% matched_sets])
  genes_i <- intersect(genes_i, rownames(expr))

  if (length(genes_i) > 0) {
    det <- Matrix::rowMeans(expr[genes_i, trajectory_cells, drop = FALSE] > 0)
    genes_i <- names(det)[det >= MIN_GENE_DETECTION]
  }

  gene_sets[[pathway_id]] <- genes_i

  gene_set_qc[[pathway_id]] <- data.frame(
    PathwayID = pathway_id,
    Pathway = METABOLIC_UNIVERSE$Pathway[i],
    Category = METABOLIC_UNIVERSE$Category[i],
    NGenes = length(genes_i),
    MatchedMSigDBSets = paste(matched_sets, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

gene_set_qc <- dplyr::bind_rows(gene_set_qc)
eligible_ids <- gene_set_qc %>%
  dplyr::filter(NGenes >= MIN_MATCHED_GENES) %>%
  dplyr::pull(PathwayID)

write.csv(
  gene_set_qc,
  file.path(TABLE_DIR, "Metabolic_gene_set_QC.csv"),
  row.names = FALSE
)

if (length(eligible_ids) < 10) {
  stop(
    "Only ", length(eligible_ids),
    " pathways had >= ", MIN_MATCHED_GENES,
    " expressed genes. Inspect Metabolic_gene_set_QC.csv."
  )
}

# -----------------------------------------------------------------------------
# Transparent pathway scoring: mean log-normalized expression, then z-score
# -----------------------------------------------------------------------------

zsafe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

score_mat <- matrix(
  NA_real_,
  nrow = length(trajectory_cells),
  ncol = length(eligible_ids),
  dimnames = list(trajectory_cells, eligible_ids)
)

for (pathway_id in eligible_ids) {
  genes_i <- gene_sets[[pathway_id]]
  raw_score <- Matrix::colMeans(expr[genes_i, trajectory_cells, drop = FALSE])
  score_mat[, pathway_id] <- zsafe(raw_score)
}

# -----------------------------------------------------------------------------
# Sample x pseudotime-bin aggregation
# -----------------------------------------------------------------------------

MIN_CELLS_PER_SAMPLE_BIN <- 5

score_df <- as.data.frame(score_mat, check.names = FALSE)
score_df$Cell <- rownames(score_df)

sample_bin_scores <- score_df %>%
  dplyr::left_join(
    traj %>% dplyr::select(Cell, Sample, Condition, PT_Bin, PT_Bin_Mid),
    by = "Cell"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(eligible_ids),
    names_to = "PathwayID",
    values_to = "Score"
  ) %>%
  dplyr::group_by(PathwayID, Sample, Condition, PT_Bin, PT_Bin_Mid) %>%
  dplyr::summarise(
    Cells = dplyr::n(),
    Score = mean(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    Cells >= MIN_CELLS_PER_SAMPLE_BIN,
    is.finite(Score),
    is.finite(PT_Bin_Mid)
  ) %>%
  dplyr::mutate(Sample = factor(Sample))

# -----------------------------------------------------------------------------
# GAM helpers
# -----------------------------------------------------------------------------

PT_GRID <- seq(0, 1, length.out = 101)

transition_midpoint <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 5) return(NA_real_)

  d <- abs(diff(y))
  total <- sum(d, na.rm = TRUE)
  if (!is.finite(total) || total <= 1e-8) return(NA_real_)

  cumulative <- c(0, cumsum(d) / total)
  idx <- which(cumulative >= 0.5)[1]
  if (length(idx) == 0 || is.na(idx)) return(NA_real_)
  x[idx]
}

fit_one_gam <- function(data_i) {
  data_i <- data_i %>%
    dplyr::filter(
      is.finite(Score),
      is.finite(PT_Bin_Mid),
      !is.na(Sample),
      !is.na(Condition)
    )

  n_pt <- length(unique(data_i$PT_Bin_Mid))
  n_samples <- length(unique(data_i$Sample))

  if (nrow(data_i) < 20 || n_pt < 5 || n_samples < 4) return(NULL)

  k_use <- min(5, n_pt - 1)

  fit <- tryCatch(
    mgcv::gam(
      Score ~ Condition +
        s(PT_Bin_Mid, k = k_use, bs = "tp") +
        s(Sample, bs = "re"),
      data = data_i,
      method = "REML"
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  sfit <- summary(fit)
  smooth_row <- grep("^s\\(PT_Bin_Mid\\)", rownames(sfit$s.table))

  p_value <- if (length(smooth_row) > 0) {
    as.numeric(sfit$s.table[smooth_row[1], "p-value"])
  } else {
    NA_real_
  }

  ref_condition <- if ("CIT" %in% levels(data_i$Condition)) {
    "CIT"
  } else {
    levels(data_i$Condition)[1]
  }

  ref_sample <- levels(data_i$Sample)[1]

  nd <- data.frame(
    Condition = factor(rep(ref_condition, length(PT_GRID)), levels = levels(data_i$Condition)),
    PT_Bin_Mid = PT_GRID,
    Sample = factor(rep(ref_sample, length(PT_GRID)), levels = levels(data_i$Sample))
  )

  pred <- tryCatch(
    as.numeric(
      stats::predict(
        fit,
        newdata = nd,
        type = "response",
        exclude = "s(Sample)"
      )
    ),
    error = function(e) rep(NA_real_, length(PT_GRID))
  )

  list(
    fit = fit,
    PValue = p_value,
    DevianceExplained = sfit$dev.expl,
    DynamicRange = diff(range(pred, na.rm = TRUE)),
    LateMinusEarly =
      mean(pred[PT_GRID >= 0.80], na.rm = TRUE) -
      mean(pred[PT_GRID <= 0.20], na.rm = TRUE),
    TransitionMidpoint = transition_midpoint(PT_GRID, pred),
    Prediction = data.frame(
      Pseudotime = PT_GRID,
      PredictedScore = pred
    )
  )
}

# -----------------------------------------------------------------------------
# Unbiased metabolic screen
# -----------------------------------------------------------------------------

fits <- list()
stats_list <- list()
pred_list <- list()

for (pathway_id in eligible_ids) {
  message("Fitting ", pathway_id)

  fit_i <- fit_one_gam(
    sample_bin_scores %>% dplyr::filter(PathwayID == pathway_id)
  )

  fits[[pathway_id]] <- fit_i

  if (!is.null(fit_i)) {
    stats_list[[pathway_id]] <- data.frame(
      PathwayID = pathway_id,
      PValue = fit_i$PValue,
      DevianceExplained = fit_i$DevianceExplained,
      DynamicRange = fit_i$DynamicRange,
      LateMinusEarly = fit_i$LateMinusEarly,
      TransitionMidpoint = fit_i$TransitionMidpoint,
      stringsAsFactors = FALSE
    )

    pred_i <- fit_i$Prediction
    pred_i$PathwayID <- pathway_id
    pred_list[[pathway_id]] <- pred_i
  }
}

screen <- dplyr::bind_rows(stats_list) %>%
  dplyr::left_join(
    METABOLIC_UNIVERSE %>%
      dplyr::select(PathwayID, Pathway, Category),
    by = "PathwayID"
  ) %>%
  dplyr::mutate(
    FDR = stats::p.adjust(PValue, method = "BH"),
    NegLog10FDR = -log10(pmax(FDR, 1e-300)),
    DiscoveryScore =
      DynamicRange *
      sqrt(pmax(DevianceExplained, 0)) *
      pmin(NegLog10FDR, 10)
  ) %>%
  dplyr::arrange(FDR, dplyr::desc(DiscoveryScore)) %>%
  dplyr::mutate(Rank = dplyr::row_number())

predictions <- dplyr::bind_rows(pred_list) %>%
  dplyr::left_join(
    screen %>%
      dplyr::select(
        PathwayID,
        Pathway,
        Category,
        FDR,
        DiscoveryScore,
        TransitionMidpoint
      ),
    by = "PathwayID"
  )

write.csv(
  screen,
  file.path(TABLE_DIR, "Unbiased_metabolic_pathway_screen.csv"),
  row.names = FALSE
)

write.csv(
  screen %>% dplyr::filter(Category == "Lipid"),
  file.path(TABLE_DIR, "Lipid_pathways_in_unbiased_ranking.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Panel A | pathway ranking
# -----------------------------------------------------------------------------

CATEGORY_COLORS <- c(
  "Carbohydrate" = "#4C78A8",
  "Energy" = "#72A0C1",
  "Lipid" = "#D95F02",
  "Amino acid" = "#59A14F",
  "Nucleotide" = "#B279A2",
  "Cofactor / redox" = "#E3A321",
  "Glycan" = "#8C8C8C",
  "Other metabolism" = "#7A7A7A"
)

plot_A_df <- screen %>%
  dplyr::filter(is.finite(DiscoveryScore)) %>%
  dplyr::arrange(dplyr::desc(DiscoveryScore)) %>%
  dplyr::slice_head(n = 18) %>%
  dplyr::mutate(
    Pathway = factor(Pathway, levels = rev(Pathway))
  )

pA <- ggplot2::ggplot(
  plot_A_df,
  ggplot2::aes(x = DiscoveryScore, y = Pathway)
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = DiscoveryScore,
      yend = Pathway
    ),
    color = "#CCD4D9",
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    ggplot2::aes(size = NegLog10FDR, fill = Category),
    shape = 21,
    color = "#263640",
    stroke = 0.5
  ) +
  ggplot2::scale_fill_manual(
    values = CATEGORY_COLORS,
    drop = FALSE,
    name = "Metabolic class"
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.8, 7.0),
    name = "-log10(FDR)"
  ) +
  ggplot2::labs(
    x = "Dynamic discovery score",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(face = "bold", size = 8.2),
    legend.title = ggplot2::element_text(face = "bold"),
    panel.grid.major.x = ggplot2::element_line(
      color = "#E8EDF0",
      linewidth = 0.35
    )
  )

# -----------------------------------------------------------------------------
# Panel B | top dynamic metabolic-program heatmap
# -----------------------------------------------------------------------------

top_B_ids <- screen %>%
  dplyr::filter(
    is.finite(DiscoveryScore),
    is.finite(TransitionMidpoint)
  ) %>%
  dplyr::arrange(dplyr::desc(DiscoveryScore)) %>%
  dplyr::slice_head(n = 12) %>%
  dplyr::pull(PathwayID)

plot_B_df <- predictions %>%
  dplyr::filter(PathwayID %in% top_B_ids) %>%
  dplyr::group_by(PathwayID) %>%
  dplyr::mutate(SmoothZ = zsafe(PredictedScore)) %>%
  dplyr::ungroup()

heat_order <- screen %>%
  dplyr::filter(PathwayID %in% top_B_ids) %>%
  dplyr::arrange(TransitionMidpoint) %>%
  dplyr::pull(Pathway)

plot_B_df$Pathway <- factor(
  plot_B_df$Pathway,
  levels = rev(heat_order)
)

pB <- ggplot2::ggplot(
  plot_B_df,
  ggplot2::aes(
    x = Pseudotime,
    y = Pathway,
    fill = SmoothZ
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Relative\nactivity"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    x = "Normalized inferred pseudotime",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(face = "bold", size = 8.0),
    axis.ticks.y = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(face = "bold")
  )

# -----------------------------------------------------------------------------
# Panel C | reference response programs and relative transition ordering
# -----------------------------------------------------------------------------

if (!"ICAM1" %in% rownames(expr)) stop("ICAM1 is not present in RNA matrix.")

ref_cell <- traj %>%
  dplyr::select(
    Cell,
    Sample,
    Condition,
    PT_Bin,
    PT_Bin_Mid,
    SERPINE1,
    Inflammatory_SASP
  )

ref_cell$ICAM1 <- zsafe(
  as.numeric(expr["ICAM1", ref_cell$Cell, drop = TRUE])
)

ref_bin <- ref_cell %>%
  tidyr::pivot_longer(
    cols = c(SERPINE1, Inflammatory_SASP, ICAM1),
    names_to = "Feature",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    Feature = dplyr::recode(
      Feature,
      Inflammatory_SASP = "Inflammatory / SASP"
    )
  ) %>%
  dplyr::group_by(
    Feature,
    Sample,
    Condition,
    PT_Bin,
    PT_Bin_Mid
  ) %>%
  dplyr::summarise(
    Cells = dplyr::n(),
    Score = mean(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(Cells >= MIN_CELLS_PER_SAMPLE_BIN) %>%
  dplyr::mutate(Sample = factor(Sample))

ref_stats <- list()

for (feature_i in unique(ref_bin$Feature)) {
  fit_i <- fit_one_gam(
    ref_bin %>% dplyr::filter(Feature == feature_i)
  )

  if (!is.null(fit_i)) {
    ref_stats[[feature_i]] <- data.frame(
      Feature = feature_i,
      FeatureType = "Endothelial response",
      TransitionMidpoint = fit_i$TransitionMidpoint,
      LateMinusEarly = fit_i$LateMinusEarly,
      stringsAsFactors = FALSE
    )
  }
}

ref_stats <- dplyr::bind_rows(ref_stats)

metabolic_order <- screen %>%
  dplyr::filter(
    is.finite(TransitionMidpoint),
    is.finite(DiscoveryScore)
  ) %>%
  dplyr::arrange(dplyr::desc(DiscoveryScore)) %>%
  dplyr::slice_head(n = 6) %>%
  dplyr::transmute(
    Feature = Pathway,
    FeatureType = "Metabolic program",
    TransitionMidpoint,
    LateMinusEarly
  )

order_df <- dplyr::bind_rows(
  metabolic_order,
  ref_stats
) %>%
  dplyr::filter(is.finite(TransitionMidpoint)) %>%
  dplyr::arrange(TransitionMidpoint)

order_df$Feature <- factor(
  order_df$Feature,
  levels = rev(order_df$Feature)
)

pC <- ggplot2::ggplot(
  order_df,
  ggplot2::aes(
    x = TransitionMidpoint,
    y = Feature
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = TransitionMidpoint,
      yend = Feature,
      color = FeatureType
    ),
    linewidth = 0.85,
    alpha = 0.5
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = FeatureType),
    shape = 21,
    size = 4.0,
    color = "#263640",
    stroke = 0.5
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Metabolic program" = "#2A788E",
      "Endothelial response" = "#B2182B"
    ),
    guide = "none"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Metabolic program" = "#2A788E",
      "Endothelial response" = "#B2182B"
    ),
    name = NULL
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  ggplot2::labs(
    x = "Normalized inferred pseudotime",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(face = "bold", size = 8.0),
    legend.position = "top",
    panel.grid.major.x = ggplot2::element_line(
      color = "#E8EDF0",
      linewidth = 0.35
    )
  )

write.csv(
  order_df,
  file.path(TABLE_DIR, "Relative_program_transition_order.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Assemble and save
# -----------------------------------------------------------------------------

Figure7_ABC <- pA /
  (pB | pC) +
  patchwork::plot_layout(
    heights = c(1.12, 1)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = expression(
      "Single-cell metabolic screening reveals dynamic remodeling along the " *
        italic(SERPINE1) *
        "-associated endothelial injury continuum"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 15
      ),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 13
      )
    )
  )

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

save_one <- function(p, name, width, height) {
  ggplot2::ggsave(
    filename = file.path(PANEL_DIR, paste0(name, ".pdf")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = pdf_device,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = file.path(PANEL_DIR, paste0(name, ".png")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = ragg::agg_png,
    bg = "white"
  )
}

save_one(pA, "Figure7A_Unbiased_metabolic_screen", 9.5, 6.5)
save_one(pB, "Figure7B_Dynamic_metabolic_heatmap", 8.0, 5.8)
save_one(pC, "Figure7C_Program_transition_order", 8.0, 5.8)

ggplot2::ggsave(
  filename = file.path(FINAL_DIR, "Figure7_STAGE1_A_C_PREVIEW.pdf"),
  plot = Figure7_ABC,
  width = 15,
  height = 12.5,
  units = "in",
  device = pdf_device,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(FINAL_DIR, "Figure7_STAGE1_A_C_PREVIEW.png"),
  plot = Figure7_ABC,
  width = 15,
  height = 12.5,
  units = "in",
  dpi = 500,
  device = ragg::agg_png,
  bg = "white"
)

saveRDS(
  list(
    trajectory_metadata = traj,
    pathway_universe = METABOLIC_UNIVERSE,
    gene_set_qc = gene_set_qc,
    screen = screen,
    predictions = predictions,
    transition_order = order_df,
    panels = list(A = pA, B = pB, C = pC),
    assembled = Figure7_ABC
  ),
  file.path(
    OBJECT_DIR,
    "Figure7_STAGE1_SC_METABOLIC_DISCOVERY_BUNDLE.rds"
  ),
  compress = FALSE
)

cat(
  "\n============================================================\n",
  "TOP 15 METABOLIC PROGRAMS\n",
  "============================================================\n"
)

print(
  screen %>%
    dplyr::select(
      Rank,
      Pathway,
      Category,
      FDR,
      DynamicRange,
      LateMinusEarly,
      TransitionMidpoint,
      DiscoveryScore
    ) %>%
    dplyr::slice_head(n = 15),
  row.names = FALSE
)

cat(
  "\n============================================================\n",
  "LIPID PATHWAYS IN THE UNBIASED RANKING\n",
  "============================================================\n"
)

print(
  screen %>%
    dplyr::filter(Category == "Lipid") %>%
    dplyr::select(
      Rank,
      Pathway,
      FDR,
      DynamicRange,
      LateMinusEarly,
      TransitionMidpoint,
      DiscoveryScore
    ),
  row.names = FALSE
)

cat(
  "\n============================================================\n",
  "RELATIVE TRANSITION ORDER\n",
  "============================================================\n"
)

print(
  order_df %>%
    dplyr::arrange(TransitionMidpoint),
  row.names = FALSE
)

cat(
  "\nNo Monocle2 / igraph was loaded.\n",
  "Interpret timing only as ordering along inferred pseudotime.\n",
  sep = ""
)

print(Figure7_ABC)



# =============================================================================
# FIG8 Stage 1 | FINAL LOGIC A-C | NO RERUN
#
# A: unbiased metabolic screen -> lipid pathways emerge naturally
# B: top dynamic metabolic programs along established Fig6 pseudotime
# C: focused ordering of lipid-remodeling programs versus SASP/ICAM1/SERPINE1
#
# No Monocle2, no igraph, no metabolic analysis rerun.
# "Earlier" means earlier along inferred pseudotime only.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ragg)
})

options(stringsAsFactors = FALSE)

BUNDLE_FILE <- file.path(repo_output("Figure7/Metabolic"),"03_objects","Figure7_STAGE1_SC_METABOLIC_DISCOVERY_BUNDLE.rds")

OUTDIR <- repo_output("Figure7/Metabolic")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(BUNDLE_FILE)) {
  stop("Cannot find analysis bundle:\n", BUNDLE_FILE)
}

x <- readRDS(BUNDLE_FILE)
screen <- x$screen
predictions <- x$predictions
transition_order_old <- x$transition_order

required_cols <- c(
  "PathwayID", "Pathway", "Category", "FDR", "DynamicRange",
  "LateMinusEarly", "TransitionMidpoint", "DiscoveryScore"
)

if (!all(required_cols %in% colnames(screen))) {
  stop("Saved screen table is missing required columns.")
}

if (is.null(predictions) || nrow(predictions) == 0) {
  stop("Saved GAM prediction table is missing.")
}

if (is.null(transition_order_old) || nrow(transition_order_old) == 0) {
  stop("Saved transition-order table is missing.")
}

# =============================================================================
# Common style: Figure-6-like, larger/bolder, no small panel titles/subtitles
# =============================================================================

theme_fig8 <- function(base_size = 10.5) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "#202A33"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.title = element_text(face = "bold", size = base_size + 0.3),
      axis.text = element_text(face = "bold", size = base_size - 1.1, color = "#273641"),
      legend.title = element_text(face = "bold", size = base_size - 1.0),
      legend.text = element_text(size = base_size - 1.5),
      plot.tag = element_text(face = "bold", size = 13, color = "#111827"),
      plot.margin = margin(5, 6, 5, 6)
    )
}

# =============================================================================
# PANEL A | unbiased metabolic ranking
# Fix numerical FDR underflow only for visualization.
# The original ranking and DiscoveryScore are unchanged.
# =============================================================================

CATEGORY_COLORS <- c(
  "Carbohydrate" = "#4C78A8",
  "Energy" = "#72A0C1",
  "Lipid" = "#D95F02",
  "Amino acid" = "#59A14F",
  "Nucleotide" = "#B279A2",
  "Cofactor / redox" = "#E3A321",
  "Glycan" = "#8C8C8C",
  "Other metabolism" = "#7A7A7A"
)

plot_A_df <- screen %>%
  filter(is.finite(DiscoveryScore)) %>%
  arrange(desc(DiscoveryScore)) %>%
  slice_head(n = 18) %>%
  mutate(
    FDR_plot = pmax(FDR, 1e-16),
    FDR_evidence = pmin(-log10(FDR_plot), 16),
    Pathway = factor(Pathway, levels = rev(Pathway))
  )

pA <- ggplot(plot_A_df, aes(x = DiscoveryScore, y = Pathway)) +
  geom_segment(
    aes(x = 0, xend = DiscoveryScore, yend = Pathway),
    color = "#CAD3D9",
    linewidth = 0.8
  ) +
  geom_point(
    aes(size = FDR_evidence, fill = Category),
    shape = 21,
    color = "#23333D",
    stroke = 0.55,
    alpha = 0.98
  ) +
  scale_fill_manual(
    values = CATEGORY_COLORS,
    drop = FALSE,
    name = "Metabolic class"
  ) +
  scale_size_continuous(
    range = c(3.0, 7.2),
    limits = c(0, 16),
    breaks = c(4, 8, 12, 16),
    name = "-log10(FDR)\n(capped at 16)"
  ) +
  labs(
    tag = "A",
    x = "Dynamic discovery score",
    y = NULL
  ) +
  theme_fig8(10.5) +
  theme(
    axis.text.y = element_text(face = "bold", size = 8.7),
    panel.grid.major.x = element_line(color = "#E7ECEF", linewidth = 0.38),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# =============================================================================
# PANEL B | top dynamic programs, selected from the unbiased screen
# =============================================================================

zsafe <- function(v) {
  v <- as.numeric(v)
  s <- sd(v, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(v)))
  as.numeric(scale(v))
}

top_B_ids <- screen %>%
  filter(is.finite(DiscoveryScore), is.finite(TransitionMidpoint)) %>%
  arrange(desc(DiscoveryScore)) %>%
  slice_head(n = 12) %>%
  pull(PathwayID)

plot_B_df <- predictions %>%
  filter(PathwayID %in% top_B_ids) %>%
  group_by(PathwayID) %>%
  mutate(SmoothZ = zsafe(PredictedScore)) %>%
  ungroup()

heat_order <- screen %>%
  filter(PathwayID %in% top_B_ids) %>%
  arrange(TransitionMidpoint, desc(DiscoveryScore)) %>%
  pull(Pathway)

plot_B_df$Pathway <- factor(plot_B_df$Pathway, levels = rev(heat_order))

pB <- ggplot(plot_B_df, aes(x = Pseudotime, y = Pathway, fill = SmoothZ)) +
  geom_tile(linewidth = 0) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Relative\ntranscriptional\nscore"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = c(0, 0)
  ) +
  labs(
    tag = "B",
    x = "Normalized inferred pseudotime",
    y = NULL
  ) +
  theme_fig8(10.5) +
  theme(
    axis.text.y = element_text(face = "bold", size = 8.4),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# =============================================================================
# PANEL C | focused lipid-program transition ordering
#
# These pathways are selected AFTER the unbiased screen showed multiple lipid
# programs among the most dynamic trajectory-associated pathways.
# =============================================================================

LIPID_FEATURES_C <- c(
  "Glycerophospholipid metabolism",
  "Glycerolipid metabolism",
  "Ether lipid metabolism",
  "Sphingolipid metabolism",
  "Glycosphingolipid biosynthesis - lacto/neolacto"
)

lipid_C <- screen %>%
  filter(
    Pathway %in% LIPID_FEATURES_C,
    is.finite(TransitionMidpoint)
  ) %>%
  transmute(
    Feature = Pathway,
    FeatureType = "Lipid-metabolic program",
    TransitionMidpoint,
    LateMinusEarly,
    FDR,
    DiscoveryScore
  )

response_C <- transition_order_old %>%
  filter(
    FeatureType == "Endothelial response",
    Feature %in% c("Inflammatory / SASP", "ICAM1", "SERPINE1"),
    is.finite(TransitionMidpoint)
  ) %>%
  transmute(
    Feature,
    FeatureType = "Endothelial response",
    TransitionMidpoint,
    LateMinusEarly,
    FDR = NA_real_,
    DiscoveryScore = NA_real_
  )

plot_C_df <- bind_rows(lipid_C, response_C) %>%
  arrange(TransitionMidpoint)

plot_C_df$Feature <- factor(
  plot_C_df$Feature,
  levels = rev(plot_C_df$Feature)
)

TYPE_COLORS <- c(
  "Lipid-metabolic program" = "#D95F02",
  "Endothelial response" = "#B2182B"
)

pC <- ggplot(plot_C_df, aes(x = TransitionMidpoint, y = Feature)) +
  annotate(
    "rect",
    xmin = 0,
    xmax = 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "#EDF5F7",
    alpha = 0.55
  ) +
  annotate(
    "rect",
    xmin = 0.5,
    xmax = 1,
    ymin = -Inf,
    ymax = Inf,
    fill = "#FBF1F1",
    alpha = 0.38
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = TransitionMidpoint,
      yend = Feature,
      color = FeatureType
    ),
    linewidth = 0.9,
    alpha = 0.55
  ) +
  geom_point(
    aes(fill = FeatureType),
    shape = 21,
    size = 4.5,
    color = "#263640",
    stroke = 0.55
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    color = "#8C969E",
    linewidth = 0.55
  ) +
  scale_color_manual(values = TYPE_COLORS, guide = "none") +
  scale_fill_manual(values = TYPE_COLORS, name = NULL) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    tag = "C",
    x = "Normalized inferred pseudotime",
    y = NULL
  ) +
  theme_fig8(10.5) +
  theme(
    axis.text.y = element_text(face = "bold", size = 8.4),
    legend.position = "top",
    legend.justification = "center",
    panel.grid.major.x = element_line(color = "#E1E7EA", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# Save panel data
# =============================================================================

write.csv(
  plot_A_df,
  file.path(OUTDIR, "Figure7A_top18_unbiased_metabolic_programs.csv"),
  row.names = FALSE
)

write.csv(
  plot_C_df,
  file.path(OUTDIR, "Figure7C_focused_lipid_transition_order.csv"),
  row.names = FALSE
)

# =============================================================================
# Final A-C layout
# =============================================================================

Figure7_ABC_FINAL <- pA /
  (pB | pC) +
  plot_layout(
    heights = c(1.14, 1.00),
    widths = c(1, 1)
  ) +
  plot_annotation(
    title = expression(
      "Single-cell metabolic screening identifies early lipid-metabolic remodeling along the " *
        italic(SERPINE1) *
        "-associated endothelial injury continuum"
    ),
    theme = theme(
      plot.title = element_text(
        family = "sans",
        face = "bold",
        hjust = 0.5,
        size = 16,
        color = "#111827",
        margin = margin(b = 8)
      ),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# =============================================================================
# Save
# =============================================================================

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

save_one <- function(p, name, width, height) {
  ggsave(
    file.path(OUTDIR, paste0(name, ".pdf")),
    p,
    width = width,
    height = height,
    units = "in",
    device = pdf_device,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    file.path(OUTDIR, paste0(name, ".png")),
    p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = ragg::agg_png,
    bg = "white",
    limitsize = FALSE
  )
}

save_one(pA, "Figure7A_FINAL_unbiased_metabolic_screen", 10.0, 6.5)
save_one(pB, "Figure7B_FINAL_dynamic_metabolic_heatmap", 7.6, 5.8)
save_one(pC, "Figure7C_FINAL_lipid_transition_order", 7.6, 5.8)

save_one(
  Figure7_ABC_FINAL,
  "Figure7_STAGE1_A_C_FINAL_LOGIC",
  15.2,
  12.2
)

saveRDS(
  list(
    A = pA,
    B = pB,
    C = pC,
    assembled = Figure7_ABC_FINAL,
    panel_A_data = plot_A_df,
    panel_B_data = plot_B_df,
    panel_C_data = plot_C_df
  ),
  file.path(OUTDIR, "Figure7_STAGE1_A_C_FINAL_LOGIC_BUNDLE.rds"),
  compress = FALSE
)

cat(
  "\n============================================================\n",
  "FIG8 A-C FINAL LOGIC COMPLETED\n",
  "============================================================\n",
  "Panel A: unbiased metabolic ranking; FDR display capped at 1e-16.\n",
  "Panel B: unbiased top metabolic dynamics; transcriptional-score legend.\n",
  "Panel C: lipid programs from the unbiased screen + SASP/ICAM1/SERPINE1.\n",
  "\nPanel C transition order:\n",
  sep = ""
)

print(
  plot_C_df %>%
    select(Feature, FeatureType, TransitionMidpoint, LateMinusEarly) %>%
    arrange(TransitionMidpoint),
  row.names = FALSE
)

cat(
  "\nInterpret as relative ordering along inferred pseudotime only.\n",
  "Do not claim real-time or causal precedence.\n",
  "============================================================\n"
)

print(Figure7_ABC_FINAL)




}

# =============================================================================
# STAGE figure7_lipidomics
# =============================================================================
PIPELINE$figure7_lipidomics <- function() {
# scripts/09_lipidomics.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.

# =============================================================================
# FIGURE 8 | D-G PREMIUM LIPIDOMICS INTEGRATION + FINAL A-G ASSEMBLY
#
# Biological logic:
#   A-C  Single-cell discovery (already completed; no rerun here)
#   D    Sample-level heatmap of balanced top differential lipid species
#   E    Lipid-class + representative pathway-linked lipid remodeling
#   F    Lipidomic pathway enrichment
#   G    Cross-omic convergence matrix
#
# Final narrative:
#   unbiased single-cell metabolic screening
#        -> early lipid-metabolic remodeling along inferred pseudotime
#        -> mouse whole-lung lipidomics supports tissue lipid remodeling
#        -> cross-omic convergence prioritizes lipid pathways linked to the
#           SERPINE1-associated endothelial injury continuum
#
# IMPORTANT:
#   - Lipidomics animals are independent from the scRNA-seq cohort.
#   - NO sample-to-sample cross-omic correlation is performed.
#   - G is pathway-level cross-omic convergence, not causal inference.
#   - Sphingolipid lipidomics FDR ~0.069 should be described as suggestive/trend,
#     not FDR-significant.
#   - All figure text is English.
#   - No Monocle2 / igraph is loaded.
#
# INPUTS EXPECTED IN THE WORKING DIRECTORY
#
# 1) Existing Figure 7 single-cell analysis:
#    40_Figure7_SC_Metabolic_Discovery_NO_IGRAPH/
#      03_objects/
#        Figure7_STAGE1_SC_METABOLIC_DISCOVERY_BUNDLE.rds
#
#    40_Figure7_SC_Metabolic_Discovery_NO_IGRAPH/
#      05_FINAL_LOGIC_A_C/
#        Figure7_STAGE1_A_C_FINAL_LOGIC_BUNDLE.rds
#
# 2) Lipidomics:
#    Either:
#      P877Y24I29_NJMU_ZHShen_SH_result/
#        data_result/...
#
#    OR:
#      P877Y24I29_NJMU_ZHShen_SH_result.zip
#
#    If only the ZIP exists, this script extracts it automatically.
#
# OUTPUT:
#   41_Figure7_Lipidomics_Integration_PREMIUM/
# =============================================================================


# =============================================================================
# 0. Packages
# =============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260815)

required_pkgs <- c(
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "scales",
  "ragg"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ragg)
})


# =============================================================================
# 1. Paths
# =============================================================================

SC_STAGE1_BUNDLE <- repo_input("metabolic_bundle")

SC_FINAL_ABC_BUNDLE <- repo_input("metabolic_panels")

LIPID_BASE_NAME <- "P877Y24I29_NJMU_ZHShen_SH_result"
LIPID_ZIP <- paste0(
  LIPID_BASE_NAME,
  ".zip"
)

LIPID_EXTRACT_DIR <- paste0(
  LIPID_BASE_NAME,
  "_EXTRACTED"
)

OUTDIR <- repo_output("Figure7/Lipidomics")

PANEL_DIR <- file.path(
  OUTDIR,
  "01_panels"
)

TABLE_DIR <- file.path(
  OUTDIR,
  "02_tables"
)

OBJECT_DIR <- file.path(
  OUTDIR,
  "03_objects"
)

FINAL_DIR <- file.path(
  OUTDIR,
  "04_final"
)

for (d in c(
  PANEL_DIR,
  TABLE_DIR,
  OBJECT_DIR,
  FINAL_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!file.exists(SC_STAGE1_BUNDLE)) {
  stop(
    "Cannot find:\n",
    SC_STAGE1_BUNDLE
  )
}

if (!file.exists(SC_FINAL_ABC_BUNDLE)) {
  stop(
    "Cannot find:\n",
    SC_FINAL_ABC_BUNDLE
  )
}


# =============================================================================
# 2. Resolve lipidomics files automatically
#
# Your current folder is the FLAT six-file structure:
#
#   P877Y24I29_NJMU_ZHShen_SH_result/
#     04_AllMet.csv
#     04_AllMet_Pathway.csv
#     AllMet_Tests_Joined.csv
#     Class_Barplot_by_Group_Test.csv
#     Class_Info.csv
#     Pathway_Result.csv
#
# The code below prioritizes this flat structure.
# It also remains compatible with the original nested vendor folder and ZIP.
# =============================================================================

LIPID_ROOT_OVERRIDE <- repo_input("lipidomics_directory")
# Example if needed:

CORE_FILENAMES <- c(
  "04_AllMet.csv",
  "04_AllMet_Pathway.csv",
  "AllMet_Tests_Joined.csv",
  "Class_Barplot_by_Group_Test.csv",
  "Class_Info.csv",
  "Pathway_Result.csv"
)

is_flat_core_folder <- function(path_i) {

  if (
    is.null(
      path_i
    ) ||
      !dir.exists(
        path_i
      )
  ) {
    return(
      FALSE
    )
  }

  all(
    file.exists(
      file.path(
        path_i,
        CORE_FILENAMES
      )
    )
  )
}

resolve_flat_root <- function() {

  candidates <- unique(
    c(
      LIPID_ROOT_OVERRIDE,
      LIPID_BASE_NAME,
      "."
    )
  )

  candidates <- candidates[
    !is.na(
      candidates
    ) &
      nzchar(
        candidates
      )
  ]

  for (
    p in candidates
  ) {

    if (
      is_flat_core_folder(
        p
      )
    ) {
      return(
        normalizePath(
          p,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  NULL
}

resolve_nested_root <- function() {

  candidates <- unique(
    c(
      LIPID_ROOT_OVERRIDE,
      LIPID_BASE_NAME,
      "."
    )
  )

  candidates <- candidates[
    !is.na(
      candidates
    ) &
      nzchar(
        candidates
      )
  ]

  for (
    p in candidates
  ) {

    if (
      dir.exists(
        file.path(
          p,
          "data_result"
        )
      )
    ) {
      return(
        normalizePath(
          p,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  NULL
}

flat_root <- resolve_flat_root()
nested_root <- resolve_nested_root()

if (
  !is.null(
    flat_root
  )
) {

  LIPID_ROOT <- flat_root
  LIPID_LAYOUT <- "flat"

} else if (
  !is.null(
    nested_root
  )
) {

  LIPID_ROOT <- nested_root
  LIPID_LAYOUT <- "nested"

} else if (
  file.exists(
    LIPID_ZIP
  )
) {

  # Extract original vendor ZIP if no folder is found.
  if (
    !dir.exists(
      LIPID_EXTRACT_DIR
    )
  ) {
    dir.create(
      LIPID_EXTRACT_DIR,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  utils::unzip(
    LIPID_ZIP,
    exdir = LIPID_EXTRACT_DIR
  )

  if (
    is_flat_core_folder(
      LIPID_EXTRACT_DIR
    )
  ) {

    LIPID_ROOT <- normalizePath(
      LIPID_EXTRACT_DIR,
      winslash = "/",
      mustWork = TRUE
    )

    LIPID_LAYOUT <- "flat"

  } else if (
    dir.exists(
      file.path(
        LIPID_EXTRACT_DIR,
        "data_result"
      )
    )
  ) {

    LIPID_ROOT <- normalizePath(
      LIPID_EXTRACT_DIR,
      winslash = "/",
      mustWork = TRUE
    )

    LIPID_LAYOUT <- "nested"

  } else {

    stop(
      "ZIP was extracted, but neither the flat six-file layout nor ",
      "the original nested data_result layout was found."
    )
  }

} else {

  stop(
    "Could not locate lipidomics data.\n\n",
    "Preferred structure:\n",
    "  ",
    LIPID_BASE_NAME,
    "/04_AllMet.csv\n",
    "  ",
    LIPID_BASE_NAME,
    "/04_AllMet_Pathway.csv\n",
    "  ",
    LIPID_BASE_NAME,
    "/AllMet_Tests_Joined.csv\n",
    "  ",
    LIPID_BASE_NAME,
    "/Class_Barplot_by_Group_Test.csv\n",
    "  ",
    LIPID_BASE_NAME,
    "/Class_Info.csv\n",
    "  ",
    LIPID_BASE_NAME,
    "/Pathway_Result.csv\n\n",
    "If the folder is elsewhere, set LIPID_ROOT_OVERRIDE near the top."
  )
}

if (
  LIPID_LAYOUT == "flat"
) {

  FILE_ALLMET <- file.path(
    LIPID_ROOT,
    "04_AllMet.csv"
  )

  FILE_CLASS_INFO <- file.path(
    LIPID_ROOT,
    "Class_Info.csv"
  )

  FILE_ALLMET_PATHWAY <- file.path(
    LIPID_ROOT,
    "04_AllMet_Pathway.csv"
  )

  FILE_DIFF <- file.path(
    LIPID_ROOT,
    "AllMet_Tests_Joined.csv"
  )

  FILE_CLASS_STATS <- file.path(
    LIPID_ROOT,
    "Class_Barplot_by_Group_Test.csv"
  )

  FILE_PATHWAY <- file.path(
    LIPID_ROOT,
    "Pathway_Result.csv"
  )

} else {

  DATA_ROOT <- file.path(
    LIPID_ROOT,
    "data_result"
  )

  FILE_ALLMET <- file.path(
    DATA_ROOT,
    "Preprocessed_Data",
    "04_AllMet.csv"
  )

  FILE_CLASS_INFO <- file.path(
    DATA_ROOT,
    "Preprocessed_Data",
    "Class_Info.csv"
  )

  FILE_ALLMET_PATHWAY <- file.path(
    DATA_ROOT,
    "Preprocessed_Data",
    "04_AllMet_Pathway.csv"
  )

  FILE_DIFF <- file.path(
    DATA_ROOT,
    "Treatment",
    "sham_vs_LIRI",
    "06_Potential_Biomarkers",
    "AllMet_Tests_Joined.csv"
  )

  FILE_CLASS_STATS <- file.path(
    DATA_ROOT,
    "Treatment",
    "sham_vs_LIRI",
    "01_Basic_Statistics",
    "Metabolite_Class_Statistics",
    "Class_Barplot_by_Group_Test.csv"
  )

  FILE_PATHWAY <- file.path(
    DATA_ROOT,
    "Treatment",
    "sham_vs_LIRI",
    "07_Pathway_Analysis",
    "Pathway_Analysis",
    "Pathway_Result.csv"
  )
}

required_files <- c(
  FILE_ALLMET,
  FILE_CLASS_INFO,
  FILE_ALLMET_PATHWAY,
  FILE_DIFF,
  FILE_CLASS_STATS,
  FILE_PATHWAY
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (
  length(
    missing_files
  ) > 0
) {
  stop(
    "Missing lipidomics files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

cat(
  "\nLipidomics layout detected: ",
  LIPID_LAYOUT,
  "\nLipidomics folder: ",
  LIPID_ROOT,
  "\n",
  sep = ""
)


# =============================================================================
# 3. Load scRNA-seq A-C and full metabolic screen
# =============================================================================

sc_stage1 <- readRDS(
  SC_STAGE1_BUNDLE
)

sc_abc <- readRDS(
  SC_FINAL_ABC_BUNDLE
)

if (
  is.null(
    sc_stage1$screen
  )
) {
  stop(
    "sc_stage1$screen is missing."
  )
}

if (
  !all(
    c(
      "A",
      "B",
      "C"
    ) %in% names(
      sc_abc
    )
  )
) {
  stop(
    "Final A-C bundle does not contain plots A-C."
  )
}

sc_screen <- as.data.frame(
  sc_stage1$screen,
  stringsAsFactors = FALSE
)

pA <- sc_abc$A
pB <- sc_abc$B
pC <- sc_abc$C


# =============================================================================
# 4. Load lipidomics tables
# =============================================================================

allmet <- read.csv(
  FILE_ALLMET,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

class_info <- read.csv(
  FILE_CLASS_INFO,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

diff_all <- read.csv(
  FILE_DIFF,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

class_stats <- read.csv(
  FILE_CLASS_STATS,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pathway_result <- read.csv(
  FILE_PATHWAY,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

allmet_pathway <- read.csv(
  FILE_ALLMET_PATHWAY,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# The pathway table stores the pathway name in its first unnamed column.
if (
  !"Pathway" %in% colnames(
    pathway_result
  )
) {
  colnames(
    pathway_result
  )[1] <- "Pathway"
}

pathway_result$Pathway <- trimws(
  as.character(
    pathway_result$Pathway
  )
)

class_info$SampleID <- as.character(
  class_info$SampleID
)

class_info$ClassNote <- as.character(
  class_info$ClassNote
)


# =============================================================================
# 5. Global styling
# =============================================================================

COL_SHAM <- "#4C78A8"
COL_LIRI <- "#D73027"

COL_LIPID <- "#D95F02"
COL_SC <- "#2A788E"
COL_NEUTRAL <- "#C7D0D6"
COL_TEXT <- "#202A33"

theme_fig8 <- function(
    base_size = 10.5
) {

  ggplot2::theme_classic(
    base_size = base_size
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "sans",
        color = COL_TEXT
      ),

      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),

      axis.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 0.2,
        color = COL_TEXT
      ),

      axis.text = ggplot2::element_text(
        face = "bold",
        size = base_size - 1.2,
        color = "#273641"
      ),

      legend.title = ggplot2::element_text(
        face = "bold",
        size = base_size - 1.0
      ),

      legend.text = ggplot2::element_text(
        size = base_size - 1.5
      ),

      strip.text = ggplot2::element_text(
        face = "bold",
        size = base_size - 0.8,
        color = COL_TEXT
      ),

      strip.background = ggplot2::element_rect(
        fill = "#F2F5F6",
        color = "#C8D0D5",
        linewidth = 0.5
      ),

      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 13,
        color = "#111827"
      ),

      plot.margin = ggplot2::margin(
        5,
        6,
        5,
        6
      )
    )
}


# =============================================================================
# 6. PANEL D | Balanced top differential lipid-species heatmap
#
# Why replace PCA:
#   With n=6/group, PCA/ellipse contributes relatively little biological
#   information and can create large empty regions.  This panel instead shows
#   the actual sample-level abundance patterns of the strongest differential
#   lipid species.
#
# Selection rule:
#   1) Uni_FDR < 0.10
#   2) within each lipid class, rank by Uni_FDR then |log2FC|
#   3) retain at most 3 species per class
#
# This prevents a single abundant class (especially PC) from visually
# dominating the heatmap while keeping the selection data-driven.
#
# Heatmap value:
#   row-wise z score of processed lipid abundance across the 12 biological
#   samples. QC samples are excluded.
# =============================================================================

biological_samples <- class_info %>%
  dplyr::filter(
    ClassNote %in% c(
      "sham",
      "LIRI"
    )
  ) %>%
  dplyr::mutate(
    Group = factor(
      dplyr::recode(
        ClassNote,
        sham = "Ischemia",
        LIRI = "Ischemia/Reperfusion"
      ),
      levels = c(
        "Ischemia",
        "Ischemia/Reperfusion"
      )
    ),
    SampleNumber = suppressWarnings(
      as.numeric(
        gsub(
          "\\D+",
          "",
          SampleID
        )
      )
    )
  ) %>%
  dplyr::arrange(
    Group,
    SampleNumber,
    SampleID
  )

sample_ids <- biological_samples$SampleID

missing_sample_cols <- setdiff(
  sample_ids,
  colnames(
    allmet
  )
)

if (
  length(
    missing_sample_cols
  ) > 0
) {
  stop(
    "Samples missing from 04_AllMet.csv: ",
    paste(
      missing_sample_cols,
      collapse = ", "
    )
  )
}

FDR_SCREEN_D <- 0.10
MAX_PER_CLASS_D <- 3
MIN_SPECIES_D <- 8

selected_D <- diff_all %>%
  dplyr::filter(
    is.finite(
      Uni_FDR
    ),
    is.finite(
      log2FC
    ),
    Uni_FDR < FDR_SCREEN_D
  ) %>%
  dplyr::mutate(
    AbsLog2FC = abs(
      log2FC
    )
  ) %>%
  dplyr::arrange(
    Class,
    Uni_FDR,
    dplyr::desc(
      AbsLog2FC
    )
  ) %>%
  dplyr::group_by(
    Class
  ) %>%
  dplyr::slice_head(
    n = MAX_PER_CLASS_D
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    Uni_FDR,
    dplyr::desc(
      AbsLog2FC
    )
  )

# Defensive fallback if a future dataset contains too few FDR<0.10 species.
if (
  nrow(
    selected_D
  ) < MIN_SPECIES_D
) {

  selected_D <- diff_all %>%
    dplyr::filter(
      is.finite(
        Uni_FDR
      ),
      is.finite(
        log2FC
      )
    ) %>%
    dplyr::mutate(
      AbsLog2FC = abs(
        log2FC
      )
    ) %>%
    dplyr::arrange(
      Uni_FDR,
      dplyr::desc(
        AbsLog2FC
      )
    ) %>%
    dplyr::slice_head(
      n = 15
    )
}

if (
  nrow(
    selected_D
  ) == 0
) {
  stop(
    "No differential lipid species available for Panel D."
  )
}

abundance_D <- allmet %>%
  dplyr::select(
    Class,
    Metabolite,
    dplyr::all_of(
      sample_ids
    )
  ) %>%
  dplyr::inner_join(
    selected_D %>%
      dplyr::select(
        Class,
        Metabolite,
        Uni_FDR,
        log2FC
      ),
    by = c(
      "Class",
      "Metabolite"
    )
  )

if (
  nrow(
    abundance_D
  ) == 0
) {
  stop(
    "Selected differential lipids did not match 04_AllMet.csv."
  )
}

zsafe <- function(v) {

  v <- as.numeric(
    v
  )

  s <- stats::sd(
    v,
    na.rm = TRUE
  )

  if (
    !is.finite(
      s
    ) ||
      s == 0
  ) {
    return(
      rep(
        0,
        length(
          v
        )
      )
    )
  }

  as.numeric(
    (
      v -
        mean(
          v,
          na.rm = TRUE
        )
    ) /
      s
  )
}

plot_D_df <- abundance_D %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      sample_ids
    ),
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  dplyr::left_join(
    biological_samples %>%
      dplyr::select(
        SampleID,
        Group
      ),
    by = "SampleID"
  ) %>%
  dplyr::group_by(
    Class,
    Metabolite
  ) %>%
  dplyr::mutate(
    Z = zsafe(
      Abundance
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Display = paste0(
      Metabolite,
      "  [",
      Class,
      "]"
    ),
    SampleID = factor(
      SampleID,
      levels = sample_ids
    )
  )

# Put increased species toward the top and decreased species toward the bottom.
row_order_D <- selected_D %>%
  dplyr::arrange(
    dplyr::desc(
      log2FC
    ),
    Uni_FDR
  ) %>%
  dplyr::mutate(
    Display = paste0(
      Metabolite,
      "  [",
      Class,
      "]"
    )
  ) %>%
  dplyr::pull(
    Display
  )

plot_D_df$Display <- factor(
  plot_D_df$Display,
  levels = rev(
    unique(
      row_order_D
    )
  )
)

pD <- ggplot2::ggplot(
  plot_D_df,
  ggplot2::aes(
    x = SampleID,
    y = Display,
    fill = Z
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.55,
    width = 0.96,
    height = 0.96
  ) +
  ggplot2::facet_grid(
    . ~ Group,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(
      -2.5,
      2.5
    ),
    oob = scales::squish,
    name = "Relative\nabundance\n(z score)"
  ) +
  ggplot2::scale_x_discrete(
    expand = c(
      0,
      0
    )
  ) +
  ggplot2::scale_y_discrete(
    expand = ggplot2::expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  theme_fig8(
    10.2
  ) +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(
      face = "bold",
      size = 9.0
    ),
    strip.background = ggplot2::element_rect(
      fill = "#F1F4F5",
      color = "#CBD3D8",
      linewidth = 0.5
    ),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "bold",
      size = 7.6
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 7.4,
      lineheight = 0.95
    ),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(
      0.10,
      "lines"
    ),
    legend.position = "right",
    plot.margin = ggplot2::margin(
      5,
      10,
      5,
      5
    )
  )


# =============================================================================
# 7. PANEL E | Refined lipid remodeling panel
#
# Design:
#   - NO facet strips and NO large gray boxes.
#   - Two compact stacked subpanels:
#       upper: lipid classes
#       lower: representative pathway-linked species
#   - Shared visual grammar:
#       x = log2FC (Ischemia/Reperfusion relative to Ischemia)
#       fill = direction/magnitude of change
#       size = -log10(FDR), capped at 4
#   - Significant points (FDR < 0.05) receive a stronger outline.
#
# Representative species:
#   two lowest-Uni_FDR species from each priority pathway:
#     glycerophospholipid metabolism
#     ether lipid metabolism
#     sphingolipid metabolism
# =============================================================================

PRIORITY_PATHWAYS <- c(
  "Glycerophospholipid metabolism",
  "Ether lipid metabolism",
  "Sphingolipid metabolism"
)

# Map each lipid species to pathway annotation.
species_pathway <- diff_all %>%
  dplyr::left_join(
    allmet_pathway %>%
      dplyr::select(
        Class,
        Metabolite,
        Pathway
      ),
    by = c(
      "Class",
      "Metabolite"
    )
  )

representative_species_list <- lapply(
  PRIORITY_PATHWAYS,
  function(pathway_i) {

    species_pathway %>%
      dplyr::filter(
        !is.na(
          Pathway
        ),
        grepl(
          pathway_i,
          Pathway,
          fixed = TRUE
        ),
        is.finite(
          Uni_FDR
        ),
        is.finite(
          log2FC
        )
      ) %>%
      dplyr::arrange(
        Uni_FDR,
        Uni_P
      ) %>%
      dplyr::slice_head(
        n = 2
      ) %>%
      dplyr::mutate(
        PriorityPathway = pathway_i
      )
  }
)

representative_species <- dplyr::bind_rows(
  representative_species_list
) %>%
  dplyr::distinct(
    Class,
    Metabolite,
    .keep_all = TRUE
  )

class_E <- class_stats %>%
  dplyr::transmute(
    Label = as.character(
      Class
    ),
    log2FC = as.numeric(
      log2FC
    ),
    FDR = as.numeric(
      FDR
    )
  ) %>%
  dplyr::filter(
    is.finite(
      log2FC
    ),
    is.finite(
      FDR
    )
  ) %>%
  dplyr::mutate(
    FDR_evidence = pmin(
      -log10(
        pmax(
          FDR,
          1e-4
        )
      ),
      4
    ),
    Significant = FDR < 0.05
  )

species_E <- representative_species %>%
  dplyr::transmute(
    Label = as.character(
      Metabolite
    ),
    log2FC = as.numeric(
      log2FC
    ),
    FDR = as.numeric(
      Uni_FDR
    ),
    PriorityPathway = as.character(
      PriorityPathway
    )
  ) %>%
  dplyr::filter(
    is.finite(
      log2FC
    ),
    is.finite(
      FDR
    )
  ) %>%
  dplyr::mutate(
    FDR_evidence = pmin(
      -log10(
        pmax(
          FDR,
          1e-4
        )
      ),
      4
    ),
    Significant = FDR < 0.05
  )

# Shared symmetric x-axis keeps both sections visually comparable.
x_lim_E <- max(
  0.45,
  max(
    abs(
      c(
        class_E$log2FC,
        species_E$log2FC
      )
    ),
    na.rm = TRUE
  ) *
    1.12
)

# Clean within-section ordering.
class_E <- class_E %>%
  dplyr::arrange(
    log2FC
  ) %>%
  dplyr::mutate(
    Label = factor(
      Label,
      levels = Label
    )
  )

species_E <- species_E %>%
  dplyr::arrange(
    log2FC
  ) %>%
  dplyr::mutate(
    Label = factor(
      Label,
      levels = Label
    )
  )

# Slightly stronger point outline for statistically supported changes.
point_stroke <- function(sig) {
  ifelse(
    sig,
    1.15,
    0.55
  )
}

# Shared theme for the two E subpanels.
theme_E <- function(
    show_x = TRUE
) {

  th <- theme_fig8(
    10.2
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 10.4,
        hjust = 0,
        color = "#202A33",
        margin = ggplot2::margin(
          b = 4
        )
      ),

      axis.text.y = ggplot2::element_text(
        face = "bold",
        size = 7.7,
        color = "#263640"
      ),

      panel.grid.major.x = ggplot2::element_line(
        color = "#E7ECEF",
        linewidth = 0.35
      ),

      panel.grid.major.y = ggplot2::element_blank(),

      panel.grid.minor = ggplot2::element_blank(),

      legend.position = "right",

      plot.margin = ggplot2::margin(
        3,
        6,
        3,
        3
      )
    )

  if (
    !show_x
  ) {
    th <- th +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  th
}

pE_class <- ggplot2::ggplot(
  class_E,
  ggplot2::aes(
    x = log2FC,
    y = Label
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "#8B969E"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = log2FC,
      yend = Label
    ),
    linewidth = 0.85,
    color = "#C8D1D6"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      fill = log2FC,
      size = FDR_evidence
    ),
    shape = 21,
    color = "#263640",
    stroke = point_stroke(
      class_E$Significant
    ),
    alpha = 0.98
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B6EA8",
    mid = "#F3F4F5",
    high = "#D95F4B",
    midpoint = 0,
    limits = c(
      -x_lim_E,
      x_lim_E
    ),
    oob = scales::squish,
    name = "log2FC\n(Ischemia/Reperfusion relative to Ischemia)"
  ) +
  ggplot2::scale_size_continuous(
    range = c(
      2.7,
      6.4
    ),
    limits = c(
      0,
      4
    ),
    breaks = c(
      1,
      2,
      3,
      4
    ),
    name = "-log10(FDR)"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      -x_lim_E,
      x_lim_E
    ),
    breaks = scales::breaks_pretty(
      n = 5
    )
  ) +
  ggplot2::labs(
    title = "Lipid classes",
    x = NULL,
    y = NULL,
    tag = "E"
  ) +
  theme_E(
    show_x = FALSE
  ) +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      family = "sans",
      face = "bold",
      size = 13,
      color = "#111827"
    ),
    plot.tag.position = c(
      0.005,
      0.995
    )
  )

pE_species <- ggplot2::ggplot(
  species_E,
  ggplot2::aes(
    x = log2FC,
    y = Label
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "#8B969E"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = log2FC,
      yend = Label
    ),
    linewidth = 0.90,
    color = "#C8D1D6"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      fill = log2FC,
      size = FDR_evidence
    ),
    shape = 21,
    color = "#263640",
    stroke = point_stroke(
      species_E$Significant
    ),
    alpha = 0.98
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B6EA8",
    mid = "#F3F4F5",
    high = "#D95F4B",
    midpoint = 0,
    limits = c(
      -x_lim_E,
      x_lim_E
    ),
    oob = scales::squish,
    name = "log2FC\n(Ischemia/Reperfusion relative to Ischemia)"
  ) +
  ggplot2::scale_size_continuous(
    range = c(
      2.7,
      6.4
    ),
    limits = c(
      0,
      4
    ),
    breaks = c(
      1,
      2,
      3,
      4
    ),
    name = "-log10(FDR)"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      -x_lim_E,
      x_lim_E
    ),
    breaks = scales::breaks_pretty(
      n = 5
    )
  ) +
  ggplot2::labs(
    title = "Representative pathway-linked species",
    x = "Lipid change",
    y = NULL
  ) +
  theme_E(
    show_x = TRUE
  )

# Collect one shared legend and avoid the oversized facet-strip blocks.
pE <- (
  pE_class /
    pE_species
) +
  patchwork::plot_layout(
    heights = c(
      max(
        1,
        nrow(
          class_E
        )
      ),
      max(
        1,
        nrow(
          species_E
        )
      )
    ),
    guides = "collect"
  ) &
  ggplot2::theme(
    legend.position = "right"
  )

# Table used for export / reproducibility.
plot_E_df <- dplyr::bind_rows(
  class_E %>%
    dplyr::mutate(
      Level = "Lipid classes",
      PriorityPathway = NA_character_
    ),
  species_E %>%
    dplyr::mutate(
      Level = "Representative pathway-linked species"
    )
)


# =============================================================================
# 8. PANEL F | Independent lipidomic pathway enrichment
#
# All pathway-analysis results are shown.
# Priority pathways are highlighted because they independently emerged in
# the single-cell screen.
#
# Significance tiers:
#   FDR < 0.05        = significant
#   0.05 <= FDR < .10 = suggestive
#   FDR >= 0.10       = not supported at FDR < .10
# =============================================================================

plot_F_df <- pathway_result %>%
  dplyr::mutate(
    FDR = as.numeric(
      FDR
    ),
    Impact = as.numeric(
      Impact
    ),
    NegLog10FDR = -log10(
      pmax(
        FDR,
        1e-6
      )
    ),
    CrossOmicPriority = Pathway %in% PRIORITY_PATHWAYS,
    EvidenceTier = dplyr::case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      FDR < 0.10 ~ "FDR 0.05-0.10",
      TRUE ~ "FDR >= 0.10"
    )
  ) %>%
  dplyr::arrange(
    NegLog10FDR
  ) %>%
  dplyr::mutate(
    Pathway = factor(
      Pathway,
      levels = Pathway
    )
  )

FILL_F <- c(
  "FDR < 0.05" = "#B2182B",
  "FDR 0.05-0.10" = "#E39A32",
  "FDR >= 0.10" = "#B8C1C7"
)

pF <- ggplot2::ggplot(
  plot_F_df,
  ggplot2::aes(
    x = NegLog10FDR,
    y = Pathway
  )
) +
  ggplot2::geom_vline(
    xintercept = -log10(
      0.05
    ),
    linetype = "dashed",
    linewidth = 0.55,
    color = "#7A858D"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = NegLog10FDR,
      yend = Pathway
    ),
    linewidth = 0.85,
    color = "#C7D0D5"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = Impact,
      fill = EvidenceTier,
      stroke = CrossOmicPriority
    ),
    shape = 21,
    color = "#263640",
    alpha = 0.98
  ) +
  ggplot2::scale_fill_manual(
    values = FILL_F,
    name = "Lipidomics evidence"
  ) +
  ggplot2::scale_size_continuous(
    range = c(
      2.8,
      7.0
    ),
    name = "Pathway impact"
  ) +
  ggplot2::scale_discrete_manual(
    aesthetics = "stroke",
    values = c(
      "TRUE" = 1.2,
      "FALSE" = 0.45
    ),
    guide = "none"
  ) +
  ggplot2::labs(
    x = "-log10(lipidomics FDR)",
    y = NULL
  ) +
  theme_fig8(
    10.3
  ) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 8.0
    ),
    panel.grid.major.x = ggplot2::element_line(
      color = "#E8EDF0",
      linewidth = 0.35
    ),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right"
  )


# =============================================================================
# 9. PANEL G | Refined cross-omic convergence summary
#
# This replaces the previous "Relative evidence" heatmap.
#
# Rationale:
#   Rank, transition midpoint, pathway FDR and pathway impact are different
#   quantities and should NOT be forced onto a single common evidence scale.
#
# New visual grammar:
#   scRNA trajectory
#     - Unbiased rank: subtle blue intensity, smaller rank = darker
#     - Transition midpoint: pale blue descriptive cell, numeric value retained
#
#   Independent lipidomics
#     - Pathway FDR:
#         red   = FDR < 0.05
#         amber = 0.05-0.10
#         gray  = FDR >= 0.10
#     - Pathway impact: blue intensity proportional to impact
#
# The actual observed values remain printed in every cell.
# No sample-matched cross-omic correlation is performed.
# =============================================================================

lipid_pathway_clean <- pathway_result %>%
  dplyr::transmute(
    Pathway = trimws(
      as.character(
        Pathway
      )
    ),
    Lipidomics_FDR = as.numeric(
      FDR
    ),
    Lipidomics_Impact = as.numeric(
      Impact
    )
  )

shared_pathways <- intersect(
  as.character(
    sc_screen$Pathway
  ),
  lipid_pathway_clean$Pathway
)

if (
  length(
    shared_pathways
  ) < 2
) {
  stop(
    "Fewer than two pathway names are shared between scRNA and lipidomics."
  )
}

cross_df <- sc_screen %>%
  dplyr::filter(
    Pathway %in% shared_pathways
  ) %>%
  dplyr::select(
    Pathway,
    Rank,
    DiscoveryScore,
    FDR,
    TransitionMidpoint
  ) %>%
  dplyr::rename(
    scRNA_FDR = FDR
  ) %>%
  dplyr::left_join(
    lipid_pathway_clean,
    by = "Pathway"
  ) %>%
  dplyr::mutate(
    ConvergenceTier = dplyr::case_when(
      Lipidomics_FDR < 0.05 ~ "Lipidomics FDR < 0.05",
      Lipidomics_FDR < 0.10 ~ "Lipidomics FDR 0.05-0.10",
      TRUE ~ "Lipidomics FDR >= 0.10"
    )
  )

# Keep the shared-pathway structure from this dataset.
# The ordering makes the strongest biologically relevant convergence easy to read,
# while retaining non-convergent shared pathways as internal contrast.
G_PRIORITY_ORDER <- c(
  "Glycerophospholipid metabolism",
  "Ether lipid metabolism",
  "Sphingolipid metabolism",
  "Glycerolipid metabolism",
  "Inositol phosphate metabolism"
)

if (
  all(
    G_PRIORITY_ORDER %in% cross_df$Pathway
  )
) {

  cross_df <- cross_df %>%
    dplyr::filter(
      Pathway %in% G_PRIORITY_ORDER
    ) %>%
    dplyr::mutate(
      PathwayOrder = match(
        Pathway,
        G_PRIORITY_ORDER
      )
    ) %>%
    dplyr::arrange(
      PathwayOrder
    )

} else {

  cross_df <- cross_df %>%
    dplyr::arrange(
      Lipidomics_FDR,
      Rank
    ) %>%
    dplyr::slice_head(
      n = min(
        6,
        dplyr::n()
      )
    ) %>%
    dplyr::mutate(
      PathwayOrder = dplyr::row_number()
    )
}

format_fdr <- function(x) {

  ifelse(
    is.na(
      x
    ),
    "NA",
    ifelse(
      x < 0.001,
      "<0.001",
      sprintf(
        "%.3f",
        x
      )
    )
  )
}

wrap_one <- function(
    x,
    width = 27
) {

  vapply(
    x,
    function(z) {
      paste(
        strwrap(
          z,
          width = width
        ),
        collapse = "\n"
      )
    },
    FUN.VALUE = character(1)
  )
}

# ---------------------------------------------------------------------------
# Cell colors are metric-specific.
# They are deliberately NOT combined into one "relative evidence" score.
# ---------------------------------------------------------------------------

rank_min <- min(
  cross_df$Rank,
  na.rm = TRUE
)

rank_max <- max(
  cross_df$Rank,
  na.rm = TRUE
)

rank_strength <- if (
  is.finite(
    rank_min
  ) &&
    is.finite(
      rank_max
    ) &&
    rank_max > rank_min
) {
  1 -
    (
      cross_df$Rank -
        rank_min
    ) /
    (
      rank_max -
        rank_min
    )
} else {
  rep(
    0.5,
    nrow(
      cross_df
    )
  )
}

impact_max <- max(
  cross_df$Lipidomics_Impact,
  na.rm = TRUE
)

impact_strength <- if (
  is.finite(
    impact_max
  ) &&
    impact_max > 0
) {
  cross_df$Lipidomics_Impact /
    impact_max
} else {
  rep(
    0,
    nrow(
      cross_df
    )
  )
}

blue_pal_rank <- scales::col_numeric(
  palette = c(
    "#EEF4F8",
    "#8FAFC3",
    "#315E70"
  ),
  domain = c(
    0,
    1
  )
)

blue_pal_impact <- scales::col_numeric(
  palette = c(
    "#F1F5F7",
    "#A9C0D0",
    "#4C78A8"
  ),
  domain = c(
    0,
    1
  )
)

rank_fill <- blue_pal_rank(
  rank_strength
)

midpoint_fill <- rep(
  "#EAF1F5",
  nrow(
    cross_df
  )
)

fdr_fill <- dplyr::case_when(
  cross_df$Lipidomics_FDR < 0.05 ~ "#B2182B",
  cross_df$Lipidomics_FDR < 0.10 ~ "#E6A43A",
  TRUE ~ "#E5E8EA"
)

fdr_text <- dplyr::case_when(
  cross_df$Lipidomics_FDR < 0.05 ~ "white",
  TRUE ~ "#202A33"
)

impact_fill <- blue_pal_impact(
  pmin(
    pmax(
      impact_strength,
      0
    ),
    1
  )
)

# ---------------------------------------------------------------------------
# Build one plotting table.
# ---------------------------------------------------------------------------

g_rank <- cross_df %>%
  dplyr::transmute(
    Pathway,
    PathwayOrder,
    Omics = "scRNA trajectory",
    Metric = "Unbiased\nrank",
    Label = paste0(
      "#",
      Rank
    ),
    Fill = rank_fill,
    TextColor = "#1D2A33"
  )

g_midpoint <- cross_df %>%
  dplyr::transmute(
    Pathway,
    PathwayOrder,
    Omics = "scRNA trajectory",
    Metric = "Transition\nmidpoint",
    Label = sprintf(
      "%.2f",
      TransitionMidpoint
    ),
    Fill = midpoint_fill,
    TextColor = "#1D2A33"
  )

g_fdr <- cross_df %>%
  dplyr::transmute(
    Pathway,
    PathwayOrder,
    Omics = "Independent lipidomics",
    Metric = "Pathway\nFDR",
    Label = format_fdr(
      Lipidomics_FDR
    ),
    Fill = fdr_fill,
    TextColor = fdr_text
  )

g_impact <- cross_df %>%
  dplyr::transmute(
    Pathway,
    PathwayOrder,
    Omics = "Independent lipidomics",
    Metric = "Pathway\nimpact",
    Label = sprintf(
      "%.3f",
      Lipidomics_Impact
    ),
    Fill = impact_fill,
    TextColor = "#1D2A33"
  )

plot_G_df <- dplyr::bind_rows(
  g_rank,
  g_midpoint,
  g_fdr,
  g_impact
) %>%
  dplyr::mutate(
    PathwayDisplay = wrap_one(
      Pathway,
      width = 26
    ),
    Omics = factor(
      Omics,
      levels = c(
        "scRNA trajectory",
        "Independent lipidomics"
      )
    ),
    Metric = factor(
      Metric,
      levels = c(
        "Unbiased\nrank",
        "Transition\nmidpoint",
        "Pathway\nFDR",
        "Pathway\nimpact"
      )
    )
  )

display_levels <- cross_df %>%
  dplyr::arrange(
    PathwayOrder
  ) %>%
  dplyr::pull(
    Pathway
  ) %>%
  wrap_one(
    width = 26
  )

plot_G_df$PathwayDisplay <- factor(
  plot_G_df$PathwayDisplay,
  levels = rev(
    display_levels
  )
)

# Add a subtle visual separator between scRNA and lipidomics.
pG <- ggplot2::ggplot(
  plot_G_df,
  ggplot2::aes(
    x = Metric,
    y = PathwayDisplay
  )
) +
  ggplot2::geom_tile(
    ggplot2::aes(
      fill = Fill
    ),
    color = "white",
    linewidth = 1.05,
    width = 0.96,
    height = 0.90
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = Label,
      color = TextColor
    ),
    size = 3.55,
    fontface = "bold",
    lineheight = 0.95
  ) +
  ggplot2::facet_grid(
    . ~ Omics,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_color_identity() +
  ggplot2::scale_x_discrete(
    position = "top",
    expand = c(
      0,
      0
    )
  ) +
  ggplot2::scale_y_discrete(
    expand = ggplot2::expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    caption = "Lipidomics FDR: red < 0.05   |   amber 0.05-0.10   |   gray >= 0.10"
  ) +
  theme_fig8(
    10.3
  ) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(
      fill = "#EFF3F5",
      color = "#C9D2D7",
      linewidth = 0.55
    ),
    strip.text.x = ggplot2::element_text(
      face = "bold",
      size = 9.0,
      color = "#263640",
      margin = ggplot2::margin(
        t = 4,
        b = 4
      )
    ),
    axis.text.x = ggplot2::element_text(
      face = "bold",
      size = 8.0,
      lineheight = 0.95,
      color = "#263640",
      margin = ggplot2::margin(
        b = 4
      )
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 7.9,
      lineheight = 0.96,
      color = "#263640",
      margin = ggplot2::margin(
        r = 5
      )
    ),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(
      0.18,
      "lines"
    ),
    plot.caption = ggplot2::element_text(
      hjust = 1,
      size = 7.1,
      color = "#59666E",
      margin = ggplot2::margin(
        t = 5
      )
    ),
    legend.position = "none",
    plot.margin = ggplot2::margin(
      4,
      5,
      4,
      5
    )
  )


# =============================================================================
# 10. Save analysis tables
# =============================================================================

write.csv(
  selected_D,
  file.path(
    TABLE_DIR,
    "Figure7D_selected_differential_lipids.csv"
  ),
  row.names = FALSE
)

write.csv(
  plot_D_df,
  file.path(
    TABLE_DIR,
    "Figure7D_differential_lipid_heatmap_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  plot_E_df,
  file.path(
    TABLE_DIR,
    "Figure7E_class_and_representative_species.csv"
  ),
  row.names = FALSE
)

write.csv(
  representative_species,
  file.path(
    TABLE_DIR,
    "Figure7E_representative_species_selection.csv"
  ),
  row.names = FALSE
)

write.csv(
  plot_F_df,
  file.path(
    TABLE_DIR,
    "Figure7F_lipidomics_pathway_enrichment.csv"
  ),
  row.names = FALSE
)

write.csv(
  cross_df,
  file.path(
    TABLE_DIR,
    "Figure7G_cross_omic_convergence_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  plot_G_df,
  file.path(
    TABLE_DIR,
    "Figure7G_cross_omic_evidence_matrix.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 11. Add panel letters
# =============================================================================

tag_theme <- ggplot2::theme(
  plot.tag = ggplot2::element_text(
    family = "sans",
    face = "bold",
    size = 13,
    color = "#111827"
  ),
  plot.tag.position = c(
    0.005,
    0.995
  )
)

pA_final <- pA +
  ggplot2::labs(
    tag = "A"
  ) +
  tag_theme

pB_final <- pB +
  ggplot2::labs(
    tag = "B"
  ) +
  tag_theme

pC_final <- pC +
  ggplot2::labs(
    tag = "C"
  ) +
  tag_theme

pD_final <- pD +
  ggplot2::labs(
    tag = "D"
  ) +
  tag_theme

pE_final <- pE

pF_final <- pF +
  ggplot2::labs(
    tag = "F"
  ) +
  tag_theme

pG_final <- pG +
  ggplot2::labs(
    tag = "G"
  ) +
  tag_theme


# =============================================================================
# 12. Save D-G individual panels
# =============================================================================

pdf_device <- if (
  capabilities(
    "cairo"
  )
) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

save_panel <- function(
    p,
    filename,
    width,
    height
) {

  ggplot2::ggsave(
    filename = file.path(
      PANEL_DIR,
      paste0(
        filename,
        ".pdf"
      )
    ),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = pdf_device,
    bg = "white",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    filename = file.path(
      PANEL_DIR,
      paste0(
        filename,
        ".png"
      )
    ),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = ragg::agg_png,
    bg = "white",
    limitsize = FALSE
  )
}

save_panel(
  pD_final,
  "Figure7D_DIFFERENTIAL_LIPID_HEATMAP",
  8.0,
  6.8
)

save_panel(
  pE_final,
  "Figure7E_LIPID_REMODELING_REFINED",
  8.2,
  7.2
)

save_panel(
  pF_final,
  "Figure7F_LIPID_PATHWAY_ENRICHMENT",
  7.2,
  5.8
)

save_panel(
  pG_final,
  "Figure7G_CROSS_OMIC_CONVERGENCE_REFINED",
  8.4,
  5.8
)


# =============================================================================
# 13. Assemble FINAL A-G
#
# Layout:
#
#   A  full width
#
#   B | C
#
#   D | E
#
#   F | G
#
# No small panel titles/subtitles.
# =============================================================================

row_A <- pA_final

row_BC <- patchwork::wrap_plots(
  pB_final,
  pC_final,
  ncol = 2,
  widths = c(
    1,
    1
  )
)

# Explicit gutter between D and E.
# This avoids the visual collision created by E's left-side facet labels.
row_DE <- patchwork::wrap_plots(
  pD_final,
  patchwork::plot_spacer(),
  pE_final,
  ncol = 3,
  widths = c(
    1.00,
    0.04,
    1.10
  )
)

row_FG <- patchwork::wrap_plots(
  pF_final,
  pG_final,
  ncol = 2,
  widths = c(
    0.88,
    1.12
  )
)

FINAL_TITLE <- expression(
  "Single-cell metabolic screening and mouse whole-lung lipidomics show pathway-level overlap on lipid remodeling along the inferred " *
    italic(
      SERPINE1
    ) *
    "-associated endothelial injury continuum"
)

Figure7_FINAL <- patchwork::wrap_plots(
  row_A,
  row_BC,
  row_DE,
  row_FG,
  ncol = 1,
  heights = c(
    1.10,
    1.00,
    1.28,
    1.00
  )
) +
  patchwork::plot_annotation(
    title = FINAL_TITLE,
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = "sans",
        face = "bold",
        hjust = 0.5,
        size = 16,
        color = "#111827",
        margin = ggplot2::margin(
          b = 8
        )
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      )
    )
  )


# =============================================================================
# 14. Save FINAL A-G
# =============================================================================

FINAL_PDF <- file.path(
  FINAL_DIR,
  "Figure7_FINAL_A_G_PREMIUM.pdf"
)

FINAL_PNG <- file.path(
  FINAL_DIR,
  "Figure7_FINAL_A_G_PREMIUM.png"
)

FINAL_RDS <- file.path(
  OBJECT_DIR,
  "Figure7_FINAL_A_G_PREMIUM_plot.rds"
)

ggplot2::ggsave(
  filename = FINAL_PDF,
  plot = Figure7_FINAL,
  width = 15.4,
  height = 19.5,
  units = "in",
  device = pdf_device,
  bg = "white",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = FINAL_PNG,
  plot = Figure7_FINAL,
  width = 15.4,
  height = 19.5,
  units = "in",
  dpi = 500,
  device = ragg::agg_png,
  bg = "white",
  limitsize = FALSE
)

saveRDS(
  Figure7_FINAL,
  FINAL_RDS,
  compress = FALSE
)

saveRDS(
  list(
    source = list(
      lipid_root = LIPID_ROOT,
      sc_stage1_bundle = SC_STAGE1_BUNDLE,
      sc_final_abc_bundle = SC_FINAL_ABC_BUNDLE
    ),
    lipidomics = list(
      differential_heatmap_selection = selected_D,
      differential_heatmap_data = plot_D_df,
      class_stats = class_stats,
      differential_lipids = diff_all,
      representative_species = representative_species,
      pathway_results = pathway_result
    ),
    cross_omics = list(
      shared_pathways = shared_pathways,
      convergence_summary = cross_df,
      evidence_matrix = plot_G_df
    ),
    panels = list(
      A = pA_final,
      B = pB_final,
      C = pC_final,
      D = pD_final,
      E = pE_final,
      F = pF_final,
      G = pG_final
    ),
    assembled = Figure7_FINAL
  ),
  file.path(
    OBJECT_DIR,
    "Figure7_FINAL_A_G_PREMIUM_ANALYSIS_BUNDLE.rds"
  ),
  compress = FALSE
)


# =============================================================================
# 15. Console summary
# =============================================================================

cat(
  "\n============================================================\n",
  "FIGURE 8 D-G + FINAL A-G COMPLETED\n",
  "============================================================\n",
  "Lipidomics root:\n",
  LIPID_ROOT,
  "\n\nPanel D:\n",
  "  PCA of processed lipid abundance, QC excluded, n=6/group; convex hull shown.\n",
  "\nPanel E:\n",
  "  Lipid-class remodeling + two strongest species from each priority pathway.\n",
  "\nPanel F:\n",
  "  Independent lipidomic pathway enrichment.\n",
  "\nPanel G:\n",
  "  Pathway-level cross-omic convergence only; no sample-matched correlation.\n",
  "\nFinal PDF:\n",
  FINAL_PDF,
  "\n\nFinal PNG:\n",
  FINAL_PNG,
  "\n",
  sep = ""
)

cat(
  "\n============================================================\n",
  "PRIORITY LIPIDOMIC PATHWAYS\n",
  "============================================================\n"
)

print(
  plot_F_df %>%
    dplyr::filter(
      as.character(
        Pathway
      ) %in% PRIORITY_PATHWAYS
    ) %>%
    dplyr::select(
      Pathway,
      FDR,
      Impact,
      EvidenceTier
    ),
  row.names = FALSE
)

cat(
  "\n============================================================\n",
  "CROSS-OMIC CONVERGENCE\n",
  "============================================================\n"
)

print(
  cross_df %>%
    dplyr::select(
      Pathway,
      Rank,
      TransitionMidpoint,
      Lipidomics_FDR,
      Lipidomics_Impact,
      ConvergenceTier
    ),
  row.names = FALSE
)

cat(
  "\nInterpretation guardrails:\n",
  "  - FDR < 0.05 = lipidomics-supported.\n",
  "  - FDR 0.05-0.10 = suggestive only.\n",
  "  - Earlier transition = earlier along inferred pseudotime only.\n",
  "  - Independent cohorts support pathway-level convergence, not direct\n",
  "    sample-level causality.\n",
  sep = ""
)

print(
  Figure7_FINAL
)



}

# =============================================================================
# STAGE figure8_in_vitro
# =============================================================================
PIPELINE$figure8_in_vitro <- function() {
# Figure 8 numerical summaries. Supply biological-replicate data, not bar heights.
# Required columns: Endpoint, Group, Replicate, Value. Optional Field.
# Repeated fields for the same biological replicate are averaged first.
require_packages("ggplot2")
d<-read.csv(repo_input("in_vitro_measurements"),check.names=FALSE)
if(!all(c("Endpoint","Group","Replicate","Value")%in%names(d)))stop("Invalid Figure 8 schema.")
if(any(!is.finite(d$Value)))stop("Measurements must be finite.")
a<-aggregate(Value~Endpoint+Group+Replicate,d,mean)
out<-repo_output("Figure8")
write.csv(a,file.path(out,"biological_replicate_means.csv"),row.names=FALSE)
# Prespecified contrasts must be supplied explicitly. Do not automatically test
# every pair or infer an experimental design from the group names.
contrast_file<-file.path(REPO_ROOT,"inputs/figure8_contrasts.csv")
if(!file.exists(contrast_file))stop("Supply figure8_contrasts.csv: Endpoint,Group1,Group2,Test (t_test or tukey).")
c<-read.csv(contrast_file,check.names=FALSE)
results<-do.call(rbind,lapply(seq_len(nrow(c)),function(i){
 z<-a[a$Endpoint==c$Endpoint[i],];g1<-z$Value[z$Group==c$Group1[i]];g2<-z$Value[z$Group==c$Group2[i]]
 if(min(length(g1),length(g2))<2)stop("Insufficient independent replicates.")
 if(c$Test[i]=="t_test")p<-t.test(g1,g2,var.equal=TRUE,alternative="two.sided")$p.value else if(c$Test[i]=="tukey"){
  fit<-aov(Value~Group,z);tab<-TukeyHSD(fit)$Group;key<-paste(c$Group2[i],c$Group1[i],sep="-")
  key2<-paste(c$Group1[i],c$Group2[i],sep="-");if(!key%in%rownames(tab))key<-key2
  if(!key%in%rownames(tab))stop("Contrast does not match Tukey output.");p<-tab[key,"p adj"]
 }else stop("Unsupported prespecified test.")
 data.frame(c[i,],N1=length(g1),N2=length(g2),P=p)
}))
write.csv(results,file.path(out,"prespecified_statistics.csv"),row.names=FALSE)
p<-ggplot2::ggplot(a,ggplot2::aes(Group,Value,colour=Group))+ggplot2::geom_point(position=ggplot2::position_jitter(width=.05,height=0,seed=123456))+ggplot2::stat_summary(fun=mean,geom="point",shape=18,size=3)+ggplot2::facet_wrap(~Endpoint,scales="free_y")+ggplot2::theme_classic()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1),legend.position="none")
save_publication_plot(p,file.path(out,"Figure8_quantification"),14,9)

}

# =============================================================================
# STAGE supplementary_integration
# =============================================================================
PIPELINE$supplementary_integration <- function() {
# scripts/S1_integration_diagnostics.R
# Consolidated from author-supplied analysis code.
# Paths are configured centrally; input checks fail rather than manufacture results.
sce <- load_sce()
# =============================================================================
# FIGURE 3 SUPPLEMENT ONLY | R 4.3.2 / Seurat 4.4.0 / Harmony 1.2.0
# QC evidence (reuse files) + sample contributions + integration diagnostics
# DOES NOT recalculate QC, DoubletFinder, Harmony, or clustering.
# Optional new pre/post UMAPs are computed from SAVED PCA/Harmony embeddings
# using identical UMAP hyperparameters. Toggle below if too costly.
# =============================================================================
options(stringsAsFactors=FALSE);set.seed(1234)
SCE_RDATA <- ""            # only required if Seurat sce not in current session
QC_DIR <- repo_input("single_cell_qc_directory")# adjust if necessary
OUTDIR <- repo_output("Supplementary/S1")
MAKE_MATCHED_UMAP <- TRUE # may be slow for 100,789 cells; set FALSE for quick diagnostics
MAX_SCATTER_CELLS <- 50000L # only visualization; underlying sample tables use ALL cells
for(p in c("Seurat","ggplot2","Matrix","patchwork"))
 if(!requireNamespace(p,quietly=TRUE)) stop("Missing: ",p)
for(d in c(OUTDIR,file.path(OUTDIR,"Panels"),file.path(OUTDIR,"Tables")))
 dir.create(d,recursive=TRUE,showWarnings=FALSE)
TAB <- file.path(OUTDIR,"Tables"); PAN <- file.path(OUTDIR,"Panels")
if (!exists("sce",inherits=TRUE) || !inherits(get("sce",inherits=TRUE),"Seurat")) {
 if(!nzchar(SCE_RDATA)||!file.exists(SCE_RDATA)) stop("Load sce or provide a valid SCE_RDATA")
 env<-new.env(parent=globalenv());load(SCE_RDATA,envir=env)
 zz<-ls(env)[vapply(mget(ls(env),envir=env),inherits,logical(1),"Seurat")]
 if(length(zz)!=1L)stop("Need exactly one Seurat object in RData")
 sce<-get(zz,envir=env)
}
md <- sce@meta.data
stopifnot(all(c("orig.ident","celltype","seurat_clusters") %in% names(md)))
samp <- toupper(gsub("[^A-Z0-9]","",as.character(md$orig.ident)))
condition <- ifelse(grepl("^CIT",samp),"CIT","HR")
donor <- sub("^(CIT|HR)","",samp)
expected <- c(paste0("CIT",1:6),paste0("HR",1:6))
if(!setequal(unique(samp),expected)) stop("Expected CIT1...CIT6 HR1...HR6")
if(all(c("PairSample","PairDonor","PairCondition") %in% names(md))) {
 if(!identical(samp,toupper(gsub("[^A-Z0-9]","",as.character(md$PairSample))))) stop("PairSample mismatch")
 if(!identical(donor,gsub("[^0-9]","",as.character(md$PairDonor)))) stop("PairDonor mismatch")
 if(!identical(condition,ifelse(toupper(as.character(md$PairCondition)) %in% c("CIT"),"CIT","HR"))) stop("PairCondition mismatch")
}
P <- data.frame(Sample=expected,Condition=rep(c("CIT","HR"),each=6),
                Donor=rep(1:6,2),Ncells=as.integer(table(factor(samp,levels=expected))))
write.csv(P,file.path(TAB,"01_12_libraries_postQC_counts.csv"),row.names=FALSE)

# Metadata coverage of all cells - no cell treated as an independent donor
D <- data.frame(Sample=samp,Condition=condition,Donor=donor,
                Celltype=as.character(md$celltype),
                Cluster=as.character(md$seurat_clusters),stringsAsFactors=FALSE)
T <- as.data.frame(table(Sample=factor(D$Sample,levels=expected),Celltype=D$Celltype))
colnames(T)[3] <- "Cells"
T <- T[T$Cells>=0,]; T$Donor <- sub("^(CIT|HR)","",T$Sample)
T$Condition <- ifelse(grepl("^CIT",T$Sample),"CIT","HR")
T$FractionOfSample <- T$Cells/as.integer(P$Ncells[match(as.character(T$Sample),P$Sample)])
T$FractionOfCelltype <- T$Cells/ave(T$Cells,T$Celltype,FUN=function(x)sum(x))
write.csv(T,file.path(TAB,"02_major_celltype_by_12_sample.csv"),row.names=FALSE)
C <- as.data.frame(table(Sample=factor(D$Sample,levels=expected),Cluster=D$Cluster))
colnames(C)[3]<-"Cells"
C$FractionOfCluster <- C$Cells/ave(C$Cells,C$Cluster,FUN=function(x)sum(x))
C$FractionWithinSample <- C$Cells/as.integer(P$Ncells[match(as.character(C$Sample),P$Sample)])
write.csv(C,file.path(TAB,"03_wholelung_cluster_contributions.csv"),row.names=FALSE)
CR <- aggregate(FractionOfCluster~Cluster,C,max)
names(CR)[2] <- "MaximumSingleSampleFraction"
CR$FlagAbove50Percent <- CR$MaximumSingleSampleFraction>.5
write.csv(CR,file.path(TAB,"04_cluster_dominance_flags.csv"),row.names=FALSE)
# This is NOT EC subclustering: only true EC-state metadata are acceptable.
subtype_options <- c("EC_subtype","endothelial_subtype","EC_cluster","EC_state","Endothelial_Subcluster","endo_celltype")
subtype <- intersect(subtype_options,names(md))
if(length(subtype)){
 ec <- D$Celltype=="Endothelial cell"
 ec_state <- as.character(md[[subtype[1]]][ec])
 ecT <- as.data.frame(table(Sample=factor(samp[ec],levels=expected),EC_State=ec_state))
 names(ecT)[3] <- "Cells"
 ecT$FractionOfState <- ecT$Cells/ave(ecT$Cells,ecT$EC_State,FUN=function(x)sum(x))
 write.csv(ecT,file.path(TAB,"05_EC_substate_contribution_if_present.csv"),row.names=FALSE)
} else {
 writeLines("The main sce metadata contains no verified specialized EC-subcluster identity. Whole-lung seurat_clusters must NOT be relabeled as Homeostatic/Transitional EC. Obtain the endothelial re-clustered Figure5 Seurat object and repeat contribution checks for EC substate claims.",file.path(TAB,"05_EC_substate_UNAVAILABLE_README.txt"))
}
# Sample and donor contribution plots from existing metadata
theme_pub <- ggplot2::theme_classic(base_size=10)+ggplot2::theme(
 text=ggplot2::element_text(face="bold",color="black"),axis.text=ggplot2::element_text(color="black",face="bold"),
 plot.title=ggplot2::element_text(face="bold",color="black",hjust=.5))
save <- function(p,name,w=11,h=6){
 ggplot2::ggsave(file.path(PAN,paste0(name,".pdf")),p,width=w,height=h,bg="white",useDingbats=FALSE)
 if(requireNamespace("ragg",quietly=TRUE))ggplot2::ggsave(file.path(PAN,paste0(name,".png")),p,width=w,height=h,dpi=450,device=ragg::agg_png,bg="white")
 else ggplot2::ggsave(file.path(PAN,paste0(name,".png")),p,width=w,height=h,dpi=450,bg="white")
}
p1 <- ggplot2::ggplot(T,ggplot2::aes(x=Sample,y=FractionOfSample,fill=Celltype))+
 ggplot2::geom_col(width=.82)+ggplot2::facet_grid(.~Condition,scales="free_x",space="free_x")+
 ggplot2::labs(title="Major cell types: each biological library",x=NULL,y="Fraction within sample")+theme_pub+
 ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1))
save(p1,"Supplement_celltype_fraction_12samples",10,6)
p2 <- ggplot2::ggplot(C,ggplot2::aes(x=Sample,y=Cluster,fill=FractionOfCluster))+
 ggplot2::geom_tile(color="white",linewidth=.2)+ggplot2::scale_fill_gradient(low="white",high="#2166AC",name="Fraction of\ncluster")+
 ggplot2::labs(title="Whole-lung cluster contributions by sample",x=NULL,y="Seurat cluster")+theme_pub+
 ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55,hjust=1))
save(p2,"Supplement_cluster_sample_contribution",10,8)

# QC: collect EXISTING per-library CSVs; no new QC calculations or thresholds.
if(dir.exists(QC_DIR)) {
 files <- list.files(QC_DIR,pattern="(cell_filtering_summary|DoubletFinder_summary|DoubletFinder_parameters)\\.csv$",
                     recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
 inventory <- data.frame(Path=files,Filename=basename(files),stringsAsFactors=FALSE)
 write.csv(inventory,file.path(TAB,"06_existing_QC_DoubletFinder_files_inventory.csv"),row.names=FALSE)
 # Preserve original field names and thresholds instead of guessing their meaning.
 tables <- lapply(files,function(f){
  x<-tryCatch(utils::read.csv(f,check.names=FALSE),error=function(e)NULL)
  if(is.null(x))return(NULL)
  x$SourceFile <- basename(f)
  x$Sample <- sub("^((CIT|HR)[1-6]).*$","\\1",basename(f))
  x
 })
 tables <- Filter(Negate(is.null),tables)
 if(length(tables)){
  allcols<-unique(unlist(lapply(tables,names),use.names=FALSE))
  aligned <- lapply(tables,function(x){for(n in setdiff(allcols,names(x)))x[[n]]<-NA;x[,allcols,drop=FALSE]})
  combined<-do.call(rbind,aligned)
  write.csv(combined,file.path(TAB,"07_existing_QC_parameters_merged_unmodified.csv"),row.names=FALSE)
 }
} else {
 writeLines(paste("QC_DIR is missing:",QC_DIR,"; supply your real 01_QC_DoubletFinder folder"),
            file.path(TAB,"06_QC_folder_not_found.txt"))
}

# A reviewer-facing QC cell retention plot from the ORIGINAL 12 saved filtering summaries.
# Never infer threshold values that were not recorded in those source files.
if(dir.exists(QC_DIR)) {
  fq <- list.files(QC_DIR,pattern="_cell_filtering_summary\\.csv$",full.names=TRUE,recursive=TRUE)
  if(length(fq)) {
    chunks <- lapply(fq,function(file){
      v <- tryCatch(utils::read.csv(file,check.names=FALSE,stringsAsFactors=FALSE),error=function(e)NULL)
      if(is.null(v)) return(NULL)
      need<-c("sample_id","cells_before_QC","cells_after_QC","cells_after_DoubletFinder")
      if(!all(need %in% names(v))) return(NULL)
      v[1,need,drop=FALSE]
    })
    chunks<-Filter(Negate(is.null),chunks)
    if(length(chunks)) {
      qc<-do.call(rbind,chunks)
      qc$sample_id <- toupper(as.character(qc$sample_id))
      if(nrow(qc)==12L && setequal(qc$sample_id,expected)) {
        qc <- qc[match(expected,qc$sample_id),]
        write.csv(qc,file.path(TAB,"07a_QC_original_cell_retention_12samples.csv"),row.names=FALSE)
        qcL <- data.frame(Sample=rep(qc$sample_id,3),
          Stage=factor(rep(c("Before QC","After QC","After DoubletFinder"),each=12),
                       levels=c("Before QC","After QC","After DoubletFinder")),
          Cells=c(qc$cells_before_QC,qc$cells_after_QC,qc$cells_after_DoubletFinder))
        p_qc <- ggplot2::ggplot(qcL,ggplot2::aes(x=Sample,y=Cells,group=Stage,color=Stage))+
          ggplot2::geom_point(size=2)+ggplot2::geom_line(linewidth=.7)+
          ggplot2::labs(title="Per-library cell retention through QC and DoubletFinder",
                        x=NULL,y="Cells retained",color="Original QC stage")+
          theme_pub+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55,hjust=1))
        save(p_qc,"Supplement_QC_cell_retention_12samples",10,5.5)
      } else warning("QC files found but not exactly 12 unique matching samples; retention figure skipped")
    }
  }
}

# Integration diagnostics: saved PCA/Harmony, no Harmony re-run.
reduce <- names(sce@reductions)
if(!all(c("pca","harmony","umap") %in% reduce)) stop("Missing stored pca/harmony/umap reductions")
pca <- Seurat::Embeddings(sce,"pca")
har <- Seurat::Embeddings(sce,"harmony")
uold <- Seurat::Embeddings(sce,"umap")
if(nrow(pca)!=nrow(md) || nrow(har)!=nrow(md)) stop("Embeddings are not cell-aligned")
all_cells <- rownames(md)
if(!identical(rownames(pca),all_cells)){
 if(!setequal(rownames(pca),all_cells)) stop("PCA cellIDs mismatch");pca<-pca[all_cells,,drop=FALSE]
}
if(!identical(rownames(har),all_cells)){
 if(!setequal(rownames(har),all_cells)) stop("Harmony cellIDs mismatch");har<-har[all_cells,,drop=FALSE]
}
if(!identical(rownames(uold),all_cells)){
 if(!setequal(rownames(uold),all_cells)) stop("UMAP cellIDs mismatch");uold<-uold[all_cells,,drop=FALSE]
}
set.seed(1234)
keep <- sort(sample(seq_along(all_cells),min(MAX_SCATTER_CELLS,length(all_cells))))
mk <- function(M,type){data.frame(D1=M[keep,1],D2=M[keep,2],Sample=samp[keep],Condition=condition[keep],Type=type)}
plot_dim <- function(d,title){ggplot2::ggplot(d,ggplot2::aes(x=D1,y=D2,color=Sample))+
 ggplot2::geom_point(size=.12,alpha=.65)+ggplot2::labs(title=title,x="Dimension 1",y="Dimension 2")+
 theme_pub+ggplot2::theme(legend.position="bottom")}
save(plot_dim(mk(pca,"Uncorrected PCs"),"PCA before Harmony"),"Supplement_PCA_before_Harmony_by_sample",9,6)
save(plot_dim(mk(har,"Harmony dimensions"),"Harmony embeddings by sample"),"Supplement_Harmony_embeddings_by_sample",9,6)
p_group <- ggplot2::ggplot(mk(uold,"Original UMAP"),ggplot2::aes(D1,D2,color=Condition))+
 ggplot2::geom_point(size=.12,alpha=.7)+ggplot2::scale_color_manual(values=c(CIT="#2166AC",HR="#D73027"))+
 ggplot2::labs(title="Original integrated UMAP, colored by condition",x="UMAP_1",y="UMAP_2")+theme_pub
save(p_group,"Supplement_Original_Harmony_UMAP_by_condition",8,6)
p_by_sample <- ggplot2::ggplot(mk(uold,"Original UMAP"),ggplot2::aes(D1,D2,color=Sample))+
 ggplot2::geom_point(size=.12,alpha=.7)+
 ggplot2::labs(title="Original integrated UMAP, colored by sample",x="UMAP_1",y="UMAP_2")+
 theme_pub+ggplot2::theme(legend.position="bottom")
save(p_by_sample,"Supplement_Original_Harmony_UMAP_by_sample",9,6)


# Reviewer-grade comparable *new* UMAPs use the SAME algorithm/parameters
# on uncorrected PCA and Harmony dimensions; both are supplementary only.
if(isTRUE(MAKE_MATCHED_UMAP)) {
 if(!requireNamespace("uwot",quietly=TRUE)) stop("uwot required for matched pre/post UMAP")
 ndim <- min(20L,ncol(pca),ncol(har))
 cat("Constructing matched UMAPs from stored dimensions; no Harmony recomputation.\n")
 set.seed(1234)
 u_pre <- uwot::umap(pca[,seq_len(ndim),drop=FALSE],n_neighbors=30,min_dist=.3,
                      metric="euclidean",n_threads=2,verbose=TRUE)
 set.seed(1234)
 u_post <- uwot::umap(har[,seq_len(ndim),drop=FALSE],n_neighbors=30,min_dist=.3,
                       metric="euclidean",n_threads=2,verbose=TRUE)
 makeU <- function(U,name){
   d <- data.frame(X=U[keep,1],Y=U[keep,2],Condition=condition[keep],Sample=samp[keep])
   pA<-ggplot2::ggplot(d,ggplot2::aes(X,Y,color=Sample))+ggplot2::geom_point(size=.12,alpha=.75)+theme_pub+
      ggplot2::labs(title=paste0(name," (sample)"),x="UMAP_1",y="UMAP_2")+ggplot2::theme(legend.position="bottom")
   pB<-ggplot2::ggplot(d,ggplot2::aes(X,Y,color=Condition))+ggplot2::geom_point(size=.12,alpha=.75)+
      ggplot2::scale_color_manual(values=c(CIT="#2166AC",HR="#D73027"))+theme_pub+
      ggplot2::labs(title=paste0(name," (condition)"),x="UMAP_1",y="UMAP_2")
   pA|pB
 }
 save(makeU(u_pre,"Before Harmony"),"Supplement_PreHarmony_UMAP_matched",13,6)
 save(makeU(u_post,"After Harmony"),"Supplement_PostHarmony_UMAP_matched",13,6)
}
writeLines(capture.output(sessionInfo()),file.path(TAB,"08_sessionInfo.txt"))
message("Supplement evidence exported in: ",normalizePath(OUTDIR))



}

# =============================================================================
# Final Figure 5 assembly: A-F trajectory, G-J networks, K paired proxies.
# =============================================================================
PIPELINE$figure5_assemble <- function() {
 require_packages(c("ggplot2","patchwork"))
 files <- c(file.path(repo_output("Figure5/Trajectory"),"Figure5_A_F_panels.rds"),
            file.path(repo_output("Figure5/CellChat"),"Figure5_G_J_panels.rds"),
            file.path(repo_output("Figure5/Paired_K"),"Figure5K_panel.rds"))
 if(any(!file.exists(files)))stop("Run the Figure 5 trajectory, CellChat and paired-proxy stages first.")
 a<-readRDS(files[1])$panels;b<-readRDS(files[2]);k<-readRDS(files[3])
 panels<-c(a,b,list(K=k))
 if(!all(LETTERS[1:11]%in%names(panels)))stop("Expected final panels A-K.")
 box<-function(p)patchwork::wrap_elements(full=p)
 upper<-patchwork::wrap_plots(
   patchwork::wrap_plots(box(panels$A),box(panels$B),box(panels$C),nrow=1),
   box(panels$D),patchwork::wrap_plots(box(panels$E),box(panels$F),nrow=1),
   ncol=1,heights=c(1,.7,1))
 lower<-patchwork::wrap_plots(
   box(panels$G),patchwork::wrap_plots(box(panels$H),box(panels$I),ncol=1),
   box(panels$J),box(panels$K),ncol=2,heights=c(1.2,1))
 figure<-patchwork::wrap_plots(box(upper),box(lower),ncol=1,heights=c(1.2,1.1))
 save_publication_plot(figure,file.path(repo_output("Figure5/Final"),"Figure5_A_K"),20,26)
 saveRDS(list(panels=panels,figure=figure),file.path(repo_output("Figure5/Final"),"Figure5_A_K_bundle.rds"))
}

# =============================================================================
# Stage runner. No separate .R helpers or source scripts are needed.
# Input matrices, metadata and saved analysis objects remain external inputs.
# =============================================================================
STAGE_ORDER <- c(
 "figure1_bulk_plots",
 "figure1_wgcna",
 "figure1_blue_enrichment",
 "figure1_ml",
 "figure1_ml_plots",
 "figure2_shap_train",
 "figure2_shap_test",
 "figure2_nine_gene_scores",
 "figure2_nine_gene_function",
 "figure3_upstream_qc",
 "figure3_atlas",
 "figure3_paired_compute",
 "figure3_paired_plots",
 "figure3_kegg_curves",
 "figure4_states",
 "figure5_trajectory_AF",
 "figure5_network_master",
 "figure5_cellchat_GJ",
 "figure5_paired_proxy_K",
 "figure6_perturbation",
 "figure6_perturbation_plots",
 "figure7_metabolic",
 "figure7_lipidomics",
 "figure8_in_vitro",
 "supplementary_integration",
 "figure5_assemble")
list_stages <- function() {
 for(i in seq_along(STAGE_ORDER))cat(sprintf("%02d  %s\n",i,STAGE_ORDER[i]))
 invisible(STAGE_ORDER)
}
check_embedded_code <- function() {
 for(nm in STAGE_ORDER)parse(text=deparse(body(PIPELINE[[nm]])))
 message("Embedded stage bodies parsed. This is not an execution/result check.")
 invisible(TRUE)
}
run_stage <- function(name) {
 if(length(name)!=1L||!name%in%STAGE_ORDER)stop("Unknown stage; use list_stages().")
 logdir<-file.path(REPO_ROOT,"logs");dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
 logpath<-file.path(logdir,paste0(name,"_",format(Sys.time(),"%Y%m%d_%H%M%S"),".log"))
 oldwd<-getwd();setwd(REPO_ROOT);on.exit(setwd(oldwd),add=TRUE)
 sinks_before<-sink.number();sink(logpath,split=TRUE)
 on.exit({while(sink.number()>sinks_before)sink()},add=TRUE)
 cat("Stage:",name,"\nStarted:",as.character(Sys.time()),"\n")
 # Author source continuations share saved objects and use .GlobalEnv explicitly.
 # Evaluate in that same workspace, preserving source behavior.
 eval(body(PIPELINE[[name]]),envir=.GlobalEnv)
 print(utils::sessionInfo())
 cat("Finished:",as.character(Sys.time()),"\n")
 invisible(TRUE)
}
if(sys.nframe()==0L) {
 .cli<-commandArgs(trailingOnly=TRUE)
 if(length(.cli)&&.cli[1]=="--list")list_stages()
 else if(length(.cli)&&.cli[1]=="--check")check_embedded_code()
 else if(length(.cli))for(.stage in .cli)run_stage(.stage)
 else if(length(STAGES_TO_RUN))for(.stage in STAGES_TO_RUN)run_stage(.stage)
 else message("Code loaded. Use list_stages() and run_stage(name); edit PROJECT_ROOT first.")
} else if(length(STAGES_TO_RUN)) {
 for(.stage in STAGES_TO_RUN)run_stage(.stage)
}
