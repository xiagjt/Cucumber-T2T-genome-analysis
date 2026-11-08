#!/usr/bin/env Rscript

# Signed, two-tailed ICA gene-module assignment for picaplot::runICA().
#
# Usage:
# Rscript strategy1_ica_signed_modules.R \
#   eQTL.expression.filtered.xls ICA_signed_results \
#   20 2.0 0.5 0.30 20 1 80 15 10 1234875
#
# Positional arguments:
#   1  expression_file   genes x samples, first column contains gene IDs
#   2  output_dir
#   3  n_cores
#   4  loading_z_cutoff  minimum absolute robust loading Z-score
#   5  min_delta_z       best minus second-best absolute loading Z-score
#   6  min_kME           minimum positive kME after signed-module PCA
#   7  min_module_size
#   8  log_offset        expression transform is log2(x + log_offset)
#   9  var_cutoff        passed to picaplot::runICA
#  10  n_runs           passed to picaplot::runICA
#  11  max_iter         passed to picaplot::runICA
#  12  seed

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(picaplot)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(i, default, mode = c("character", "numeric", "integer")) {
  mode <- match.arg(mode)
  value <- if (length(args) >= i && nzchar(args[i])) args[i] else default
  if (mode == "numeric") return(as.numeric(value))
  if (mode == "integer") return(as.integer(value))
  as.character(value)
}

expression_file <- get_arg(1, "eQTL.expression.filtered.xls")
output_dir <- get_arg(2, "ICA_signed_results")
n_cores <- get_arg(3, 20L, "integer")
loading_z_cutoff <- get_arg(4, 2.0, "numeric")
min_delta_z <- get_arg(5, 0.5, "numeric")
min_kME <- get_arg(6, 0.30, "numeric")
min_module_size <- get_arg(7, 20L, "integer")
log_offset <- get_arg(8, 1.0, "numeric")
var_cutoff <- get_arg(9, 80, "numeric")
n_runs <- get_arg(10, 15L, "integer")
max_iter <- get_arg(11, 10L, "integer")
seed <- get_arg(12, 1234875L, "integer")

if (!file.exists(expression_file)) {
  stop("Expression file not found: ", expression_file, call. = FALSE)
}
if (!is.finite(n_cores) || n_cores < 1L) stop("n_cores must be >= 1")
if (!is.finite(loading_z_cutoff) || loading_z_cutoff <= 0) {
  stop("loading_z_cutoff must be > 0")
}
if (!is.finite(min_delta_z) || min_delta_z < 0) stop("min_delta_z must be >= 0")
if (!is.finite(min_kME) || min_kME < -1 || min_kME > 1) {
  stop("min_kME must be between -1 and 1")
}
if (!is.finite(min_module_size) || min_module_size < 2L) {
  stop("min_module_size must be >= 2")
}
if (!is.finite(log_offset) || log_offset <= 0) stop("log_offset must be > 0")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path, row_names = FALSE, col_names = TRUE) {
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = row_names,
    col.names = col_names, na = "NA"
  )
}

robust_z <- function(x) {
  center <- median(x, na.rm = TRUE)
  spread <- mad(x, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) spread <- sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
  (x - center) / spread
}

calculate_eigengene <- function(gene_by_sample) {
  sample_by_gene <- t(gene_by_sample)
  if (ncol(sample_by_gene) == 1L) {
    score <- as.numeric(scale(sample_by_gene[, 1]))
    kme <- 1
    names(kme) <- colnames(sample_by_gene)
    names(score) <- rownames(sample_by_gene)
    return(list(score = score, kME = kme))
  }

  pc <- prcomp(sample_by_gene, center = TRUE, scale. = TRUE, rank. = 1)
  score <- pc$x[, 1]
  kme <- as.numeric(cor(
    sample_by_gene, score, use = "pairwise.complete.obs",
    method = "pearson"
  ))
  names(kme) <- colnames(sample_by_gene)

  # PC sign is arbitrary. Orient it so that most genes have positive kME.
  if (median(kme, na.rm = TRUE) < 0) {
    score <- -score
    kme <- -kme
  }
  names(score) <- rownames(sample_by_gene)
  list(score = score, kME = kme)
}

cat("Reading expression matrix ...\n")
expr_df <- read.delim(
  expression_file, header = TRUE, row.names = 1, check.names = FALSE,
  quote = "", comment.char = "", na.strings = c("NA", "NaN", "", ".")
)

if (anyDuplicated(rownames(expr_df))) stop("Duplicated gene IDs were found")
if (anyDuplicated(colnames(expr_df))) stop("Duplicated sample IDs were found")

raw_matrix <- as.matrix(expr_df)
raw_na <- sum(is.na(raw_matrix))
suppressWarnings(storage.mode(raw_matrix) <- "numeric")
if (sum(is.na(raw_matrix)) > raw_na) {
  stop("Expression matrix contains non-numeric values", call. = FALSE)
}
if (anyNA(raw_matrix)) {
  stop("Expression matrix contains NA values; handle them before ICA", call. = FALSE)
}
if (any(!is.finite(raw_matrix))) stop("Expression matrix contains Inf values")
if (any(raw_matrix < 0)) stop("Expression values must be non-negative")
if (ncol(raw_matrix) < 3L) stop("At least three samples are required")

transformed_matrix <- log2(raw_matrix + log_offset)
gene_sd <- apply(transformed_matrix, 1, sd)
keep_gene <- is.finite(gene_sd) & gene_sd > 0

if (any(!keep_gene)) {
  excluded <- data.frame(
    Gene = rownames(transformed_matrix)[!keep_gene],
    Reason = "zero_or_invalid_variance"
  )
  write_tsv(excluded, file.path(output_dir, "excluded_before_ICA.tsv"))
  transformed_matrix <- transformed_matrix[keep_gene, , drop = FALSE]
}

cat("Genes used:", nrow(transformed_matrix), "\n")
cat("Samples:", ncol(transformed_matrix), "\n")
cat("Transform: log2(x +", log_offset, ")\n")
cat("Running ICA ...\n")

set.seed(seed)
ica_result <- runICA(
  transformed_matrix,
  var_cutoff = var_cutoff,
  max_iter = max_iter,
  n_runs = n_runs,
  n_cores = n_cores
)

gene_weights <- as.matrix(ica_result$S)
component_activity <- as.matrix(ica_result$A)

if (nrow(gene_weights) != nrow(transformed_matrix)) {
  stop("Unexpected S dimensions: rows do not match genes")
}
if (ncol(gene_weights) != nrow(component_activity)) {
  stop("Unexpected ICA dimensions: ncol(S) != nrow(A)")
}
if (ncol(component_activity) != ncol(transformed_matrix)) {
  stop("Unexpected A dimensions: columns do not match samples")
}

n_components <- ncol(gene_weights)
component_names <- paste0("IC", seq_len(n_components))
rownames(gene_weights) <- rownames(transformed_matrix)
colnames(gene_weights) <- component_names
rownames(component_activity) <- component_names
colnames(component_activity) <- colnames(transformed_matrix)

# Make component polarity reproducible: the largest absolute loading is positive.
# Flip S and the corresponding row of A together so reconstruction is unchanged.
orientation_anchor <- character(n_components)
for (j in seq_len(n_components)) {
  anchor <- which.max(abs(gene_weights[, j]))
  orientation_anchor[j] <- rownames(gene_weights)[anchor]
  if (gene_weights[anchor, j] < 0) {
    gene_weights[, j] <- -gene_weights[, j]
    component_activity[j, ] <- -component_activity[j, ]
  }
}

# Robustly standardize each component so loading magnitudes are comparable.
loading_z <- apply(gene_weights, 2, robust_z)
loading_z <- as.matrix(loading_z)
rownames(loading_z) <- rownames(gene_weights)
colnames(loading_z) <- colnames(gene_weights)

absolute_z <- abs(loading_z)
best_component_index <- max.col(absolute_z, ties.method = "first")
best_z <- loading_z[cbind(seq_len(nrow(loading_z)), best_component_index)]
best_abs_z <- abs(best_z)

if (n_components > 1L) {
  second_abs_z <- apply(absolute_z, 1, function(v) {
    sort(v, decreasing = TRUE)[2L]
  })
} else {
  second_abs_z <- rep(0, nrow(loading_z))
}

delta_z <- best_abs_z - second_abs_z
assigned_by_loading <- best_abs_z >= loading_z_cutoff & delta_z >= min_delta_z
polarity <- ifelse(best_z >= 0, "pos", "neg")
initial_module <- paste0(
  component_names[best_component_index], "_", polarity
)
initial_module[!assigned_by_loading] <- "Unassigned"

assignment_reason <- ifelse(
  best_abs_z < loading_z_cutoff,
  "weak_loading",
  ifelse(delta_z < min_delta_z, "ambiguous_loading", "assigned")
)

quality <- ifelse(
  !assigned_by_loading,
  "Unassigned",
  ifelse(
    best_abs_z >= max(3, loading_z_cutoff + 1) &
      delta_z >= max(1, 2 * min_delta_z),
    "High_confidence",
    ifelse(
      best_abs_z >= max(2.5, loading_z_cutoff + 0.5) &
        delta_z >= max(0.5, min_delta_z),
      "Medium_confidence",
      "Low_confidence"
    )
  )
)

selected_weight <- gene_weights[
  cbind(seq_len(nrow(gene_weights)), best_component_index)
]

assignment <- data.frame(
  Gene = rownames(gene_weights),
  Best_Component = component_names[best_component_index],
  Direction = ifelse(best_z >= 0, "positive", "negative"),
  ICA_Weight = selected_weight,
  Loading_Z = best_z,
  Abs_Loading_Z = best_abs_z,
  Second_Abs_Loading_Z = second_abs_z,
  Delta_Z = delta_z,
  Initial_Module = initial_module,
  Final_Module = initial_module,
  Assignment_Quality = quality,
  Assignment_Reason = assignment_reason,
  kME = NA_real_,
  stringsAsFactors = FALSE
)

# Refine each signed ICA tail using expression coherence. Genes whose signed
# module membership is below min_kME are removed, then PC1 is recalculated.
cat("Refining signed modules using kME ...\n")
max_refinement_rounds <- 10L
for (round_i in seq_len(max_refinement_rounds)) {
  changed <- FALSE
  module_names <- setdiff(unique(assignment$Final_Module), "Unassigned")

  for (module_name in module_names) {
    idx <- which(assignment$Final_Module == module_name)

    if (length(idx) < min_module_size) {
      assignment$Final_Module[idx] <- "Unassigned"
      assignment$Assignment_Reason[idx] <- "small_module"
      assignment$Assignment_Quality[idx] <- "Unassigned"
      changed <- TRUE
      next
    }

    genes <- assignment$Gene[idx]
    eig <- calculate_eigengene(transformed_matrix[genes, , drop = FALSE])
    assignment$kME[idx] <- eig$kME[genes]
    remove_idx <- idx[!is.finite(assignment$kME[idx]) |
                        assignment$kME[idx] < min_kME]

    if (length(remove_idx) > 0L) {
      assignment$Final_Module[remove_idx] <- "Unassigned"
      assignment$Assignment_Reason[remove_idx] <- "low_signed_kME"
      assignment$Assignment_Quality[remove_idx] <- "Unassigned"
      assignment$kME[remove_idx] <- NA_real_
      changed <- TRUE
    }
  }

  cat("  refinement round", round_i,
      "assigned genes:", sum(assignment$Final_Module != "Unassigned"), "\n")
  if (!changed) break
}

# Final eigengenes and kME after all removals.
final_modules <- sort(setdiff(unique(assignment$Final_Module), "Unassigned"))
eigengene_matrix <- matrix(
  NA_real_, nrow = ncol(transformed_matrix), ncol = length(final_modules),
  dimnames = list(colnames(transformed_matrix), final_modules)
)

for (module_name in final_modules) {
  idx <- which(assignment$Final_Module == module_name)
  genes <- assignment$Gene[idx]
  eig <- calculate_eigengene(transformed_matrix[genes, , drop = FALSE])
  eigengene_matrix[, module_name] <- eig$score[colnames(transformed_matrix)]
  assignment$kME[idx] <- eig$kME[genes]
}

assignment <- assignment[order(
  assignment$Final_Module == "Unassigned",
  assignment$Final_Module,
  -assignment$Abs_Loading_Z,
  assignment$Gene
), , drop = FALSE]

cat("Writing output files ...\n")
write_tsv(
  assignment,
  file.path(output_dir, "ICA_signed_gene_assignment.tsv")
)
write_tsv(
  assignment[assignment$Final_Module == "Unassigned", , drop = FALSE],
  file.path(output_dir, "ICA_unassigned_genes.tsv")
)
write_tsv(
  data.frame(Gene = rownames(gene_weights), gene_weights, check.names = FALSE),
  file.path(output_dir, "ICA_oriented_gene_weights.tsv")
)
write_tsv(
  data.frame(Gene = rownames(loading_z), loading_z, check.names = FALSE),
  file.path(output_dir, "ICA_robust_loading_zscores.tsv")
)
write_tsv(
  data.frame(
    Sample = colnames(component_activity),
    t(component_activity), check.names = FALSE
  ),
  file.path(output_dir, "ICA_component_activities.samples_by_components.tsv")
)
write_tsv(
  data.frame(
    Sample = rownames(eigengene_matrix),
    eigengene_matrix, check.names = FALSE
  ),
  file.path(output_dir, "ICA_signed_module_eigengenes.samples_by_modules.tsv")
)

orientation_table <- data.frame(
  Component = component_names,
  Positive_Anchor_Gene = orientation_anchor,
  stringsAsFactors = FALSE
)
write_tsv(
  orientation_table,
  file.path(output_dir, "ICA_component_orientation.tsv")
)

if (!is.null(ica_result$percent_var)) {
  percent_var <- as.numeric(ica_result$percent_var)
  write_tsv(
    data.frame(
      Component = component_names[seq_along(percent_var)],
      Percent_Variance = percent_var
    ),
    file.path(output_dir, "ICA_component_percent_variance.tsv")
  )
}

# Per-module gene files. Weight remains signed; no absolute value is substituted.
for (module_name in final_modules) {
  module_df <- assignment[
    assignment$Final_Module == module_name,
    c(
      "Gene", "ICA_Weight", "Loading_Z", "Abs_Loading_Z",
      "Delta_Z", "kME", "Assignment_Quality"
    ),
    drop = FALSE
  ]
  module_df <- module_df[order(-module_df$Abs_Loading_Z), , drop = FALSE]

  write_tsv(
    module_df,
    file.path(output_dir, paste0("Module_", module_name, "_genes.tsv"))
  )

  simple_df <- data.frame(
    Gene = module_df$Gene,
    Weight = module_df$ICA_Weight,
    Quality = module_df$Assignment_Quality,
    Loading_Z = module_df$Loading_Z,
    kME = module_df$kME
  )
  write_tsv(
    simple_df,
    file.path(output_dir, paste0("Module_", module_name, "_genes_simple.txt"))
  )
}

if (length(final_modules) > 0L) {
  module_summary <- do.call(rbind, lapply(final_modules, function(m) {
    z <- assignment[assignment$Final_Module == m, , drop = FALSE]
    data.frame(
      Module = m,
      Component = z$Best_Component[1],
      Direction = z$Direction[1],
      Genes = nrow(z),
      Median_Abs_Loading_Z = median(z$Abs_Loading_Z, na.rm = TRUE),
      Median_kME = median(z$kME, na.rm = TRUE),
      Min_kME = min(z$kME, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
} else {
  module_summary <- data.frame(
    Module = character(), Component = character(), Direction = character(),
    Genes = integer(), Median_Abs_Loading_Z = numeric(),
    Median_kME = numeric(), Min_kME = numeric()
  )
}
write_tsv(module_summary, file.path(output_dir, "ICA_signed_module_summary.tsv"))

settings <- data.frame(
  Parameter = c(
    "expression_file", "genes_input", "genes_used", "samples",
    "transform", "var_cutoff", "n_runs", "max_iter", "n_cores", "seed",
    "loading_standardization", "loading_z_cutoff", "min_delta_z",
    "min_kME", "min_module_size", "components", "signed_modules",
    "assigned_genes", "unassigned_genes"
  ),
  Value = c(
    normalizePath(expression_file), nrow(raw_matrix), nrow(transformed_matrix),
    ncol(transformed_matrix), paste0("log2(x+", log_offset, ")"),
    var_cutoff, n_runs, max_iter, n_cores, seed, "median/MAD", loading_z_cutoff,
    min_delta_z, min_kME, min_module_size, n_components,
    length(final_modules), sum(assignment$Final_Module != "Unassigned"),
    sum(assignment$Final_Module == "Unassigned")
  )
)
write_tsv(settings, file.path(output_dir, "ICA_signed_analysis_settings.tsv"))

# Compact diagnostic PDF using ASCII text to avoid encoding warnings.
pdf(
  file.path(output_dir, "ICA_signed_diagnostics.pdf"),
  width = 14, height = 9, useDingbats = FALSE
)
par(mfrow = c(2, 3), mar = c(5, 5, 3, 1))

barplot(
  table(assignment$Assignment_Reason), las = 2, col = "#4DBBD5",
  main = "Assignment outcomes", ylab = "Genes"
)
hist(
  assignment$Abs_Loading_Z, breaks = 60, col = "#00A087",
  border = "white", main = "Best absolute loading Z", xlab = "Absolute robust Z"
)
abline(v = loading_z_cutoff, col = "#E64B35", lwd = 2, lty = 2)
hist(
  assignment$Delta_Z, breaks = 60, col = "#3C5488",
  border = "white", main = "Loading dominance", xlab = "Best Z - second-best Z"
)
abline(v = min_delta_z, col = "#E64B35", lwd = 2, lty = 2)

if (nrow(module_summary) > 0L) {
  module_colors <- ifelse(module_summary$Direction == "positive", "#E64B35", "#3C5488")
  barplot(
    module_summary$Genes, names.arg = module_summary$Module,
    las = 2, cex.names = 0.65, col = module_colors,
    main = "Final signed module sizes", ylab = "Genes"
  )
  plot(
    module_summary$Median_Abs_Loading_Z, module_summary$Median_kME,
    pch = 21, bg = module_colors, col = "#333333", cex = 1.3,
    xlab = "Median absolute loading Z", ylab = "Median kME",
    main = "Module loading and coherence"
  )
} else {
  plot.new(); title("No final signed modules")
  plot.new()
}

heatmap(
  component_activity,
  scale = "row",
  Colv = NA,
  labCol = NA,
  col = colorRampPalette(c("#3C5488", "white", "#E64B35"))(101),
  main = "ICA component activities"
)
dev.off()

sink(file.path(output_dir, "ICA_signed_report.txt"))
cat("Signed two-tailed ICA report\n")
cat("=============================\n")
cat("Genes supplied:", nrow(raw_matrix), "\n")
cat("Genes used:", nrow(transformed_matrix), "\n")
cat("Samples:", ncol(transformed_matrix), "\n")
cat("ICA components:", n_components, "\n")
cat("Final signed modules:", length(final_modules), "\n")
cat("Assigned genes:", sum(assignment$Final_Module != "Unassigned"), "\n")
cat("Unassigned genes:", sum(assignment$Final_Module == "Unassigned"), "\n")
cat("Loading criterion: absolute robust Z >=", loading_z_cutoff, "\n")
cat("Dominance criterion: delta Z >=", min_delta_z, "\n")
cat("Expression coherence criterion: kME >=", min_kME, "\n")
cat("Minimum module size:", min_module_size, "\n\n")
cat("Important interpretation:\n")
cat("1. IC*_pos and IC*_neg are opposite tails of one ICA component.\n")
cat("2. Component activities are tested once per IC for trait association.\n")
cat("3. Signed module PC1 eigengenes can be used for arm-specific analyses.\n")
cat("4. Unassigned genes are retained in the assignment table, not discarded.\n")
sink()

cat("\nAnalysis completed\n")
cat("Output directory:", normalizePath(output_dir), "\n")
cat("ICA components:", n_components, "\n")
cat("Final signed modules:", length(final_modules), "\n")
cat("Assigned genes:", sum(assignment$Final_Module != "Unassigned"), "\n")
cat("Unassigned genes:", sum(assignment$Final_Module == "Unassigned"), "\n")
