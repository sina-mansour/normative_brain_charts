# ============================================================================
# SigClust2 hierarchical significance testing (Kimes et al. 2017)
# doi:10.1111/biom.12647
#
# Runs the SHC procedure across representations of the AD deviation maps.
# Writes per-node test outcomes, root-node null distributions, cluster
# assignments, a cross-configuration summary, and a sessionInfo record.
#
# Usage: Rscript sigclust2_analysis.R [data_dir]
# ============================================================================

library(arrow)
library(sigclust2)

args <- commandArgs(trailingOnly = TRUE)
sigclust_dir <- if (length(args) > 0) args[1] else
  "/mountpoint/code/projects/normative_brain_charts/data/sigclust"

# --- Testing parameters -----------------------------------------------------
KS      <- 2:10     # cluster counts exported for visualisation
N_SIM   <- 10000    # simulated null datasets per node
ALPHA   <- 0.05     # < 1 enables the sequential FWER-controlling procedure
ICOVEST <- 1        # soft thresholding (recommended, Huang et al. 2014)

# --- Analysis configurations ------------------------------------------------
# null_alg: "hclust" gives greater power for rotation-invariant metric/linkage
#           combinations; "2means" is recommended only for non-rotation-invariant
#           metrics and yields conservative p-values with ci = "2CI"
configs <- list(
  list(name = "mds", file = "sigclust2_input.parquet", col_re = "^MDS",
       null_alg = "hclust",
       label = "2D MDS landscape, hclust null"),
  list(name = "pca", file = "sigclust2_pca_input.parquet", col_re = "^PC",
       null_alg = "hclust",
       label = "20-component PCA, hclust null"),
  # Conservative variant, reported as a robustness check
  list(name = "mds_2means", file = "sigclust2_input.parquet", col_re = "^MDS",
       null_alg = "2means",
       label = "2D MDS landscape, 2-means null"),
  list(name = "pca_2means", file = "sigclust2_pca_input.parquet", col_re = "^PC",
       null_alg = "2means",
       label = "20-component PCA, 2-means null")
)

summary_rows <- list()

for (cfg in configs) {

  cat("\n\n==========================================================\n")
  cat("SigClust2:", cfg$label, "\n")
  cat("==========================================================\n")

  # --- Load ---------------------------------------------------------------
  df <- read_parquet(file.path(sigclust_dir, cfg$file))
  cols <- grep(cfg$col_re, names(df), value = TRUE)
  X <- as.matrix(df[, cols])
  cat("Data:", nrow(X), "observations x", ncol(X), "features\n")
  cat("Settings: icovest =", ICOVEST,
      "| null_alg =", cfg$null_alg,
      "| n_sim =", N_SIM, "| alpha =", ALPHA, "\n")

  # --- Run SHC ------------------------------------------------------------
  set.seed(42)
  t0 <- Sys.time()
  shc_res <- shc(
    X,
    metric   = "euclidean",
    linkage  = "ward.D2",     # equivalent to scipy's method = "ward"
    alpha    = ALPHA,
    n_sim    = N_SIM,
    icovest  = ICOVEST,
    ci       = "2CI",
    null_alg = cfg$null_alg,
    rcpp     = FALSE          # Rclusterpp is unmaintained
  )
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  cat("Elapsed:", round(as.numeric(elapsed), 2), "minutes\n")

  print(summary(shc_res))

  # --- Extract node-level results -----------------------------------------
  # NOTE: p_norm and p_emp are ordered DESCENDING from the top of the
  # dendrogram, so entry 1 is the root. This is the reverse of hclust's
  # merge-matrix ordering, and the two must not be bound together.
  # p_norm and p_emp are alternative p-value estimates (Gaussian
  # approximation and empirical), not raw and corrected values; p_norm is
  # the statistic tested against the FWER cutoffs by default.
  p_norm_vec <- as.vector(shc_res$p_norm)
  p_emp_vec  <- as.vector(shc_res$p_emp)
  nd_type    <- shc_res$nd_type

  node_results <- data.frame(
    node_rank = seq_along(nd_type),   # 1 = root, increasing down the tree
    nd_type   = nd_type,
    p_norm    = p_norm_vec,
    p_emp     = p_emp_vec
  )

  tested <- nd_type %in% c("sig", "not_sig")

  cat("\n--- Node status ---\n")
  print(table(nd_type))
  cat("Nodes tested      :", sum(tested), "/", length(nd_type), "\n")
  cat("Significant nodes :", sum(nd_type == "sig"), "\n")
  cat("Root node         : nd_type =", nd_type[1],
      "| p_norm =", round(p_norm_vec[1], 4),
      "| p_emp =", round(p_emp_vec[1], 4), "\n")

  # --- Root-node cluster index and its null distribution ------------------
  stopifnot(nrow(shc_res$ci_dat) == length(nd_type))
  ci_obs  <- shc_res$ci_dat[1, 1]
  ci_null <- shc_res$ci_sim[1, , 1]
  stopifnot(length(ci_null) == N_SIM)
  cat("Root cluster index:", round(ci_obs, 5),
      "| null mean =", round(mean(ci_null), 5),
      "| null sd =", round(sd(ci_null), 5), "\n")

  summary_rows[[cfg$name]] <- data.frame(
    config       = cfg$name,
    label        = cfg$label,
    n_features   = ncol(X),
    icovest      = ICOVEST,
    null_alg     = cfg$null_alg,
    ci_obs       = ci_obs,
    ci_null_mean = mean(ci_null),
    ci_null_sd   = sd(ci_null),
    ci_null_q05  = as.numeric(quantile(ci_null, 0.05)),
    p_norm       = p_norm_vec[1],
    p_emp        = p_emp_vec[1],
    nd_type      = nd_type[1],
    n_tested     = sum(tested),
    n_sig        = sum(nd_type == "sig"),
    stringsAsFactors = FALSE
  )

  # --- Cluster assignments for visualisation ------------------------------
  assign_mat <- sapply(KS, function(k) cutree(shc_res$hc_dat, k = k))
  colnames(assign_mat) <- paste0("k", KS)
  cluster_df <- data.frame(subject_idx = seq_len(nrow(X)) - 1,  # 0-indexed
                           assign_mat)

  # --- Write --------------------------------------------------------------
  res_path <- file.path(sigclust_dir, sprintf("sigclust2_%s_results.csv", cfg$name))
  clu_path <- file.path(sigclust_dir, sprintf("sigclust2_%s_clusters.csv", cfg$name))
  cin_path <- file.path(sigclust_dir, sprintf("sigclust2_%s_ci_null.csv", cfg$name))
  write.csv(node_results, res_path, row.names = FALSE)
  write.csv(cluster_df,   clu_path, row.names = FALSE)
  write.csv(data.frame(ci_null = ci_null), cin_path, row.names = FALSE)
  cat("\nWrote:", res_path, "\n      ", clu_path, "\n      ", cin_path, "\n")
}

# --- Cross-configuration summary --------------------------------------------
summary_df <- do.call(rbind, summary_rows)
sum_path <- file.path(sigclust_dir, "sigclust2_summary.csv")
write.csv(summary_df, sum_path, row.names = FALSE)

cat("\n\n=== Summary across configurations ===\n")
print(summary_df[, c("config", "n_features", "null_alg",
                     "ci_obs", "p_norm", "p_emp", "nd_type")],
      row.names = FALSE)
cat("\nNote: ci_obs depends only on the data and dendrogram, so it is identical\n")
cat("across null algorithms within each representation.\n")
cat("\nWrote summary to", sum_path, "\n")

# --- Record environment -----------------------------------------------------
si_path <- file.path(sigclust_dir, "R_sessionInfo.txt")
writeLines(capture.output(sessionInfo()), si_path)
cat("Wrote sessionInfo to", si_path, "\n")
