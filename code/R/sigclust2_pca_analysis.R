# ============================================================================
# SigClust2 on PCA-reduced deviation maps
# ============================================================================

library(arrow)
library(sigclust2)

sigclust_dir <- "/mountpoint/code/projects/normative_brain_charts/data/sigclust"

# --- Load data ---
df <- read_parquet(file.path(sigclust_dir, "sigclust2_pca_input.parquet"))
# Drop the subject_idx column if it was included; keep only the PC columns
pc_cols <- grep("^PC", names(df), value = TRUE)
X <- as.matrix(df[, pc_cols])
cat("Loaded data with", nrow(X), "subjects and", ncol(X), "PCs.\n")

# --- Run shc ---
# Key change vs MDS version: icovest = 2 (soft-thresholded covariance).
# Kimes recommends this for high-dim data where sample covariance is ill-conditioned.
# With p = 10 and n ~ 200 we're still moderately well-conditioned, but icovest=2
# is the safer default once p > 5 or so.
set.seed(42)
shc_result <- shc(
  X,
  metric    = "euclidean",
  linkage   = "ward.D2",
  alpha     = 0.05,
  n_sim     = 1000,
  icovest   = 2,            # soft-thresholded covariance for higher-dim data
  ci        = "2CI",
  null_alg  = "2means",
  rcpp      = FALSE
)

cat("\n=== SigClust2 summary ===\n")
print(shc_result)

# --- Extract per-node results ---
p_raw       <- as.vector(shc_result$p_norm)
p_corrected <- as.vector(shc_result$p_emp)

merge_matrix <- shc_result$hc_dat$merge
heights      <- shc_result$hc_dat$height

node_results <- data.frame(
  node_id     = seq_len(nrow(merge_matrix)),
  merge_left  = merge_matrix[, 1],
  merge_right = merge_matrix[, 2],
  height      = heights,
  p_raw       = p_raw,
  p_corrected = p_corrected,
  tested      = !is.na(p_raw) & p_raw <= 1,
  significant = !is.na(p_corrected) & p_corrected < 0.05
)

cat("\n=== Per-node results (first 20) ===\n")
print(head(node_results, 20))
cat("\nTotal tested:", sum(node_results$tested), "\n")
cat("Significant (FWER < 0.05):", sum(node_results$significant), "\n")

# --- Save outputs ---
write.csv(node_results,
          file = file.path(sigclust_dir, "sigclust2_pca_node_results.csv"),
          row.names = FALSE)

# Cluster assignments at each k from 2 to 10
ks <- 2:10
cluster_assignments <- sapply(ks, function(k) cutree(shc_result$hc_dat, k = k))
colnames(cluster_assignments) <- paste0("k", ks)
cluster_df <- as.data.frame(cluster_assignments)
cluster_df$subject_idx <- seq_len(nrow(cluster_df)) - 1
cluster_df <- cluster_df[, c("subject_idx", colnames(cluster_assignments))]

write.csv(cluster_df,
          file = file.path(sigclust_dir, "sigclust2_pca_cluster_assignments.csv"),
          row.names = FALSE)

cat("\nResults written to:\n")
cat("  ", file.path(sigclust_dir, "sigclust2_pca_node_results.csv"), "\n")
cat("  ", file.path(sigclust_dir, "sigclust2_pca_cluster_assignments.csv"), "\n")
