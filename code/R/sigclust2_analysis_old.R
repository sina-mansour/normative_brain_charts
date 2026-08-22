# ============================================================================
# SigClust2 hierarchical significance testing
# Kimes et al. 2017, doi:10.1111/biom.12647
# ============================================================================

library(arrow)
library(sigclust2)

# --- Load data ---
data_path <- "/mountpoint/code/projects/normative_brain_charts/data/sigclust/sigclust2_input.parquet"
df <- read_parquet(data_path)
X <- as.matrix(df[, c("MDS1", "MDS2")])
cat("Loaded data with", nrow(X), "subjects and", ncol(X), "features.\n")

# --- Run shc ---
set.seed(42)
shc_result <- shc(
  X,
  metric    = "euclidean",
  linkage   = "ward.D2",
  alpha     = 0.05,
  n_sim     = 1000,
  icovest   = 1,
  ci        = "2CI",
  null_alg  = "2means",
  rcpp      = FALSE
)

cat("\n=== SigClust2 summary ===\n")
print(shc_result)

# --- Extract per-node results as a clean data frame ---
# p_norm and p_emp are matrices; flatten to vectors (we used a single test statistic, ci="2CI")
p_raw       <- as.vector(shc_result$p_norm)
p_corrected <- as.vector(shc_result$p_emp)

# The hclust merge matrix tells us which subjects/clusters merge at each node.
# Rows = nodes (1 to n-1), values < 0 are leaf indices (subjects), values > 0 are internal nodes.
merge_matrix <- shc_result$hc_dat$merge
heights      <- shc_result$hc_dat$height

node_results <- data.frame(
  node_id     = seq_len(nrow(merge_matrix)),
  merge_left  = merge_matrix[, 1],
  merge_right = merge_matrix[, 2],
  height      = heights,
  p_raw       = p_raw,
  p_corrected = p_corrected,
  tested      = !is.na(p_raw) & p_raw <= 1,    # p > 1 sentinel = not tested
  significant = !is.na(p_corrected) & p_corrected < 0.05
)

cat("\n=== Per-node results ===\n")
print(head(node_results, 20))
cat("...\n")
cat("Nodes actually tested:", sum(node_results$tested), "\n")
cat("Significant nodes (FWER < 0.05):", sum(node_results$significant), "\n")
cat("Top-level (root) split p-value:", node_results$p_raw[nrow(node_results)], "\n")

# --- Save outputs for Python visualisation ---
out_dir <- "/mountpoint/code/projects/normative_brain_charts/data/sigclust"

# Per-node results
write.csv(node_results,
          file = file.path(out_dir, "sigclust2_node_results.csv"),
          row.names = FALSE)

# Cluster assignments at every k from 2 to 10, so Python can colour-overlay any clustering
ks <- 2:10
cluster_assignments <- sapply(ks, function(k) cutree(shc_result$hc_dat, k = k))
colnames(cluster_assignments) <- paste0("k", ks)
cluster_df <- as.data.frame(cluster_assignments)
cluster_df$subject_idx <- seq_len(nrow(cluster_df)) - 1   # 0-indexed for Python compatibility
cluster_df <- cluster_df[, c("subject_idx", colnames(cluster_assignments))]

write.csv(cluster_df,
          file = file.path(out_dir, "sigclust2_cluster_assignments.csv"),
          row.names = FALSE)

cat("\nResults written to:\n")
cat("  ", file.path(out_dir, "sigclust2_node_results.csv"), "\n")
cat("  ", file.path(out_dir, "sigclust2_cluster_assignments.csv"), "\n")