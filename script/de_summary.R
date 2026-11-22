library(DESeq2)
library(ggplot2)

res <- read.delim(snakemake@input[["res"]], colClasses = c(gene_id = "character"),
                  check.names = FALSE)
dds <- readRDS(snakemake@input[["dds"]])
group_col <- snakemake@params[["group_column"]]
numerator <- snakemake@params[["numerator"]]
denominator <- snakemake@params[["denominator"]]
alpha <- snakemake@params[["alpha"]]
lfc_cutoff <- log2(snakemake@params[["min_fold_change"]])

# 差异基因表、火山图、热图及后续 ORA 共用这一次判定。
eligible <- is.finite(res$padj) & is.finite(res$log2FC_shrunk)
sig <- eligible & res$padj < alpha & abs(res$log2FC_shrunk) >= lfc_cutoff
res$direction <- "NS"
res$direction[sig & res$log2FC_shrunk > 0] <- "Up"
res$direction[sig & res$log2FC_shrunk < 0] <- "Down"
deg <- res[res$direction != "NS", ]
deg <- deg[order(deg$padj, -abs(deg$log2FC_shrunk), deg$gene_id), ]
readr::write_tsv(deg, snakemake@output[["deg"]])

plot_data <- res[eligible, ]
p <- ggplot(plot_data, aes(log2FC_shrunk, -log10(pmax(padj, .Machine$double.xmin)),
                           colour = direction)) +
  geom_vline(xintercept = unique(c(-lfc_cutoff, lfc_cutoff)), linetype = 2, colour = "grey55") +
  geom_hline(yintercept = -log10(alpha), linetype = 2, colour = "grey55") +
  geom_point(alpha = 0.65, size = 1.4) +
  scale_colour_manual(values = c(Up = "#BF423D", Down = "#336EA1", NS = "grey75"),
                      breaks = c("Up", "Down", "NS")) +
  labs(title = paste(numerator, "vs", denominator), x = "Shrunken log2 fold change",
       y = expression(-log[10]("adjusted P")), colour = NULL) + theme_bw()
if (nrow(deg)) {
  p <- p + ggrepel::geom_text_repel(data = head(deg, 10), aes(label = gene_id),
                                    size = 3, seed = 123, show.legend = FALSE)
}
ggsave(snakemake@output[["volcano"]], p, width = 7, height = 6, dpi = 300)

# 热图只展示当前比较的样本，列注释与表达矩阵取自同一 DDS。
keep_samples <- dds[[group_col]] %in% c(numerator, denominator)
anno <- as.data.frame(colData(dds))[keep_samples, all.vars(design(dds)), drop = FALSE]
mat <- log2(counts(dds, normalized = TRUE)[deg$gene_id, rownames(anno), drop = FALSE] + 1)
if (nrow(mat)) {
  mat <- mat[matrixStats::rowVars(mat) > 0, , drop = FALSE]
}
heatmap_file <- snakemake@output[["heatmap"]]
if (nrow(mat)) {
  pheatmap::pheatmap(mat, scale = "row", annotation_col = anno, cluster_rows = nrow(mat) > 1,
                     color = colorRampPalette(c("#366A9F", "white", "#B9413D"))(101),
                     show_rownames = FALSE, border_color = NA, filename = heatmap_file,
                     width = 8, height = 8)
} else {
  label <- if (nrow(deg)) "No DEGs with variable expression" else "No DEGs at the selected thresholds"
  p_empty <- ggplot() + annotate("text", x = 0, y = 0, label = label) + theme_void()
  ggsave(heatmap_file, p_empty, width = 8, height = 8, dpi = 300)
}
