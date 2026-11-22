library(DESeq2)
library(ggplot2)

dds <- readRDS(snakemake@input[["dds"]])
design_cols <- all.vars(design(dds))
group_col <- snakemake@params[["group_column"]]
meta <- as.data.frame(colData(dds))[, design_cols, drop = FALSE]

rld <- rlog(dds, blind = TRUE)
pca_data <- plotPCA(rld, intgroup = design_cols, returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)
# 从元数据按样本名取值，也适用于分组列名恰好叫 group 的设计。
pca_data$condition <- meta[pca_data$name, group_col]
p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = condition))
shape_cols <- setdiff(design_cols, group_col)
if (length(shape_cols) && is.factor(meta[[shape_cols[1]]]) &&
    nlevels(meta[[shape_cols[1]]]) <= 6L) {
  pca_data$block <- meta[pca_data$name, shape_cols[1]]
  p_pca <- p_pca + geom_point(data = pca_data, aes(shape = block), size = 3) +
    labs(shape = shape_cols[1])
} else {
  p_pca <- p_pca + geom_point(size = 3)
}
p_pca <- p_pca + geom_text(aes(label = name), vjust = -0.9, size = 3, check_overlap = TRUE) +
  labs(x = sprintf("PC1 (%.1f%%)", percent_var[1]), y = sprintf("PC2 (%.1f%%)", percent_var[2]),
       colour = group_col) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) + coord_fixed() + theme_bw()
ggsave(snakemake@output[["pca"]], p_pca, width = 7, height = 6, dpi = 300)

rld_cor <- cor(assay(rld), method = "pearson")
distance <- as.dist(pmax(1 - rld_cor, 0))
pheatmap::pheatmap(rld_cor, annotation_col = meta, clustering_distance_rows = distance,
                   clustering_distance_cols = distance, border_color = NA,
                   filename = snakemake@output[["correlation_heatmap"]], width = 8, height = 7)
