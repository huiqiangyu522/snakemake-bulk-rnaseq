library(DESeq2)

count_table <- read.delim(snakemake@input[["counts"]], row.names = "gene_id",
                          colClasses = c(gene_id = "character"), check.names = FALSE)
samples <- read.delim(snakemake@input[["samples"]], row.names = "sample_id",
                      colClasses = "character", check.names = FALSE)

design_formula <- as.formula(snakemake@params[["design"]])
design_cols <- all.vars(design_formula)
group_col <- snakemake@params[["group_column"]]
if (!group_col %in% design_cols) {
  stop("group_column 必须出现在设计式中。")
}
# 按分组的单一系数进行 apeglm 收缩，分组须作为有截距模型中的加性主效应。
design_terms <- terms(design_formula)
group_terms <- attr(design_terms, "factors")
if (attr(design_terms, "intercept") != 1L || !group_col %in% rownames(group_terms) ||
    any(attr(design_terms, "order")[group_terms[group_col, ] > 0] > 1L)) {
  stop("当前比较要求设计式含截距，且分组以无交互的主效应进入模型。")
}
if (!setequal(colnames(count_table), rownames(samples))) {
  stop("计数矩阵中的样本与样本表不一致。")
}
meta <- samples[colnames(count_table), design_cols, drop = FALSE]
for (col in design_cols) {
  if (col %in% c(group_col, "donor")) {
    meta[[col]] <- factor(meta[[col]])
  } else {
    meta[[col]] <- type.convert(meta[[col]], as.is = FALSE)
  }
}

# 所有比较共享计数、设计和归一化因子；参考组在每项比较拟合前设置。
dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_table), colData = meta,
                              design = design_formula)
dds <- estimateSizeFactors(dds)
saveRDS(dds, snakemake@output[["dds"]])
