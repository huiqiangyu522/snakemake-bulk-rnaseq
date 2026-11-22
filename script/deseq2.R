library(DESeq2)
# 读取
alpha <- snakemake@params[["alpha"]]
dds <- readRDS(snakemake@input[["dds"]])


group_col <- snakemake@params[["group_column"]]
numerator <- snakemake@params[["numerator"]]
denominator <- snakemake@params[["denominator"]]
groups <- dds[[group_col]]
if (!all(c(numerator, denominator) %in% levels(groups))) {
  stop("无效比较：", numerator, " vs ", denominator,
       "。样本表中的分组为：", paste(levels(groups), collapse = ", "))
}
dds[[group_col]] <- relevel(groups, ref = denominator)

# 运行DEseq
dds <- DESeq(dds, test = "Wald")

# 获取对比表
coef_name <- make.names(paste0(group_col, "_", numerator, "_vs_", denominator))
if (!coef_name %in% resultsNames(dds)) {
  stop("未找到系数 ", coef_name, "。模型中的系数为：",
       paste(resultsNames(dds), collapse = ", "))
}

res_raw <- results(dds, name = coef_name, alpha = alpha)
res_shrunk <- lfcShrink(dds, coef = coef_name, res = res_raw, type = "apeglm")

stopifnot(identical(rownames(res_raw), rownames(res_shrunk)))

res_tbl <- data.frame(
  gene_id = rownames(res_raw),
  baseMean = res_raw$baseMean,
  log2FC_raw = res_raw$log2FoldChange,
  lfcSE_raw = res_raw$lfcSE,
  stat = res_raw$stat,
  pvalue = res_raw$pvalue,
  padj = res_raw$padj,
  log2FC_shrunk = res_shrunk$log2FoldChange,
  lfcSD_shrunk = res_shrunk$lfcSE
)


write.table(
  res_tbl,
  file = snakemake@output[["results"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
