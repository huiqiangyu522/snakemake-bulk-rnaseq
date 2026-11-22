library(htmltools)

p <- snakemake@params
inputs <- snakemake@input
output <- snakemake@output[["html"]]
comparisons <- names(p$contrasts)

# 从 Snakemake 声明的文件读取；每个向量的第 i 个元素属于同一项比较。
table_html <- function(tbl) {
  HTML(as.character(knitr::kable(tbl, format = "html", row.names = FALSE, escape = TRUE)))
}
figure <- function(path, caption) {
  tags$figure(tags$img(src = base64enc::dataURI(file = path, mime = "image/png"),
                       alt = caption, loading = "lazy"), tags$figcaption(caption))
}
download <- function(path, label) {
  tags$a(class = "download", download = basename(path),
         href = base64enc::dataURI(file = path, mime = "text/tab-separated-values"), label)
}
read_result <- function(key, i) {
  classes <- if (key %in% c("de", "deg")) c(gene_id = "character") else NA
  read.delim(inputs[[key]][i], check.names = FALSE, colClasses = classes)
}
format_p <- function(x) formatC(x, format = "g", digits = 3)

samples <- read.delim(inputs$samples, colClasses = "character", check.names = FALSE)
sample_cols <- unique(c("sample_id", all.vars(as.formula(p$design))))
qc <- tags$section(tags$h2("样本与质量控制"),
                   table_html(samples[, sample_cols, drop = FALSE]),
                   tags$div(class = "plots", figure(inputs$pca, "rlog 后的样本 PCA"),
                            figure(inputs$correlation_heatmap, "样本 Pearson 相关性与层次聚类")))

ora_names <- c(custom = "本地基因集", go_bp = "GO BP", kegg = "KEGG")
gsea_name <- if (length(p$term2gene)) "本地基因集 GSEA" else "Hallmark GSEA"
pathway_panel <- function(key, title, i) {
  path <- inputs[[paste0(key, "_table")]][i]
  tbl <- read_result(paste0(key, "_table"), i)
  sig <- tbl[which(tbl$significant), , drop = FALSE]
  text <- if (!nrow(tbl)) {
    "当前输入未返回可报告的条目。"
  } else {
    sprintf("共返回 %d 个条目，其中 %d 个达到通路 FDR 阈值。", nrow(tbl), nrow(sig))
  }
  preview <- NULL
  if (nrow(sig)) {
    measure <- if (key == "gsea") "NES" else "Count"
    top <- head(sig[, c("Description", measure, "p.adjust"), drop = FALSE], 5)
    top$p.adjust <- format_p(top$p.adjust)
    if (key == "gsea") top$NES <- round(top$NES, 2)
    preview <- table_html(top)
  }
  tags$details(tags$summary(title), tags$p(text), preview,
               figure(inputs[[paste0(key, "_plot")]][i], title),
               download(path, "下载完整富集表（TSV）"))
}

summary_rows <- list()
sections <- list()
for (i in seq_along(comparisons)) {
  cmp <- p$contrasts[[i]]
  title <- paste(cmp$numerator, "vs", cmp$denominator)
  de <- read_result("de", i)
  deg <- read_result("deg", i)
  eligible <- is.finite(de$padj) & is.finite(de$log2FC_shrunk)
  summary_rows[[i]] <- data.frame(
    比较 = title, 分子组样本 = sum(samples[[p$group_column]] == cmp$numerator),
    分母组样本 = sum(samples[[p$group_column]] == cmp$denominator),
    输入基因 = nrow(de), 可判定基因 = sum(eligible),
    上调 = sum(deg$direction == "Up"), 下调 = sum(deg$direction == "Down")
  )
  top <- head(deg[, c("gene_id", "log2FC_shrunk", "padj", "direction")], 10)
  top$log2FC_shrunk <- round(top$log2FC_shrunk, 3)
  top$padj <- format_p(top$padj)
  panels <- list()
  for (source in p$ora_sets) {
    for (direction in c("up", "down")) {
      label <- if (direction == "up") "上调基因" else "下调基因"
      key <- paste(source, direction, sep = "_")
      panels[[key]] <- pathway_panel(key, paste(ora_names[[source]], "ORA ·", label), i)
    }
  }
  panels$gsea <- pathway_panel("gsea", gsea_name, i)
  sections[[i]] <- tags$section(
    tags$h2(title), tags$p(paste0("正 log2FC 表示 ", cmp$numerator, " 组表达较高。")),
    tags$div(class = "plots", figure(inputs$volcano[i], "差异表达火山图"),
             figure(inputs$heatmap[i], "差异基因热图 · 归一化计数经 log2(x + 1) 转换，按行标准化")),
    tags$h3("差异基因"),
    if (nrow(top)) tagList(tags$p("按 padj 从小到大列出前 10 个，实际展示数量取决于显著基因数。"),
                           table_html(top)) else tags$p("当前阈值下未检出差异表达基因。"),
    tags$p(download(inputs$de[i], "下载 DESeq2 全表（TSV）"), "　",
           download(inputs$deg[i], "下载差异基因表（TSV）")),
    tags$h3("功能富集"), panels)
}

source_note <- if (length(p$term2gene)) {
  tags$p("ORA 和 GSEA 使用本地基因集 ", tags$code(p$term2gene), "。",
         if (basename(dirname(p$term2gene)) == "c0_fixture")
           "当前基因和通路为项目自带的模拟数据，结果用于验证工作流。")
} else {
  tags$p("GO BP 使用 Ensembl ID，KEGG 使用双向唯一映射的 Entrez ID。GSEA 使用对应物种的 Hallmark 基因集。")
}
page <- tags$html(lang = "zh-CN",
  tags$head(tags$meta(charset = "utf-8"), tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$title("Bulk RNA-seq 分析报告"),
    tags$style(HTML("
      * { box-sizing: border-box } body { margin: 0 auto; max-width: 1060px; padding: 2rem 5%;
      font: 16px/1.65 system-ui, sans-serif; color: #1d323b; background: #fff }
      h1,h2,h3 { line-height: 1.25 } h1 { font-size: 2.1rem }
      section { margin: 2.6rem 0 } h2 { border-bottom: 2px solid #367166; padding-bottom: .45rem }
      table { border-collapse: collapse; display: block; overflow-x: auto; margin: 1rem 0 }
      th,td { padding: .3rem .8rem; border-bottom: 1px solid #d7e1de } th { background: #edf4f1 }
      .plots { display: grid; grid-template-columns: repeat(auto-fit,minmax(min(100%,320px),1fr)); gap: 1rem }
      figure { min-width: 0; margin: 1rem 0 }
      figure img { display: block; width: auto; height: auto; max-width: min(100%, 750px);
        max-height: 75vh; margin-inline: auto }
      figcaption { color: #455b63; font-size: .91rem } summary { cursor: pointer; font-weight: 650 }
      details { border-top: 1px solid #d7e1de; padding: .7rem 0 } a { color: #226058 }
      code { overflow-wrap: anywhere }"))),
  tags$body(tags$main(tags$h1("Bulk RNA-seq 分析报告"),
    tags$p(sprintf("%d 个样本 · %d 项比较", nrow(samples), length(comparisons))),
    tags$section(tags$h2("比较概览"), table_html(do.call(rbind, summary_rows)),
      tags$p("DESeq2 设计式 ", tags$code(p$design), "。差异基因要求 Wald 检验的 BH 校正 padj < ",
             p$alpha, "，且 apeglm 收缩后的 |log2FC| ≥ log2(", p$min_fold_change, ")。"),
      tags$p("ORA 按上下调分别检验，候选背景为可判定差异的基因，并与各数据库注释取交集。",
             "GSEA 按 Wald 统计量降序排列。通路显著阈值为 BH 校正后 P < ", p$pathway_fdr, "。"),
      source_note), qc, sections)))
# 完整文档需保留 head 中的视口声明和样式。
writeLines(enc2utf8(c("<!doctype html>", htmltools::doRenderTags(page))), output, useBytes = TRUE)
