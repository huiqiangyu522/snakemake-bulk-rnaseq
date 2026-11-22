library(ggplot2)

res <- read.delim(snakemake@input[["res"]], colClasses = c(gene_id = "character"),
                  check.names = FALSE)
deg <- read.delim(snakemake@input[["deg"]], colClasses = c(gene_id = "character"),
                  check.names = FALSE)
cutoff <- snakemake@params[["pathway_fdr"]]
term_file <- snakemake@input[["term2gene"]]
custom <- length(term_file) > 0L
ids <- if (custom) res$gene_id else sub("\\.[0-9]+$", "", res$gene_id)
if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
  stop("富集输入的基因 ID 含空值或重复值，请检查版本号去除后的 ID。")
}

eligible <- is.finite(res$padj) & is.finite(res$log2FC_shrunk)
universe <- ids[eligible]
deg_ids <- ids[match(deg$gene_id, res$gene_id)]
foreground <- list(up = deg_ids[deg$direction == "Up"], down = deg_ids[deg$direction == "Down"])

empty_plot <- function(title, label) {
  ggplot() + annotate("text", x = 0, y = 0, label = label) + labs(title = title) + theme_void()
}

# 输出所有返回的条目；阈值只作用于 significant 标记和绘图。
result_table <- function(result, empty, key) {
  tbl <- as.data.frame(result)
  if (!nrow(tbl)) tbl <- empty
  tbl$significant <- is.finite(tbl$p.adjust) & tbl$p.adjust < cutoff
  tbl <- tbl[order(tbl$p.adjust, tbl$ID), , drop = FALSE]
  readr::write_tsv(tbl, snakemake@output[[paste0(key, "_table")]])
  tbl
}

save_ora <- function(key, title, genes, background, enrich, ...) {
  result <- if (length(genes)) enrich(gene = genes, universe = background, ...,
                                      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  empty <- data.frame(ID = character(), Description = character(), GeneRatio = character(),
                      BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
                      qvalue = numeric(), geneID = character(), Count = integer())
  tbl <- result_table(result, empty, key)
  top <- head(tbl[tbl$significant, , drop = FALSE], 20)
  if (nrow(top)) {
    ratio <- strsplit(as.character(top$GeneRatio), "/", fixed = TRUE)
    top$ratio <- vapply(ratio, function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
    top$label <- factor(top$ID, levels = rev(top$ID))
    labels <- setNames(stringr::str_wrap(top$Description, width = 42), top$ID)
    p <- ggplot(top, aes(ratio, label, size = Count, colour = p.adjust)) +
      geom_point() + scale_y_discrete(labels = labels) +
      scale_colour_gradient(low = "#B9413D", high = "#386CA0") +
      labs(title = title, x = "Gene ratio", y = NULL, colour = "Adjusted P") + theme_bw()
  } else {
    note <- if (!length(genes)) "No foreground genes" else if (!nrow(tbl)) "No testable terms returned" else
      paste("No pathways with adjusted P <", cutoff)
    p <- empty_plot(title, note)
  }
  ggsave(snakemake@output[[paste0(key, "_plot")]], p, width = 9, height = 8, dpi = 300)
}

if (custom) {
  t2g <- read.delim(term_file, colClasses = "character", check.names = FALSE)
  if (!all(c("term", "gene_id") %in% names(t2g))) {
    stop("本地基因集表须包含 term 和 gene_id 两列。")
  }
  t2g <- unique(t2g[, c("term", "gene_id")])
  if (anyNA(t2g) || any(!nzchar(t2g$term)) || any(!nzchar(t2g$gene_id)) ||
      !any(ids %in% t2g$gene_id)) {
    stop("本地基因集存在空值或与差异表没有共同的 gene_id。请检查 enrichment.term2gene。")
  }
  for (direction in names(foreground)) {
    save_ora(paste0("custom_", direction), paste("Custom ORA |", direction),
             foreground[[direction]], universe, clusterProfiler::enricher, TERM2GENE = t2g)
  }
  gsea_title <- "Custom GSEA"
} else {
  if (all(grepl("^ENSG[0-9]{11}$", ids))) {
    org_package <- "org.Hs.eg.db"
    organism <- "hsa"
    db_species <- "HS"
    collection <- "H"
    species <- "Homo sapiens"
  } else if (all(grepl("^ENSMUSG[0-9]{11}$", ids))) {
    org_package <- "org.Mm.eg.db"
    organism <- "mmu"
    db_species <- "MM"
    collection <- "MH"
    species <- "Mus musculus"
  } else {
    stop("公共数据库模式需要统一的人或小鼠 Ensembl 基因 ID。模拟基因请配置 enrichment.term2gene。")
  }
  org_db <- getExportedValue(org_package, org_package)
  pairs <- AnnotationDbi::select(org_db, keys = ids, keytype = "ENSEMBL", columns = "ENTREZID")
  pairs <- unique(pairs[!is.na(pairs$ENTREZID), c("ENSEMBL", "ENTREZID")])
  ambiguous <- function(x) duplicated(x) | duplicated(x, fromLast = TRUE)
  mapping <- pairs[!ambiguous(pairs$ENSEMBL) & !ambiguous(pairs$ENTREZID), ]
  entrez <- setNames(as.character(mapping$ENTREZID), mapping$ENSEMBL)
  convert <- function(x) unname(entrez[intersect(x, names(entrez))])
  kegg_background <- convert(universe)
  message("物种 ", species, "；ORA 候选背景 ", length(universe),
          "；可一对一映射到 Entrez 的背景基因 ", length(kegg_background))

  kegg_readable <- function(...) {
    result <- clusterProfiler::enrichKEGG(..., organism = organism, keyType = "ncbi-geneid")
    if (!is.null(result) && nrow(as.data.frame(result))) {
      result <- clusterProfiler::setReadable(result, OrgDb = org_db, keyType = "ENTREZID")
    }
    result
  }
  for (direction in names(foreground)) {
    save_ora(paste0("go_bp_", direction), paste("GO BP |", direction),
             foreground[[direction]], universe, clusterProfiler::enrichGO,
             OrgDb = org_db, keyType = "ENSEMBL", ont = "BP", readable = TRUE)
    save_ora(paste0("kegg_", direction), paste("KEGG |", direction),
             convert(foreground[[direction]]), kegg_background, kegg_readable)
  }
  hallmark <- msigdbr::msigdbr(db_species = db_species, species = species, collection = collection)
  t2g <- unique(data.frame(term = hallmark$gs_name, gene_id = hallmark$ensembl_gene))
  t2g <- t2g[!is.na(t2g$gene_id) & nzchar(t2g$gene_id), ]
  gsea_title <- "Hallmark GSEA"
}

# GSEA 使用所有具有有效 Wald 统计量和原始 p 值的基因。
keep <- is.finite(res$stat) & is.finite(res$pvalue)
gene_list <- sort(setNames(res$stat[keep], ids[keep]), decreasing = TRUE)
sizes <- table(t2g$term[t2g$gene_id %in% names(gene_list)])
gsea <- NULL
max_size <- min(500L, length(gene_list) - 1L)
if (any(sizes >= 15 & sizes <= max_size) && length(unique(gene_list)) > 1L) {
  set.seed(123)
  gsea <- clusterProfiler::GSEA(geneList = gene_list, TERM2GENE = t2g,
                                pAdjustMethod = "BH", pvalueCutoff = 1,
                                minGSSize = 15, maxGSSize = max_size, verbose = FALSE)
}
empty <- data.frame(ID = character(), Description = character(), setSize = integer(),
                    enrichmentScore = numeric(), NES = numeric(), pvalue = numeric(),
                    p.adjust = numeric(), core_enrichment = character())
tbl <- result_table(gsea, empty, "gsea")
sig <- tbl[tbl$significant & is.finite(tbl$NES), , drop = FALSE]
if (nrow(sig)) {
  top <- rbind(head(sig[sig$NES > 0, ], 10), head(sig[sig$NES < 0, ], 10))
  top$label <- gsub("_", " ", sub("^HALLMARK_", "", top$Description))
  p <- ggplot(top, aes(NES, reorder(label, NES), fill = NES > 0)) +
    geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey50") +
    scale_fill_manual(values = c("TRUE" = "#B9413D", "FALSE" = "#386CA0"), guide = "none") +
    labs(title = gsea_title, x = "Normalized enrichment score (NES)", y = NULL) + theme_bw()
} else {
  note <- if (nrow(tbl)) paste("No pathways with adjusted P <", cutoff) else "No testable gene sets returned"
  p <- empty_plot(gsea_title, note)
}
ggsave(snakemake@output[["gsea_plot"]], p, width = 10, height = 8, dpi = 300)
