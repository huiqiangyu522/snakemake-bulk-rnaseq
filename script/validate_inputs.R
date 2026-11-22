# 读取
counts_paths <- as.character(snakemake@input[["counts"]])
samples <- readr::read_tsv(
  snakemake@input[["samples"]],
  col_types = readr::cols(.default = readr::col_character()),
  name_repair = "check_unique"
)
# 错误时终止
readr::stop_for_problems(samples)
# 检查samples表是否有缺失信息
design_cols <- all.vars(as.formula(snakemake@params[["design"]]))
required_cols <- unique(c("sample_id", "count_file", design_cols))
missing_cols <- setdiff(required_cols, names(samples))

if (length(missing_cols) > 0L) {
  stop("样本表缺少列：", paste(missing_cols, collapse = ", "))
}

if (nrow(samples) == 0L || anyNA(samples[required_cols])) {
  stop("样本表为空，或必填列中存在缺失值。")
}
if (anyDuplicated(samples$sample_id) > 0L || "gene_id" %in% samples$sample_id) {
  stop("sample_id 必须唯一，且不能使用保留列名 gene_id。")
}

# 确认传入的计数文件与样本表逐行对应
resolved_paths <- normalizePath(counts_paths, mustWork = TRUE)
if (anyDuplicated(resolved_paths)) {
  stop("多个样本引用了同一个计数文件，请检查 count_file。")
}
if (!identical(
  resolved_paths,
  normalizePath(samples$count_file, mustWork = TRUE)
)) {
  stop("Snakemake 传入的计数文件与样本表的 count_file 列不一致。")
}

# 定义函数
read_sample_counts <- function(path, sample_id) {
  tbl <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "check_unique"
  )
  readr::stop_for_problems(tbl)
  
  id_cols <- c("Geneid", "gene_name")
  count_col <- setdiff(names(tbl), id_cols)
  if (!all(id_cols %in% names(tbl)) || length(count_col) != 1L) {
    stop(sample_id, " 的计数文件应包含 Geneid、gene_name 和一个计数列。")
  }
  
  gene_ids <- tbl$Geneid
  if (length(gene_ids) == 0L || anyNA(gene_ids)) {
    stop(sample_id, " 的计数文件为空，或 Geneid 存在缺失值。")
  }
  
  duplicate_ids <- unique(gene_ids[duplicated(gene_ids)])
  if (length(duplicate_ids) > 0L) {
    stop(sample_id, " 的 Geneid 重复：",
         paste(head(duplicate_ids, 5), collapse = ", "))
  }
  
  counts <- suppressWarnings(as.numeric(tbl[[count_col]]))
  if (any(!is.finite(counts) | counts < 0 | counts > .Machine$integer.max |
          counts != floor(counts))) {
    stop(sample_id, " 含有无效计数，计数应为非负整数。")
  }
  
  list(gene_id = gene_ids, counts = counts)
}
# 读取计数表
counts_list <- Map(read_sample_counts, counts_paths, samples$sample_id)
# 检查计数表的基因id是否一致并合并
reference_ids <- counts_list[[1]]$gene_id
merged_counts <- tibble::tibble(gene_id = reference_ids)

for (i in seq_along(counts_list)) {
  sample_counts <- counts_list[[i]]
  sample_id <- samples$sample_id[[i]]
  
  if (!setequal(reference_ids, sample_counts$gene_id)) {
    stop(sample_id, " 的 Geneid 集合与 ", samples$sample_id[[1]], " 不一致。")
  }
  
  row_index <- match(reference_ids, sample_counts$gene_id)
  merged_counts[[sample_id]] <- sample_counts$counts[row_index]
}

readr::write_tsv(merged_counts, snakemake@output[["counts"]])
readr::write_tsv(samples, snakemake@output[["samples"]])

message("校验通过：", nrow(merged_counts), " 个基因，", nrow(samples), " 个样本。")
