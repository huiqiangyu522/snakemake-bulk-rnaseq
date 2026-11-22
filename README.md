# Bulk RNA-seq 差异表达与富集分析工作流

从每个样本的**基因水平原始计数**出发，用 Snakemake 串联输入校验、样本质控、DESeq2 差异表达、功能富集和 HTML 报告。一次运行可以处理多项组间比较。已完成比对和定量的 bulk RNA-seq 数据可以直接从这一步接入。

[查看示例报告](results/report.html) · [查看默认配置](config.yaml) · [查看示例数据说明](c0_fixture/README.md)

## 快速运行

在包含 `Snakemake` 和 `config.yaml` 的项目根目录执行。首次使用请先按下一节安装依赖。入口文件名为 `Snakemake`，用 `-s Snakemake` 指定。

```bash
conda activate snakemake
snakemake -s Snakemake --cores 2 -n
snakemake -s Snakemake --cores 2
```

第一条 Snakemake 命令预览运行计划，第二条执行分析。完成后用浏览器打开 `results/report.html`。图片与可下载的 TSV 均嵌在 HTML 中，可以单独发送这份报告。报告中的富集面板可点击展开。

**示例复核。** 默认数据包含 3 个配对供体、6 个样本和 120 个模拟基因。使用随附配置和下述测试环境运行 `LPS_vs_Control`，得到 18 个上调基因和 18 个下调基因。两张 ORA 气泡图各有一个点，分别对应 `SIM_LPS_UP` 和 `SIM_LPS_DOWN`。这些模拟条目用于检查比较方向与文件衔接。

## 安装与环境

已运行示例的环境为 macOS arm64、Snakemake 9.27.1、R 4.6.1 和 Bioconductor 3.23，主要分析包为 DESeq2 1.52.0、apeglm 1.34.0 和 clusterProfiler 4.20.3。下面的命令适用于用 Conda 管理 Snakemake、用本机 R 管理分析包的安装方式。

先安装 R，并确认终端可以找到 `Rscript`。首次创建 Snakemake 环境时执行：

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake=9.27.1
conda activate snakemake
snakemake --version
Rscript --version
```

在该终端启动 `R`，将以下包装入工作流运行时使用的 R 库。已有环境可以直接使用已安装的包。

```r
install.packages(c(
  "BiocManager", "readr", "tibble", "ggplot2", "ggrepel", "pheatmap",
  "matrixStats", "stringr", "htmltools", "knitr", "base64enc"
))
BiocManager::install(c("DESeq2", "apeglm", "clusterProfiler"))
```

使用公共数据库做富集时，额外安装 `msigdbr`、`AnnotationDbi` 和对应物种的注释包。人类数据示例如下，小鼠将最后一个包名换为 `org.Mm.eg.db`。

```r
install.packages("msigdbr")
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))
```

## 换成自己的数据

### 1. 准备样本表和计数文件

样本表为带表头的 TSV，每行一个样本。`sample_id` 必须唯一，`count_file` 指向该样本独立的文件，设计式引用的每一列都要出现在表中。相对路径从**项目根目录**解析。下面是默认样本表的前两行。

```tsv
sample_id	sampletype	donor	count_file
Control_1	Control	D1	c0_fixture/raw/normal/Control_1.counts.tsv.gz
LPS_1	LPS	D1	c0_fixture/raw/normal/LPS_1.counts.tsv.gz
```

每个计数文件为 `.tsv` 或 `.tsv.gz`，含 `Geneid`、`gene_name` 和一列非负整数原始计数。下例为 `Control_1` 的文件节选。

```tsv
Geneid	gene_name	Control_1
SIM0002	SYN_GENE_0002	78
SIM0084	SYN_GENE_0084	259
```

所有样本应包含相同的 `Geneid` 集合，行顺序可以不同。流程按 `Geneid` 对齐，允许不同 ID 共用同一个 `gene_name`；同一样本内的重复 `Geneid` 会触发报错。若已有多样本合并矩阵，需先整理成上述每样本一份的格式。若手头只有 FASTQ，应先完成上游定量。

### 2. 配置设计式与比较方向

编辑 [config.yaml](config.yaml)。以下为换用真实数据的配对设计示例，需将路径和组名换成自己的值。

```yaml
samples: "metadata/samples.tsv"

deseq2:
  design: "~ donor + sampletype"
  group_column: "sampletype"
  alpha: 0.05
  min_fold_change: 1.5
  contrasts:
    LPS_vs_Control:
      numerator: "LPS"
      denominator: "Control"

enrichment:
  pathway_fdr: 0.05
  term2gene: null
```

- **比较方向**由 `numerator` 和 `denominator` 决定。上例中，正 log2FC 表示校正供体差异后，LPS 组表达较高。`LPS_vs_Control` 是用于命名结果文件的比较名。
- **独立样本**可使用 `design: "~ group"` 和 `group_column: "group"`，并在样本表中提供 `group` 列。当前工作流使用含截距的模型，分组变量以主效应进入设计式。
- **多项比较**可在 `contrasts` 下并列增加条目。所有比较共享基础 DDS 和一次样本质控；差异分析按各自的分母组设置参考水平并拟合。切换数据集时，同名比较会使用相同输出路径；需要保留两套结果时，可在不同项目目录中运行。

### 3. 选择富集基因集

| `enrichment.term2gene` | 适用输入 | 运行内容 |
|---|---|---|
| 本地 TSV 路径，示例为 `c0_fixture/term2gene.tsv` | 基因集表含 `term` 和 `gene_id` 两列，ID 与计数文件的 `Geneid` 一致 | 以该基因集分别做上调 ORA、下调 ORA 和 GSEA |
| `null` | 人或小鼠的 Ensembl 基因 ID，同一次分析使用同一物种 | GO Biological Process、KEGG 分方向做 ORA，Hallmark 做 GSEA |

公共数据库模式会去掉 Ensembl ID 的版本后缀。GO 使用 Ensembl ID，KEGG 使用一对一映射的 Entrez ID。KEGG 需要访问 KEGG REST，首次获取未缓存的 Hallmark 基因集也需要网络。用真实数据切换模式时，应同时更换样本表和 `term2gene`；示例的 `SIM...` ID 只和本地模拟基因集对应。

## 分析方法与主要输出

`validate_inputs` 先按基因 ID 合并计数，`prepare_dds` 创建共享 DDS 并估计归一化因子。随后 `sample_qc` 对全体样本做质控，`deseq2`、`de_summary` 和 `enrichment` 按比较生成结果，最后由 `report` 汇总。各规则对应的 R 脚本位于 [script/](script/)。

| 文件或目录 | 如何使用 |
|---|---|
| `work/validated/counts.tsv`、`work/validated/samples.tsv` | 核对合并后的计数矩阵与样本表。共享基础对象保存在 `work/deseq2/dds_base.rds`。 |
| `results/figures/sample_qc/` | 查看全体样本经盲法 rlog 转换后的 PCA、Pearson 相关性与层次聚类。 |
| `results/tables/deseq2/{comparison}.tsv` | 全量差异表。`stat`、`pvalue`、`padj` 来自 Wald 检验，同时保留原始和 apeglm 收缩后的 log2FC。 |
| `results/tables/de_summary/{comparison}/deg.tsv` | 在默认阈值下，选取 `padj < 0.05` 且 `abs(log2FC_shrunk) >= log2(1.5)` 的基因，`direction` 列标记 Up 或 Down。 |
| `results/figures/de_summary/{comparison}/` | 火山图和差异基因热图。热图展示当前比较的两组样本，使用归一化计数经 `log2(x + 1)` 转换后按基因标准化。 |
| `results/tables/enrichment/{comparison}/`、`results/figures/enrichment/{comparison}/` | 每项富集分析的结果表与配图。表中保留返回条目，以 `significant` 标记达到通路 FDR 阈值的条目。 |
| `results/report.html` | 汇总所有比较的样本信息、质控图、差异基因、富集结果和 TSV 下载。 |

ORA 分别使用上调和下调差异基因，候选背景取自具有有效 `padj` 和收缩 log2FC 的基因，再与对应注释取交集。GSEA 使用所有具有有效 Wald 统计量和原始 P 值的基因，按 `stat` 降序排列。没有差异基因或没有显著条目时，流程仍会生成带相应说明的报告。

## 验证

只检查输入时可运行 `snakemake -s Snakemake --cores 1 validate_inputs`。示例的合并矩阵按 `gene_id` 排序后，应与 `c0_fixture/counts.tsv` 逐值相同。另有重复基因符号和重复基因 ID 两套校验用例，运行命令及预期行为见[测试数据说明](c0_fixture/README.md)。

当前模拟示例已从计数文件运行至 HTML 报告。人和小鼠的 GO 与 Hallmark 分支已用对应的 Ensembl ID 检查；在线 KEGG 查询仍待在可访问 KEGG REST 的环境中完成验证。
