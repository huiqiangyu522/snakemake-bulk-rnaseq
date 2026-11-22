# C0 RNA-seq 测试数据

这套合成数据包含 6 个样本、3 个配对供体和 120 个模拟基因。每个供体各有一个 Control 样本和一个 LPS 样本，适合测试设计式 `~ donor + sampletype` 和比较方向 `LPS vs Control`。原始输入是 gzip 压缩的制表符文本。每个文件包含 `Geneid`、`gene_name` 和一个计数列，计数列以该样本的 `sample_id` 命名。不同文件的基因行顺序有意打乱，便于检查合并时是否按 `Geneid` 对齐。

## 使用方法

在包含 `Snakemake` 的项目根目录解压，使该目录下出现 `c0_fixture/`。下表的相对路径均以项目根目录为起点。工作流通过 `config["samples"]` 取得样本表路径。

| 样本表 | 预期校验行为 |
|---|---|
| `c0_fixture/samples.tsv` | 成功。合并结果应有 120 个基因和 6 个样本列，可与包中的 `counts.tsv` 按 `gene_id` 排序后逐值比较。|
| `c0_fixture/samples_duplicate_symbol.tsv` | 成功。`SIM0001` 和 `SIM0002` 共享 `gene_name`（`SYN_GENE_0001`），仍应保留两行。|
| `c0_fixture/samples_duplicate_geneid.tsv` | 失败。`Control_1` 的计数文件多了一行 `SIM0001`，报错应指出样本和重复 ID。|

```bash
snakemake -s Snakemake --cores 1 validate_inputs --config samples=c0_fixture/samples.tsv
snakemake -s Snakemake --cores 1 validate_inputs --config samples=c0_fixture/samples_duplicate_symbol.tsv
snakemake -s Snakemake --cores 1 validate_inputs --config samples=c0_fixture/samples_duplicate_geneid.tsv
```

若你的 C0 已以单个 counts 矩阵为入口，可直接使用 `c0_fixture/counts.tsv` 和正常样本表。`annotation.tsv` 提供正常用例的 ID 与基因名映射。

## 下游验证

- 对于 `LPS vs Control`，`SIM0001`–`SIM0018` 注入了上调信号，`SIM0019`–`SIM0036` 注入了下调信号，其余基因未注入分组效应。先核对 log2FC 方向，显著基因数量会随模型和阈值变化。
- `SIM...` 和 `SYN_GENE_...` 是合成标识。真实 GO、KEGG 或物种注释库无法映射它们。若要测试富集表和报告的文件衔接，可将正常差异结果中的 `gene_id` 与 `term2gene.tsv` 配套使用，例如传给 `clusterProfiler::enricher()` 的 `TERM2GENE`，并用参与检验的模拟基因作背景。
- `term2gene.tsv` 中 `SIM_LPS_UP` 和 `SIM_LPS_DOWN` 分别对应注入的上调和下调集合。所有数值和通路名称仅用于检查工作流。
