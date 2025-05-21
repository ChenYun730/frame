**整理数据**
将共享文档中整理的条目按照如下方式进行分析整理
```
data/
├── SARS-CoV-2/
│   ├── A549/
│   │   ├── 24h_MOI1.0_GSE123456/
│   │   │   ├── expression_matrix.csv
|   |   |   |—— deseq2_results.csv
│   │   │   ├── volcano.png
│   │   │   ├── heatmap.png
│   │   │   └── meta.json(暂时不做）
│   │   ├── 48h_MOI1.0_GSE123789/
│   │   │   └── ...
│   └── THP-1/
│       └── ...
└── IAV/
    └── ...
```
1、建立文件夹，获得count的csv文件
```
dir.create("sampledata")
setwd("sampledata")   #进行后续的子文件夹创建时要加入绝对路径以免找不到文件
dirs_to_create <- c("sars_cov_2", "sars_cov_1", "rsv", "hpiv3", "iav")
sapply(dirs_to_create, dir.create)
setwd("/mnt/alamo01/users/chenyun730/sampledata/sars_cov_1")
dir1<- c("Calu-3_2B4")
sapply(dir1, function(x) dir.create(x, recursive = TRUE))

setwd("/mnt/alamo01/users/chenyun730/sampledata/sars_cov_2")
dir2<- c("Calu-3_2B4","NHBE","A549","A549-ACE2","Calu-3")
sapply(dir2, function(x) dir.create(x, recursive = TRUE))

setwd("/mnt/alamo01/users/chenyun730/sampledata/rsv")
dir3<- c("A549")
sapply(dir3, function(x) dir.create(x, recursive = TRUE))

setwd("/mnt/alamo01/users/chenyun730/sampledata/hpiv3")
dir4<- c("A549")
sapply(dir4, function(x) dir.create(x, recursive = TRUE))

setwd("/mnt/alamo01/users/chenyun730/sampledata/iav")
dir5<- c("A549","NHBE")
sapply(dir5, function(x) dir.create(x, recursive = TRUE))

# 三个压缩包及命名
raw_data_human <- read.table("/mnt/alamo01/users/chenyun730/data/GSE147507_RawReadCounts_Human.tsv.gz",
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE)
raw_data_ferret <- read.table("/mnt/alamo01/users/chenyun730/data/GSE147507_RawReadCounts_Ferret.tsv.gz",
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE)
raw_data_all <- read.table("/mnt/alamo01/users/chenyun730/data/GSE255647_all.counts.tsv.gz",
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE)

 head(colnames(raw_data_all), 100)
 [1] "MOCK_12hpi_rep1"   "MOCK_12hpi_rep2"   "MOCK_12hpi_rep3"
 [4] "MOCK_24hpi_rep1"   "MOCK_24hpi_rep2"   "MOCK_24hpi_rep3"
 [7] "MOCK_48hpi_rep1"   "MOCK_48hpi_rep2"   "MOCK_48hpi_rep3"
[10] "SARS-1_12hpi_rep1" "SARS-1_12hpi_rep2" "SARS-1_12hpi_rep3"
[13] "SARS-1_24hpi_rep1" "SARS-1_24hpi_rep2" "SARS-1_24hpi_rep3"
[16] "SARS-1_48hpi_rep1" "SARS-1_48hpi_rep2" "SARS-1_48hpi_rep3"
[19] "SARS-2_12hpi_rep1" "SARS-2_12hpi_rep2" "SARS-2_12hpi_rep3"
[22] "SARS-2_24hpi_rep1" "SARS-2_24hpi_rep2" "SARS-2_24hpi_rep3"
[25] "SARS-2_48hpi_rep1" "SARS-2_48hpi_rep2" "SARS-2_48hpi_rep3"

 head(colnames(raw_data_human), 110)
 [1] "Series1_NHBE_Mock_1"                 "Series1_NHBE_Mock_2"
 [3] "Series1_NHBE_Mock_3"                 "Series1_NHBE_SARS-CoV-2_1"
 [5] "Series1_NHBE_SARS-CoV-2_2"           "Series1_NHBE_SARS-CoV-2_3"
 [7] "Series2_A549_Mock_1"                 "Series2_A549_Mock_2"
 [9] "Series2_A549_Mock_3"                 "Series2_A549_SARS-CoV-2_1"
[11] "Series2_A549_SARS-CoV-2_2"           "Series2_A549_SARS-CoV-2_3"
[13] "Series3_A549_Mock_1"                 "Series3_A549_Mock_2"
[15] "Series3_A549_RSV_1"                  "Series3_A549_RSV_2"
[17] "Series4_A549_Mock_1"                 "Series4_A549_Mock_2"
[19] "Series4_A549_IAV_1"                  "Series4_A549_IAV_2"
[21] "Series5_A549_Mock_1"                 "Series5_A549_Mock_2"
[23] "Series5_A549_Mock_3"                 "Series5_A549_SARS-CoV-2_1"
[25] "Series5_A549_SARS-CoV-2_2"           "Series5_A549_SARS-CoV-2_3"
[27] "Series6_A549-ACE2_Mock_1"            "Series6_A549-ACE2_Mock_2"
[29] "Series6_A549-ACE2_Mock_3"            "Series6_A549-ACE2_SARS-CoV-2_1"
[31] "Series6_A549-ACE2_SARS-CoV-2_2"      "Series6_A549-ACE2_SARS-CoV-2_3"
[33] "Series7_Calu3_Mock_1"                "Series7_Calu3_Mock_2"
[35] "Series7_Calu3_Mock_3"                "Series7_Calu3_SARS-CoV-2_1"
[37] "Series7_Calu3_SARS-CoV-2_2"          "Series7_Calu3_SARS-CoV-2_3"
[39] "Series8_A549_Mock_1"                 "Series8_A549_Mock_2"
[41] "Series8_A549_Mock_3"                 "Series8_A549_RSV_1"
[43] "Series8_A549_RSV_2"                  "Series8_A549_RSV_3"
[45] "Series8_A549_HPIV3_3"                "Series8_A549_HPIV3_2"
[47] "Series8_A549_HPIV3_1"                "Series9_NHBE_Mock_1"
[49] "Series9_NHBE_Mock_2"                 "Series9_NHBE_Mock_3"
[51] "Series9_NHBE_Mock_4"                 "Series9_NHBE_IAV_1"
[53] "Series9_NHBE_IAV_2"                  "Series9_NHBE_IAV_3"
[55] "Series9_NHBE_IAV_4"                  "Series9_NHBE_IAVdNS1_1"
[57] "Series9_NHBE_IAVdNS1_2"              "Series9_NHBE_IAVdNS1_3"
[59] "Series9_NHBE_IAVdNS1_4"              "Series9_NHBE_IFNB_4h_1"
[61] "Series9_NHBE_IFNB_4h_2"              "Series9_NHBE_IFNB_6h_1"
[63] "Series9_NHBE_IFNB_6h_2"              "Series9_NHBE_IFNB_12h_1"
[65] "Series9_NHBE_IFNB_12h_2"             "Series15_HealthyLungBiopsy_2"
[67] "Series15_HealthyLungBiopsy_1"        "Series15_COVID19Lung_2"
[69] "Series15_COVID19Lung_1"              "Series16_A549-ACE2_Mock_1"
[71] "Series16_A549-ACE2_Mock_2"           "Series16_A549-ACE2_Mock_3"
[73] "Series16_A549-ACE2_SARS-CoV-2_1"     "Series16_A549-ACE2_SARS-CoV-2_2"
[75] "Series16_A549-ACE2_SARS-CoV-2_3"     "Series16_A549-ACE2_SARS-CoV-2_Rux_1"
[77] "Series16_A549-ACE2_SARS-CoV-2_Rux_2" "Series16_A549-ACE2_SARS-CoV-2_Rux_3"

以SARS-CoV-1感染Calu-3/2B4 1 MOI 12h为例：
mock_samples <- grep("^MOCK_12hpi_rep[123]$", colnames(raw_data_all), value=TRUE)
virus_samples <- grep("^SARS[-.]1_12hpi_rep[123]$", colnames(raw_data_all), value=TRUE)
count_matrix <- data.frame(
  Gene = rownames(raw_data_all), 
  raw_data_all[, c(mock_samples, virus_samples)],
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
col_data <- data.frame(
  sample = c(mock_samples, virus_samples),
  condition = rep(c("Mock", "Virus"), each=3),
  timepoint = "12hpi",
  virus_type = c(rep(NA, 3), rep("SARS-1", 3)),
  stringsAsFactors = FALSE
)
setwd("/mnt/alamo01/users/chenyun730/sampledata/sars_cov_1/Calu-3_2B4")
write.csv(count_matrix,"12H_MOI1_GSE255647.csv",row.names=FALSE)


 mock_samples <- c("MOCK_12hpi_rep1","MOCK_12hpi_rep2","MOCK_12hpi_rep3")
virus_samples <- c("SARS-1_12hpi_rep1","SARS-1_12hpi_rep2","SARS-1_12hpi_rep3")
count_matrix <- data.frame(
  Gene = rownames(raw_data_human),
  raw_data_human[, c(mock_samples, virus_samples)],
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
col_data <- data.frame(
  sample = c(mock_samples, virus_samples),
  condition = rep(c("Mock", "Virus"), each=3),
  timepoint = "24hpi",
  virus_type = c(rep(NA, 3), rep(SARS-1", 3)),
  stringsAsFactors = FALSE
)
setwd("/mnt/alamo01/users/chenyun730/sampledata/hpiv3/A549/sars_cov_2/Calu-3_2B4/12H_MOI1_GSE255647")
write.csv(count_matrix,"expression_matrix.csv",row.names=FALSE)


```

2. DEseq2得到deseq2_results.csv
```
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(EnhancedVolcano)
library(RColorBrewer)

setwd("/mnt/alamo01/users/chenyun730/sampledata/sars_cov_1/Calu-3_2B4/12H_MOI1_GSE255647")
count_data <- read.csv("expression_matrix.csv", row.names = 1, check.names = FALSE)
col_data <- data.frame(
  condition = factor(rep(c("Mock", "Virus"), each = 3)),
  row.names = colnames(count_data)
)
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = col_data,
  design = ~ condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]  # 保留总counts≥10的基因
dds <- DESeq(dds)
results <- results(dds, contrast = c("condition", "Virus", "Mock"))
write.csv(as.data.frame(results), "/mnt/alamo01/users/chenyun730/sampledata/sars_cov_1/Calu-3_2B4/12H_MOI1_GSE255647/deseq2_results.csv")

```
3. 绘制火山图（使用EnhancedVolcano美化）

```
rld <- rlog(dds, blind = FALSE)  #数据转换（rlog用于下游可视化）

volcano_plot <- EnhancedVolcano(
  results,
  lab = rownames(results),
  x = 'log2FoldChange',
  y = 'pvalue',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 4.0,
  colAlpha = 0.7,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey50',
  title = 'Sars-Cov-1 vs Mock (12H MOI1)',
  caption = 'FC cutoff: 1; p-value cutoff: 0.05',
  gridlines.major = FALSE,
  gridlines.minor = FALSE
) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("grey30", "forestgreen", "royalblue", "red2"))
ggsave("volcano.png", plot = volcano_plot, width = 10, height = 8, dpi = 300)
```
4. 绘制热图（选取差异最显著的30个基因）
```
norm_counts <- assay(rld)
top_genes <- rownames(results)[order(results$padj)][1:30]
annotation_col <- data.frame(
  Condition = col_data$condition,
  row.names = colnames(norm_counts)
)
heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
heatmap_plot <- pheatmap(
  norm_counts[top_genes, ],
  scale = "row",
  color = heatmap_colors,  # 使用预定义的颜色梯度
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 10,
  cellwidth = 40,
  cellheight = 10,
  angle_col = 45,
  main = "Top 30 DEGs in Calu-3/2B4\n(12hpi, SARS-CoV-1 vs Mock)",
  gaps_col = 3,
  silent = FALSE
)
png(
  "heatmap.png",
  width = 8,  
  height = 8,
  units = "in",
  res = 300 
)
print(heatmap_plot)
dev.off()
```
5、meta.json
```
library(jsonlite)
meta_info <- list(
  experiment = list(
    virus = "SARS-CoV-1",
    cell_line = "Calu-3/2B4",  # 根据您的数据调整
    infection_time = "12h",     # 根据您的MOI命名调整
    MOI = "1.0",
    replicates = 3,
    control_condition = "Mock"
  ),
  data_source = list(
    accession = "GSE147507",    # 根据您的数据调整
    comparison = "SARS-CoV-1 infected vs Mock",
    organism = "Homo sapiens",
    reference = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507"
  ),
  analysis = list(
    tool = "DESeq2",
    version = paste0("R-", R.version$major, ".", R.version$minor),
    date = format(Sys.Date(), "%Y-%m-%d"),
    parameters = list(
      pvalue_cutoff = 0.05,
      log2fc_cutoff = 1,
      normalization = "rlog")))
output_file <- "/mnt/alamo01/users/chenyun730/sampledata/sars_cov_1/Calu-3_2B4/12H_MOI1_GSE255647/meta_analysis_info.json"
write_json(
  meta_info,
  path = output_file,
  pretty = TRUE,
  auto_unbox = TRUE  # 自动解包单元素向量
)
```
批量处理的模板：
```

setwd("/mnt/alamo01/users/chenyun730/sampledata/sars_cov_2/Calu-3_2B4/12H_MOI1_GSE255647/")
count_data <- read.csv("expression_matrix.csv", row.names = 1, check.names = FALSE)
head(count_data)

col_data <- data.frame(
  condition = factor(rep(c("Mock", "Virus"), each = 3)),
  row.names = colnames(count_data)
)
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = col_data,
  design = ~ condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
results <- results(dds, contrast = c("condition", "Virus", "Mock"))
write.csv(as.data.frame(results), "deseq2_results.csv")

rld <- rlog(dds, blind = FALSE)
volcano_plot <- EnhancedVolcano(
  results,
  lab = rownames(results),
  x = 'log2FoldChange',
  y = 'pvalue',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 4.0,
  colAlpha = 0.7,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey50',
  title = 'SARS-CoV-2 vs Mock (12H MOI1)',
  caption = 'FC cutoff: 1; p-value cutoff: 0.05',
  gridlines.major = FALSE,
  gridlines.minor = FALSE
) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("grey30", "forestgreen", "royalblue", "red2"))
ggsave("volcano.png", plot = volcano_plot, width = 10, height = 8, dpi = 300)

norm_counts <- assay(rld)
top_genes <- rownames(results)[order(results$padj)][1:30]
annotation_col <- data.frame(
  Condition = col_data$condition,
  row.names = colnames(norm_counts)
)

ann_colors <- list(Condition = c(Mock = "#66C2A5", Virus = "#FC8D62"))
heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
heatmap_plot <- pheatmap(
  norm_counts[top_genes, ],
  scale = "row",
  color = heatmap_colors,
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 10,
  cellwidth = 40,
  cellheight = 10,
  angle_col = 45,
  main = "Top 30 DEGs in Calu-3/2B4\n(12hpi, SARS-CoV-2 vs Mock)",
  gaps_col = 3,
  silent = FALSE
)

png(
  "heatmap.png",
  width = 8,
  height = 8,
  units = "in",
  res = 300
)
print(heatmap_plot)  
dev.off()

meta_info <- list(
  experiment = list(
    virus = "SARS-CoV-2",
    cell_line = "Calu-3/2B4", 
    infection_time = "12h", 
    MOI = "1.0",
    replicates = 3,
    control_condition = "Mock"
  ),
  data_source = list(
    accession = "GSE255647", 
    comparison = "SARS-CoV-2 infected vs Mock",
    organism = "Homo sapiens",
    reference = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255647"
  ),
  analysis = list(
    tool = "DESeq2",
    version = paste0("R-", R.version$major, ".", R.version$minor),
    date = format(Sys.Date(), "%Y-%m-%d"),
    parameters = list(
      pvalue_cutoff = 0.05,
      log2fc_cutoff = 1,
      normalization = "rlog")))
output_file <- "meta_analysis_info.json"
write_json(
  meta_info,
  path = output_file,
  pretty = TRUE,
  auto_unbox = TRUE
)
```
