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

 head(colnames(raw_data_human), 10)
 [1] "Series1_NHBE_Mock_1"       "Series1_NHBE_Mock_2"
 [3] "Series1_NHBE_Mock_3"       "Series1_NHBE_SARS-CoV-2_1"
 [5] "Series1_NHBE_SARS-CoV-2_2" "Series1_NHBE_SARS-CoV-2_3"
 [7] "Series2_A549_Mock_1"       "Series2_A549_Mock_2"
 [9] "Series2_A549_Mock_3"       "Series2_A549_SARS-CoV-2_1"

以SARS-CoV-1感染Calu-3/2B4 1 MOI 24h为例：
mock_samples <- grep("^MOCK_24hpi_rep[123]$", colnames(raw_data_all), value=TRUE)
virus_samples <- grep("^SARS[-.]1_24hpi_rep[123]$", colnames(raw_data_all), value=TRUE)
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
  timepoint = "24hpi",
  virus_type = c(rep(NA, 3), rep("SARS-1", 3)),
  stringsAsFactors = FALSE
)
setwd("/mnt/alamo01/users/chenyun730/sampledata/sars_cov_1/Calu-3_2B4")
write.csv(count_matrix,"12H_MOI1_GSE255647.csv",row.names=FALSE)




