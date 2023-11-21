## THIS FILE ANALYZES CNV USING INFERCNV PACKAGE
# official website: https://github.com/broadinstitute/infercnv


rm(list = ls())


# 0 LOADING PACKAGES -----------------------------------------------------------
library(Seurat)
library(AnnoProbe)
library(gtools)


# 1 LOADING DATA ---------------------------------------------------------------
VS.integrated <- readRDS("./outputs/VS.integrated.labeled.rds")


# 2 PREPARING FILES ------------------------------------------------------------
# get counts
counts <-
  GetAssayData(object = VS.integrated,
               layer = "counts",
               assay = "SCT")

# get annotation
cell.annotation <- as.data.frame(Idents(VS.integrated))

# get gene information
gene.infor <- annoGene(rownames(counts), "SYMBOL")
gene.infor <- gene.infor[!(gene.infor$chr %in% c("chrM", "chrX", "chrY")), ]
gene.infor <- gene.infor[with(gene.infor, order(chr, start)), c(1, 4:6)]
gene.infor <- gene.infor[with(gene.infor, mixedorder(chr)), ]
gene.infor <- gene.infor[!duplicated(gene.infor[, 1]), ]

# match counts
counts <- counts[rownames(counts) %in% gene.infor[, 1], ]
counts <- counts[match(gene.infor[, 1], rownames(counts)), ]

# view
counts[1:10, 1:5]
head(cell.annotation)
head(gene.infor)

# save
counts.file <- "./outputs/infercnv/counts.txt"
annotation.file <- "./outputs/infercnv/cell.annotation.txt"
gene.infor.file <- "./outputs/infercnv/gene.infor.txt"
write.table(counts, counts.file, sep = "\t", quote = FALSE)
write.table(
  cell.annotation,
  annotation.file,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)
write.table(
  gene.infor,
  gene.infor.file,
  sep = '\t',
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)


# 3 DOING INFERCNV -------------------------------------------------------------
# create object
counts.file <- "./outputs/infercnv/counts.txt"
annotation.file <- "./outputs/infercnv/cell.annotation.txt"
gene.infor.file <- "./outputs/infercnv/gene.infor.txt"
normal.cells <-
  c(
    "microglia",
    "fibroblast",
    "T",
    "vascular smooth muscle cell",
    "endothelial",
    "SC (II)",
    "proliferating microglia",
    "neutrophil"
  )
SC.cnv <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = counts.file,
  annotations_file = annotation.file,
  ref_group_names = normal.cells,
  gene_order_file = gene.infor.file,
  delim = "\t"
)

# do infercnv
out_dir <- "./outputs/infercnv/cnv_files"
options("Seurat.object.assay.version" = "v3")
SC.cnv.default <- infercnv::run(
  SC.cnv,
  cutoff = 0.1,
  out_dir = out_dir,
  cluster_by_groups = TRUE,
  plot_steps = FALSE,
  denoise = TRUE,
  HMM = FALSE,
  # TRUE or FALSE
  no_prelim_plot = TRUE,
  leiden_resolution = 0.00001,
  png_res = 360
)


# the end
