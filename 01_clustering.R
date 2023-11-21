## THIS FILE CLUSTERS THE CELLS USING SEURAT PACKAGE
# official website: https://satijalab.org/seurat/


rm(list = ls())


# 0 LOADING PACKAGES -----------------------------------------------------------
library(Seurat)
library(dplyr)


# 1 EXTRACTING COUNTS ----------------------------------------------------------
# read analyzed data (provided by the author)
VS.analyzed <- readRDS("./data/VS.analyzed.rds")

# split data into 3 patients
VSs <- SplitObject(VS.analyzed, split.by = "orig.ident")

# extract counts
counts.list <- list()
for (i in 1:3) {
  counts <- VSs[[i]][["RNA"]]$counts
  counts.list[i] <- counts
}

#view
counts.list[[1]][1:10, 1:5]

# save
saveRDS(counts.list, "./outputs/counts.list.rds")
counts.list <- readRDS("./outputs/counts.list.rds")


# 2 CREATING SEURAT OBJECT -----------------------------------------------------
# create and merge objects
VS.list <- list()
for (i in 1:3) {
  VS.list[i] <- CreateSeuratObject(counts.list[[i]])
}

# view
VS.list

# save
saveRDS(VS.list, "./outputs/VS.list.rds")
VS.list <- readRDS("./outputs/VS.list.rds")


# 3 NORMALIZATION AND FIND HIGH VARIABLE FEATURES ------------------------------
# use SCTransform
# preprocess and preSCTintegrate
VS.list <- lapply(X = VS.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = VS.list,
                                      nfeatures = 3000)
VS <- PrepSCTIntegration(object.list = VS.list, anchor.features = features)


# 4 INTEGRATING DATA -----------------------------------------------------------
VS.anchors <- FindIntegrationAnchors(
  object.list = VS,
  normalization.method = "SCT",
  anchor.features = features
)
VS.integrated <- IntegrateData(anchorset = VS.anchors,
                               normalization.method = "SCT")

# save
saveRDS(VS.anchors, "./outputs/VS.anchors.rds")
saveRDS(VS.integrated, "./outputs/temp/VS.integrated.rds")
VS.anchors <- readRDS("./outputs/VS.anchors.rds")
VS.integrated <- readRDS("./outputs/VS.integrated.rds")


# 5 PERFORMING AN INTEGRATED ANALYSIS ------------------------------------------
# specify that we will perform downstream analysis on the corrected data note
# that the original unmodified data still resides in the 'RNA' assay
# using SCTransform
# no ScaleData
VS.integrated <- RunPCA(VS.integrated, npcs = 30)
VS.integrated <- RunUMAP(VS.integrated, reduction = "pca", dims = 1:25)
VS.integrated <- RunTSNE(VS.integrated, reduction = "pca", dims = 1:25)
VS.integrated <- FindNeighbors(VS.integrated, reduction = "pca", dims = 1:25)
VS.integrated <- FindClusters(VS.integrated, resolution = 0.5)

# save
saveRDS(VS.integrated, "./outputs/VS.integrated.final.rds")
VS.integrated <- readRDS("./outputs/VS.integrated.final.rds")

# visualize
DimPlot(VS.integrated, label = TRUE, reduction = "umap")
DimPlot(VS.integrated, reduction = "tsne")


# 6 FINDING MARKERS ------------------------------------------------------------
DefaultAssay(VS.integrated) <- "integrated"
all.markers <-
  FindAllMarkers(
    VS.integrated,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

# view markers
# not editting all.markers
all.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

# save all.markers
write.csv(all.markers, "./outputs/markers/all.markers.csv")
all.markers <- read.csv("./outputs/markers/all.markers.csv")

# top 10 markers
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(VS.integrated, features = top10$gene, label = FALSE)
write.csv(top10, "./outputs/markers/top10.csv")
top10 <- read.csv(top10, "./outputs/markers/top10.csv")


# 7 IDENTIFYING CLUSTERS -------------------------------------------------------
DefaultAssay(VS.integrated) <- "SCT"
# according to the paper
# microglia 0, 1, 7, 8, 13, 15
VlnPlot(VS.integrated, features = c("PTPRC", "CD74", "C1QC"))
# SC (I) 2 (, 14)
VlnPlot(VS.integrated, features = c("SOX2", "PMP2", "S100B"))
# SC (II) 14
VlnPlot(VS.integrated, features = c("NCMAP", "MAG", "DRP2"))
# fibroblast 3
VlnPlot(VS.integrated, features = c("DCN", "LUM", "CYP1B1"))
# vascular smooth muscle cell 10
VlnPlot(VS.integrated, features = c("RGS5", "ACTA2", "TAGLN"))
# endothelial 11
VlnPlot(VS.integrated, features = c("PLVAP", "VWF", "FLT1"))
# proliferating microglia 15 (, 18)
VlnPlot(VS.integrated, features = c("TOP2A", "MKI67", "PTPRC"))
# T 6, 18
VlnPlot(VS.integrated, features = c("CCL5", "TRBC2", "CD3E"))
# neutrophil 16
VlnPlot(VS.integrated, features = c("S100A9", "S100A8", "G0S2"))
VlnPlot(VS.analyzed, features = "CMTM2")
# more
VlnPlot(VS.integrated, features = c("CD63", "PTPRC"))

# according to heatmap
# fibroblast
VlnPlot(VS.integrated, features = c("CFH", "ISLR", "COL1A1"))
VlnPlot(VS.integrated, features = c("DCN", "IGFBP5", "LUM"))
VlnPlot(VS.integrated, features = c("APOD", "CYP1B1"))
# microglia
VlnPlot(VS.integrated, features = c("APOC1", "CIQC", "C1QB"))
VlnPlot(VS.integrated, features = c("HLA-DPB1", "HLA-DPA1", "C3"))
VlnPlot(VS.integrated, features = c("CD74", "HLA-DRA"))
# SC (1)
VlnPlot(VS.integrated, features = c("LGI4", "PLP1", "CLU"))
VlnPlot(VS.integrated, features = c("ALDH1A1", "NRXN1", "L1CAM"))
VlnPlot(VS.integrated, features = c("GPM6B", "S100B"))
#
VlnPlot(VS.integrated, features = c("S100A9", "S100A8", "G0S2"))
VlnPlot(VS.integrated, features = c("S100A9", "S100A8", "G0S2"))
VlnPlot(VS.integrated, features = c("S100A9", "S100A8", "G0S2"))
VlnPlot(VS.integrated, features = c("S100A9", "S100A8", "G0S2"))
VlnPlot(VS.integrated, features = c("S100A9", "S100A8", "G0S2"))

# 0 : microglia
# 1 : microglia
# 2 : SC (I)
# 3 : fibroblast
# 4 : SC
# 5 : SC
# 6 : T
# 7 : microglia
# 8 : microglia
# 9 : SC
# 10: vascular smooth muscle cell
# 11: endothelial
# 12: SC
# 13: microglia
# 14: SC (II) myelinated
# 15: proliferating microglia
# 16: neutrophil
# 17: neutrophil
# 18: T

VS.integrated <- RenameIdents(
  VS.integrated,
  `0` = "microglia",
  `1` = "microglia",
  `2` = "SC (I)",
  `3` = "fibroblast",
  `4` = "SC (I)",
  `5` = "SC (I)",
  `6` = "T",
  `7` = "microglia",
  `8` = "microglia",
  `9` = "SC (I)",
  `10` = "vascular smooth muscle cell",
  `11` = "endothelial",
  `12` = "SC (I)",
  `13` = "microglia",
  `14` = "SC (II)",
  `15` = "proliferating microglia",
  `16` = "neutrophil",
  `17` = "neutrophil",
  `18` = "T"
)

DimPlot(VS.integrated,
        label = TRUE,
        repel = TRUE,
        reduction = "umap")

# save
saveRDS(VS.integrated, "./outputs/VS.integrated.labeled.rds")

# check
# mine
table(Idents(VS.integrated))
# authors
table(Idents(VS.analyzed))


# the end
