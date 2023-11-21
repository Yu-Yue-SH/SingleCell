## THIS FILE CLUSTERS THE SUBSET OF SC CELLS USING SEURAT PACKAGE


rm(list = ls())


# 0 LOADING PACKAGES -----------------------------------------------------------
library(Seurat)


# 1 LOADING DATA ---------------------------------------------------------------
VS.integrated <- readRDS("./outputs/VS.integrated.labeled.rds")


# 2 SUBSETTING SC GROUPS -------------------------------------------------------
SC.subset <- VS.integrated[, Idents(VS.integrated) %in% c("SC (I)",  "SC (II)")]
# view
SC.subset


# 3 RECLUSTERING SC.SUBSET -----------------------------------------------------
SC.list <- SplitObject(SC.subset, split.by = "orig.ident")
SC.list <- lapply(X = SC.list, FUN = SCTransform)
SC.features <- SelectIntegrationFeatures(object.list = SC.list,
                                         nfeatures = 3000)
SC.list <- PrepSCTIntegration(object.list = SC.list,
                              anchor.features = SC.features)
SC.anchors <- FindIntegrationAnchors(
  object.list = SC.list,
  normalization.method = "SCT",
  anchor.features = SC.features
)
SC.integrated <- IntegrateData(anchorset = SC.anchors,
                               normalization.method = "SCT")

# save
saveRDS(SC.integrated, "./outputs/SC.integrated.rds")
SC.integrated <- readRDS("./outputs/SC.integrated.rds")

# identify clusters
DefaultAssay(SC.integrated) <- "integrated"
ElbowPlot(SC.integrated)
SC.integrated <- RunPCA(SC.integrated, npcs = 30)
SC.integrated <- RunUMAP(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- RunTSNE(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindNeighbors(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindClusters(SC.integrated, resolution = 0.3)

# visualize
DimPlot(SC.integrated, label = TRUE, reduction = "umap")
DimPlot(SC.integrated, label = TRUE, reduction = "tsne")

DefaultAssay(SC.integrated) <- "SCT"
# all
VlnPlot(SC.integrated, features = c("MAL"))
# GFRA3+ 0, 2, 3, 6, 8
VlnPlot(SC.integrated, features = c("GFRA3"))
# PMP2+ 1
VlnPlot(SC.integrated, features = c("PMP2"))
# FOSB+ 5
VlnPlot(SC.integrated, features = c("FOSB"))
# VEGFA+ 4
VlnPlot(SC.integrated, features = c("VEGFA"))
# PRX+ 7
VlnPlot(SC.integrated, features = c("PRX"))

SC.integrated <-
  RenameIdents(
    SC.integrated,
    `0` = "GFRA3+",
    `1` = "PMP2+",
    `2` = "GFRA3+",
    `3` = "GFRA3+",
    `4` = "VEGFA+",
    `5` = "FOSB+",
    `6` = "GFRA3+",
    `7` = "PRX+",
    `8` = "GFRA3+"
  )

# visualize
DimPlot(SC.integrated, label = TRUE, reduction = "umap")
DimPlot(SC.integrated, label = TRUE, reduction = "tsne")

# save
saveRDS(SC.integrated, "./outputs/SC.integrated.labeled.rds")
SC.integrated <- readRDS("./outputs/SC.integrated.labeled.rds")


# the end
