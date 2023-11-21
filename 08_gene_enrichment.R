## THIS FILE ANALYZES CELL FUNCTIONS USING CLUSTERPROFILER PACKAGE


rm(list = ls())


# 0 LOADING PACKAGES -----------------------------------------------------------
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)


# 1 LOADING DATA ---------------------------------------------------------------
VS.integrated <- readRDS("./outputs/VS.integrated.labeled.rds")


# 2 FINDING MARKERS AND CHANGING DATA INTO REQUIRED FORMAT ---------------------
# find markers
DefaultAssay(VS.integrated) <- "integrated"
SC.markers <- FindMarkers(VS.integrated, ident.1 = "SC (I)",
                           ident.2 = "SC (II)")

up <- rownames(SC.markers[intersect(which(SC.markers [, 1] < 0.05),
                                    which(SC.markers [, 2] >= 0.25)),])
down <- rownames(SC.markers[intersect(which(SC.markers [, 1] < 0.05),
                                      which(SC.markers [, 2] <= (-0.25))),])

# change gene names into ids
SC1.markers <- bitr(
    up,
    fromType = "SYMBOL",
    toType = c("ENTREZID", "ENSEMBL"),
    OrgDb = org.Hs.eg.db
  )

SC2.markers <- bitr(
  down,
  fromType = "SYMBOL",
  toType = c("ENTREZID", "ENSEMBL"),
  OrgDb = org.Hs.eg.db
)

# view
head(SC1.markers)
head(SC2.markers)


# 3 DOING GO ANALYSIS ----------------------------------------------------------
SC1.ego <- enrichGO(
  gene          = SC1.markers$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)
head(summary(SC1.ego))

SC2.ego <- enrichGO(
  gene          = SC2.markers$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)
head(summary(SC2.ego))

# visualize
dotplot(SC1.ego, showCategory=30)
dotplot(SC2.ego, showCategory=30)
barplot(SC1.ego, showCategory=30)
barplot(SC2.ego, showCategory=30)
enrichMap(SC1.ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
cnetplot(SC1.ego)
plotGOgraph(SC1.ego)


















