## THIS FILE PLOTS FIGURE 1


rm(list = ls())


# 0 LOADING PACKAGES -----------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)


# 1 LOADING DATA ---------------------------------------------------------------
VS.integrated <- readRDS("./outputs/VS.integrated.labeled.rds")


# 2 PLOTTING -------------------------------------------------------------------
# FIGURE 1 B UMAP --------------------------------------------------------------
DimPlot(VS.integrated, label = TRUE, reduction = "umap") +
  NoLegend() +
  theme(axis.ticks = element_blank()) + # no x, y axes ticks and text
  theme(axis.text = element_blank()) +
  ggtitle("VS clusters") +
  theme(plot.title = element_text(hjust = 0.5)) # center the title


# FIGURE 1 C HEATMAP -----------------------------------------------------------
plotting <-
  c(
    "APOC1",
    "C1QC",
    "C1QB",
    "HLA-DPB1",
    "HLA-DPA1",
    "C3",
    "CD74",
    "HLA-DRA",
    "LGI4",
    "PLP1",
    "CLU",
    "ALDH1A1",
    "NRXN1",
    "L1CAM",
    "GPM6B",
    "S100B",
    "CFH",
    "ISLR",
    "COL1A1",
    "DCN",
    "IGFBP5",
    "LUM",
    "APOD",
    "CYP1B1",
    "GZMA",
    "TRAC",
    "CD52",
    "CD2",
    "IL7R",
    "IL32",
    "TRBC2",
    "CCL5",
    "CALD1",
    "IGFBP7",
    "BGN",
    "NDUFA4L2",
    "NOTCH3",
    "RGS5",
    "TAGLN",
    "ACTA2",
    "TM4SF1",
    "CD34",
    "ADGRF5",
    "VWF",
    "IGFBP3",
    "FLT1",
    "ITM2A",
    "PLVAP",
    "CXCL14",
    "MME",
    "SFRP1",
    "GLDN",
    "MLIP",
    "NCMAP",
    "FGFBP2",
    "PRX",
    "HMMR",
    "TYMS",
    "NUSAP1",
    "TPX2",
    "ASPM",
    "CENPF",
    "MKI67",
    "TOP2A",
    "S100A12",
    "NCF1",
    "RIPOR2",
    "CXCR2",
    "G0S2",
    "FCGR3B",
    "S100A8",
    "S100A9"
  )

DefaultAssay(VS.integrated) <- "integrated"
DoHeatmap(
  VS.integrated,
  slot = "scale.data",
  features = plotting,
  disp.max = 1.5,
  disp.min = -1.5,
  label = FALSE
)
## no CMTM2 in high variable features


# FIGURE 1 D FEATURE PLOT ------------------------------------------------------
DefaultAssay(VS.integrated) <- "SCT"
p <- FeaturePlot(
  VS.integrated,
  features = c(
    "C1QC",
    "S100B",
    "CCL5",
    "S100A9",
    "RGS5",
    "LUM",
    "PLVAP",
    "MKI67"
  ),
  cols = c("grey", "red"),
  reduction = "umap",
  ncol = 2
)
for (i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
p


# FIGURE 1 E BARPLOT -----------------------------------------------------------
# calculate proportion
cell.prop <- as.data.frame(prop.table(table(
  Idents(VS.integrated),
  VS.integrated$orig.ident
)))
colnames(cell.prop) <- c("clusters", "samples", "proportion")
cell.prop$samples <- c(rep("VS1", 9), rep("VS2", 9), rep("VS3", 9))

# plot
ggplot(cell.prop, aes(samples, proportion, fill = clusters)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("") +
  theme_light() +
  ylab(NULL) +
  scale_fill_brewer(palette = "RdBu")


# the end
