## THIS FILE PLOTS FIGURE 2


rm(list = ls())


# 0 LOADING PACKAGES -----------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)


# 1 LOADING DATA ---------------------------------------------------------------
SC.integrated <- readRDS("./outputs/SC.integrated.labeled.rds")


# 2 PLOTTING -------------------------------------------------------------------
# FIGURE 2 A TSNE --------------------------------------------------------------
DimPlot(SC.integrated,
        label = TRUE,
        repel = TRUE,
        reduction = "tsne") +
  NoLegend() +
  theme(axis.ticks = element_blank()) + # no x, y axes ticks and text
  theme(axis.text = element_blank()) +
  ggtitle("SC clusters") +
  theme(plot.title = element_text(hjust = 0.5)) # center the title


# FIGURE 2 B VIOLINPLOT --------------------------------------------------------
genes <- c("PMP2", "GFRA3", "VEGFA", "FOSB", "PRX", "MAL")
VlnPlot(SC.integrated,
        features = genes,
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  xlab(NULL) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


# FIGURE 2 C VOCANOPLOT --------------------------------------------------------




# ==============================================================================
# set threshold
log2FC = 2.5
padj = 0.01

# find markers
VS.integrated <- readRDS("./outputs/VS.integrated.labeled.rds")
DefaultAssay(VS.integrated) <- "integrated"

SC.markers <- FindMarkers(
  VS.integrated,
  ident.1 = "SC (I)",
  ident.2 = "SC (II)",
  slot = "data",
  logfc.threshold = 0,
  min.pct = 0
)

# set significancy
SC.markers$threshold = "Normal"
SC.markers[which(SC.markers$avg_log2FC  > log2FC
                  & SC.markers$p_val_adj < padj), ]$threshold = "Up"
SC.markers[which(SC.markers$avg_log2FC  < (-log2FC)
                  & SC.markers$p_val_adj < padj), ]$threshold = "Down"
SC.markers$threshold = factor(SC.markers$threshold,
                               levels = c("Down", "Normal", "Up"))

# plot
ggplot(data = SC.markers, aes(
  x = avg_log2FC,
  y = -log10(p_val_adj),
  color = threshold
)) +
  geom_point(alpha = 0.8, size = 0.8) +
  geom_vline(
    xintercept = c(-log2FC, log2FC),
    linetype = 2,
    color = "grey"
  ) +
  geom_hline(
    yintercept = -log10(padj),
    linetype = 2,
    color = "grey"
  ) +
  xlab(bquote(Log[2] * FoldChange)) +
  ylab(bquote(-Log[10] * italic(P.adj))) +
  theme_classic(base_size = 14) +
  scale_color_manual(
    '',
    labels = c(
      paste0("Down(",
             table(SC.markers$threshold)[[1]], ')'),
      'Normal',
      paste0("Up(",
             table(SC.markers$threshold)[[3]], ')')
    ),
    values = c("blue", "grey", "red")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  xlim(-50, 50)
















# FIGURE 2 D INFERCNV ----------------------------------------------------------
# see in 05_CNVanalysis.R


# FIGURE 2 E CELL CYCLE --------------------------------------------------------
# see in 06_cellcycle.R


# the end

