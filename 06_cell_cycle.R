## THIS FILE ANALYZES CELL CYCLE USING SEURAT PACKAGE


rm(list = ls())


# 0 LOADING PACKAGES -----------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)


# 1 LOADING DATA ---------------------------------------------------------------
SC.integrated <- readRDS("./outputs/SC.integrated.labeled.rds")


# 2 CELL CYCLE SCORING ---------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

SC.integrated <- CellCycleScoring(
  SC.integrated,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)


# 3 PLOTTING -------------------------------------------------------------------
phase <- as.data.frame(
    table(SC.integrated@meta.data[["Phase"]], Idents(SC.integrated)))
colnames(phase) <- c("phase", "cluster", "num")

ggplot(phase, aes(x = cluster, y = num, fill = phase)) +
  geom_bar(stat = 'identity', position = "fill") +
  scale_fill_brewer(palette = "Spectral")


# the end
