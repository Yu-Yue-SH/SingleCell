## THIS FILE DOES PSUDOTIME ANALYSIS USING MONOCLE3 PACKAGE
# official website: https://cole-trapnell-lab.github.io/monocle3/


rm(list = ls())


# 0 LOADING PACKAGES -----------------------------------------------------------
library(monocle3)
library(Seurat)
# library(SeuratWrappers)
# library(dplyr)
# library(ggplot2)


# 1 LOADING DATA ---------------------------------------------------------------
# load seurat object
SC.integrated <- readRDS("./outputs/SC.integrated.labeled.rds")
SC.integrated[["cell_type"]] <- Idents(SC.integrated)

# create cell_data_set
counts <- as.matrix(SC.integrated[["SCT"]]$counts)
meta.data <- SC.integrated@meta.data
gene.data <- data.frame(gene_short_name = row.names(counts),
                        row.names = row.names(counts))

SC.cells <- new_cell_data_set(counts, cell_metadata = meta.data,
                              gene_metadata = gene.data)


# 2 PREPROCESS -----------------------------------------------------------------
SC.cells <- preprocess_cds(SC.cells, num_dim = 50)
plot_pc_variance_explained(SC.cells)
SC.cells <- align_cds(SC.cells, alignment_group = "orig.ident")


# 3 REDUCING DIMENTIONALITY ----------------------------------------------------
SC.cells <- reduce_dimension(SC.cells)
plot_cells(
  SC.cells,
  label_groups_by_cluster = FALSE,
  group_label_size = 4,
  color_cells_by = "cell_type"
)


# 4 CLUSTERING CELLS -----------------------------------------------------------
SC.cells <- cluster_cells(SC.cells, resolution=1e-5)
plot_cells(SC.cells)
plot_cells(cds, color_cells_by = "partition", group_cells_by = "partition")


# 5 TRAJECTORY ANALYSIS --------------------------------------------------------
SC.cells <- learn_graph(SC.cells)
plot_cells(
  SC.cells,
  color_cells_by = "cell_type",
  label_groups_by_cluster = TRUE,
  group_label_size = 4,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

# order the cells in pseudotime
plot_cells(SC.cells,
           color_cells_by = "cell_type",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=3)




























