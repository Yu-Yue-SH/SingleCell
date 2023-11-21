## THIS FILE INSTALLS ALL THE PACKAGES NEEDED FOR THE ANALYSIS


# Monocle3 and Seurat (is also installed)
install.packages("BiocManager")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
# test
library(monocle3)

# 
packages.install("AnnoProbe")
BiocManager::install("infercnv")