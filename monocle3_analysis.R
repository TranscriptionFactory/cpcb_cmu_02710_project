# install necessary packages for monocle3
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'),force = TRUE)

devtools::install_version('speedglm', '0.3-4', repos = 'https://packagemanager.rstudio.com/cran/2023-03-31')
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

# load in appropriate libraries
library(monocle3)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)

root <- "~/Desktop/CPCB/Assignments/Spring_2023/genomics/project/analysis/T-cell-development-master/train_input/"
analysis_fold <- "tcell_analysis" # assumes data.E file is stored in some folder off of root
wkdir <- paste0(root,analysis_fold)
setwd(wkdir)

# set analysis variables to match data format
if (analysis_fold == "tcell_analysis"){
  tblfile = "tcell_data.E"
  timecol = "time_point"
  annotcol = "assignment"
} else{
  tblfile = "integrated_data.E"
  timecol = "time"
  annotcol = "bigstage"
}

# read in and modify data file(s) to fit monocle3 conventions
tcell_data <- read.table(tblfile, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
sorted_rows <- order(rownames(tcell_data))
tcell_data <- tcell_data[sorted_rows,]

cell_metadata <- tcell_data[, 1:2]
expr_matrix <- tcell_data[, -c(1, 2)]
expr_matrix <- t(expr_matrix)

if (file.exists("tcell_data.Einit_cluster.txt")){
  cell_cluster_matrix <- read.table("tcell_data.Einit_cluster.txt", sep="\t", header=TRUE, row.names=1)
  sorted_rows <- order(rownames(cell_cluster_matrix))
  cell_cluster_matrix <- cell_cluster_matrix[sorted_rows, , drop = FALSE]
  
  if (!identical(rownames(cell_metadata), rownames(cell_cluster_matrix))) {
    stop("Row names of cell_metadata and cell_cluster_matrix must match.")
  }
  cell_metadata <- cbind(cell_metadata, cell_cluster_matrix)
}
gene_metadata <- data.frame(gene_short_name = rownames(expr_matrix))

if (ncol(expr_matrix) != nrow(cell_metadata)) {
  stop("Number of columns in expression matrix must match number of rows in cell_metadata.")
}
if (nrow(expr_matrix) != nrow(gene_metadata)) {
  stop("Number of rows in expression matrix must match number of rows in gene_metadata.")
}

rownames(cell_metadata) <- colnames(expr_matrix)
rownames(gene_metadata) <- rownames(expr_matrix)

# create object for monocle analysis
cds <- new_cell_data_set(as.matrix(expr_matrix),
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

# preprocess data
cds <- preprocess_cds(cds, num_dim = 100)
# cds <- align_cds(cds, alignment_group = "batch") # no batches to correct :/
plot_pc_variance_explained(cds)

# reduce dimensionality
cds <- reduce_dimension(cds,preprocess_method = "PCA", reduction_method = "UMAP")
plot_cells(cds, color_cells_by=annotcol, label_cell_groups = FALSE)
# plot genes of interest if desired
tcell_genes = c("Cd4","Cd8a","Ctla4","Foxp3")
plot_cells(cds, genes=tcell_genes)
# batch effect check (on timepoints I guess?)
plot_cells(cds, color_cells_by=timecol, label_cell_groups=FALSE)

# cluster cells
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

# learn the trajectory graph
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = annotcol,
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           color_cells_by = timecol,
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

# select graph ordering
cds <- order_cells(cds)

# plot pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

save_monocle_objects(cds=cds, directory_path = "cds_out", comment='First pass on T cell dataset.')
