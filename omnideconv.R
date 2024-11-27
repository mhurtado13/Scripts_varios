library(Seurat)
library(tidyverse)
library(omnideconv)
library(readxl)

load('Vanderbilt_new_annot.RData')
DimPlot(seurat_anndata)
rna = read.csv('Documents/LungPredict2/Data/Vanderbilt/Raw_files/TPM_Counts_Van.txt', sep = '\t', row.names = 1)

counts.matrix <- as.matrix(seurat_anndata@assays$RNA@counts)
cell.type.annotations <- seurat_anndata$new_annotation
batch.ids = seurat_anndata$sample

########### DWLS
signature.matrix.dwls <- omnideconv::build_model(single_cell_object = counts.matrix,
                                                 cell_type_annotations = cell.type.annotations,
                                                 method = 'dwls', 
                                                 dwls_method = 'mast_optimized')

deconvolution.results.bayesprism <- deconvolute(bulk_gene_expression = rna,
                                                single_cell_object = counts.matrix,
                                                cell_type_annotations = cell.type.annotations,
                                                signature=NULL,
                                                method = 'bayesprism', 
                                                n_cores=12)

omnideconv::plot_deconvolution(list('bayesprism' = deconvolution.results.bayesprism), "bar", "method", "Spectral")

set_cibersortx_credentials("marcelo.hurtado@inserm.fr", "734212f6ad77fc4eea2bdb502792f294")
signature.matrix.cbsx = omnideconv::build_model_cibersortx(single_cell_object = counts.matrix,
                                                           cell_type_annotations = cell.type.annotations)
deconv_cbsx = deconvolute_cibersortx(rna, signature.matrix.cbsx)
omnideconv::plot_deconvolution(list('CBSX' = deconv_cbsx), "heatmap", "method", "Spectral")
