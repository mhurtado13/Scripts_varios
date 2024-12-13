#' Inference of TFs activity based on gene expression 
#' 
#' \code{compute.TFs.activity} computes TFs activty based on a gene expression matrix using VIPER algorithm and a collection of GRN (collectri) from the OmnipathR package  
#' 
#' @param RNA.tpm Gene expression matrix normalized by TPM (genes X samples).
#' @return A matrix of protein activity (samples X tfs).
#' 
compute.TFs.activity <- function(RNA.counts, TF.collection = "CollecTRI", min_targets_size = 3, tfs.pruned = FALSE, universe){
  
  tfs2viper_regulons <- function(df){
    regulon_list <- split(df, df$source)
    regulons <- lapply(regulon_list, function(regulon) {
      tfmode <- stats::setNames(regulon$mor, regulon$target)
      list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
    })
    return(regulons)}
  
  if(tfs.pruned==T){
    cat("Pruned TFs is set to TRUE. Please specify the maximun size of targets allowed/n")
    max_size_targets = as.numeric(readline(prompt = "Maximun size of TFs-targets: "))
  }
  
  if(TF.collection == "CollecTRI"){
    net_regulons = tfs2viper_regulons(universe)
  } else if(TF.collection == "Dorothea"){
    net = decoupleR::get_dorothea(organism = 'human', levels = c("A", "B", "C", "D"))
    net_regulons = tfs2viper_regulons(net)
  } 
  
  if(TF.collection == "ARACNE"){
    cat("For ARACNE analysis you need to specify the path of your network file. Remember this file should be a 3 columns text file, with regulator in the first column, target in the second and mutual information in the third column")
    network_file = readline(prompt = "Path for network file from aracne (no quotes): ")
    net_regulons <- aracne2regulon(network_file, as.matrix(RNA.counts), format = "3col")
  }
  
  if(tfs.pruned == TRUE){
    net_regulons = pruneRegulon(net_regulons, cutoff = max_size_targets)
  }
  
  sample_acts <- viper::viper(as.matrix(RNA.counts), net_regulons, minsize = min_targets_size, verbose=F, method = "scale")
  message("TFs scores computed")
  
  return(data.frame(t(sample_acts)))
  
}

create_tfs_modules = function(TF.matrix, network_tfs){
  library(dplyr)
  tfs.modules = TF.matrix %>%
    t() %>%
    data.frame() %>%
    dplyr::mutate(Module = "na")
  
  for (i in 1:length(network_tfs[[3]])) {
    tfs.modules$Module[which(rownames(tfs.modules) %in% network_tfs[[3]][[i]])] = names(network_tfs[[3]])[i]
  }
  
  tfs_colors = tfs.modules %>%
    dplyr::pull(Module)
  
  MEList = WGCNA::moduleEigengenes(TF.matrix, colors = tfs_colors, scale = F) #Data already scale
  MEs = MEList$eigengenes
  MEs =  WGCNA::orderMEs(MEs)
  
  return(MEs)
}

minMax <- function(x) {
  #columns: features
  x = data.matrix(x)
  for(i in 1:ncol(x)){
    x[,i] = (x[,i] - min(x[,i], na.rm = T)) / (max(x[,i], na.rm = T) - min(x[,i], na.rm = T))
  }
  
  return(x)
  
}