benchmarking_deconvolution = function(deconvolution, groundtruth, include_combinations = NULL, corr_type = "spearman", scatter = T, file_name){
  require(reshape2)
  require(ggpubr)
  require(dplyr)
  require(purrr)
  require(tibble)
  #deconvolution = deconv
  #groundtruth = cells_groundtruth
  
  deconvolution_combinations =  c("Quantiseq", "MCP", "XCell", "Epidish_BPRNACan_",  "Epidish_BPRNACanProMet", "Epidish_BPRNACan3DProMet", "Epidish_CBSX.HNSCC.scRNAseq", "Epidish_CBSX.Melanoma.scRNAseq",
                                  "Epidish_CBSX.NSCLC.PBMCs.scRNAseq", "Epidish_CCLE.TIL10", "Epidish_TIL10", "Epidish_LM22", "DeconRNASeq_BPRNACan_", "DeconRNASeq_BPRNACanProMet",
                                  "DeconRNASeq_BPRNACan3DProMet", "DeconRNASeq_CBSX.HNSCC.scRNAseq", "DeconRNASeq_CBSX.Melanoma.scRNAseq", "DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq", "DeconRNASeq_CCLE.TIL10",
                                  "DeconRNASeq_TIL10", "DeconRNASeq_LM22", "CBSX_BPRNACan_", "CBSX_BPRNACanProMet", "CBSX_BPRNACan3DProMet", "CBSX_CBSX.HNSCC.scRNAseq",
                                  "CBSX_CBSX.Melanoma.scRNAseq", "CBSX_CBSX.NSCLC.PBMCs.scRNAseq", "CBSX_CCLE.TIL10", "CBSX_TIL10", "CBSX_LM22", "DWLS_BPRNACan_",  "DWLS_BPRNACanProMet", "DWLS_BPRNACan3DProMet", "DWLS_CBSX.HNSCC.scRNAseq", "Epidish_CBSX.Melanoma.scRNAseq",
                                  "DWLS_CBSX.NSCLC.PBMCs.scRNAseq", "DWLS_CCLE.TIL10", "DWLS_TIL10", "DWLS_LM22")
  
  if(is.null(include_combinations)!=T){
    deconvolution_combinations = c(deconvolution_combinations, include_combinations)
  }
  

  ###Correlation function
  correlation <- function(data, corr) {
    require(tidyr)
    M <- Hmisc::rcorr(as.matrix(data), type = corr)
    Mdf <- map(M, ~data.frame(.x))
    
    corr_df = Mdf %>%
      map(~rownames_to_column(.x, var="measure1")) %>%
      map(~pivot_longer(.x, -measure1, names_to = "measure2")) %>%
      bind_rows(.id = "id") %>%
      pivot_wider(names_from = id, values_from = value) %>%
      mutate(sig_p = ifelse(P < .05, T, F),
             p_if_sig = ifelse(sig_p, P, NA),
             r_if_sig = ifelse(sig_p, r, NA)) 
    
    return(corr_df)
    
  }
  
  #####Scatter plot function
  scatter_plots = function(deconv, ground, method){
    for (i in 1:ncol(deconv)) {
      data = cbind(deconv[,i], ground)
      colnames(data) = c("x", "y")
      p <- ggplot(data, aes(x = x, y = y)) +
        geom_point(color = "blue", size = 3, alpha = 0.7) +  # Customize the points
        geom_smooth(method = "lm", se = T, color = "red") +  # Add regression line
        theme_minimal() +  # Apply a minimal theme
        labs(
          #title = paste0("Linear correlation - ", colnames(ground)),  # Set the title
          x = colnames(deconv)[i],                 # Set the x-axis label
          y = colnames(ground)                # Set the y-axis label
        ) +
        theme(
          plot.title = element_text(size = 20, hjust = 0.5),  # Adjust title font size and position
          axis.title = element_text(size = 15)                # Adjust axis label font size
        ) +
        stat_cor(method = method,label.x = summary(data[,1])[[5]]-0.02, label.y = summary(data[,2])[[5]]+0.07)  # Add correlation coefficient
      
      print(p)
    }
  }
  
  cell_clusters = colnames(groundtruth)
  
  ###Correlation matrix
  corr_matrix = data.frame(matrix(ncol = length(deconvolution_combinations), nrow = length(cell_clusters)))
  rownames(corr_matrix) = cell_clusters
  colnames(corr_matrix) = deconvolution_combinations
  cells_discard = c()
  
  ###Correlation computation
  for (i in 1:length(cell_clusters)) {
    #Consider specific cases (Mast cells including Mast activated and resting, NK including activated and resting and Macrophages including M1 and M2)
    if(cell_clusters[i]=="Mast.cells"){
      deconv = deconvolution[,grep(paste0(c(cell_clusters[i], "Mast.activated.cells", "Mast.resting.cells"), collapse = "|"), colnames(deconvolution))]
      Mast.active = deconv[,grep("Mast.activated.cells", colnames(deconv))]
      Mast.rest = deconv[,grep("Mast.resting.cells", colnames(deconv))]
      Mast_cells = Mast.active + Mast.rest
      methods_signatures = unlist(strsplit(colnames(Mast_cells), paste0(c("Mast.activated.cells", "Mast.resting.cells"), collapse = "|")))
      colnames(Mast_cells) = paste0(methods_signatures, "Mast.cells")
      deconv = deconv[,-which(colnames(deconv)%in%c(colnames(Mast.active), colnames(Mast.rest)))]
      deconv = cbind(deconv, Mast_cells)
    }else if(cell_clusters[i]=="NK.cells"){
      deconv = deconvolution[,grep(paste0(c(cell_clusters[i], "NK.activated", "NK.resting"), collapse = "|"), colnames(deconvolution))]
      NK.active = deconv[,grep("NK.activated", colnames(deconv))]
      NK.rest = deconv[,grep("NK.resting", colnames(deconv))]
      NK_cells = NK.active + NK.rest
      methods_signatures = unlist(strsplit(colnames(NK_cells), paste0(c("NK.activated", "NK.resting"), collapse = "|")))
      colnames(NK_cells) = paste0(methods_signatures, "NK.cells")
      deconv = deconv[,-which(colnames(deconv)%in%c(colnames(NK.active), colnames(NK.rest)))]
      deconv = cbind(deconv, NK_cells)
    }else if(cell_clusters[i]=="Macrophages.cells"){
      deconv = deconvolution[,grep(paste0(c(cell_clusters[i], "Macrophages.M0" ,"Macrophages.M1", "Macrophages.M2"), collapse = "|"), colnames(deconvolution))]
      m0 = deconv[,grep("Macrophages.M0", colnames(deconv))]
      m1 = deconv[,grep("Macrophages.M1", colnames(deconv))]
      m2 = deconv[,grep("Macrophages.M2", colnames(deconv))]
      m1_m2 = m1 + m2
      methods_m0 = unlist(strsplit(colnames(m0), "Macrophages.M0")) #Extract methods-signatures containing M0
      macrophages_m1_m2 = m1_m2[,-grep(paste0(methods_m0, collapse = "|"), colnames(m1_m2))] #Subset methods-signatures containing only M1 and M2
      macrophages_m0_m1_m2 = m1_m2[,grep(paste0(methods_m0, collapse = "|"), colnames(m1_m2))] #Extract methods-signatures containing M0, M1, M2
      m0_m1_m2 = macrophages_m0_m1_m2 + m0 #Methods with m0 + m1 + m2
      macrophages_all = cbind(macrophages_m1_m2, m0_m1_m2)
      methods_signatures = unlist(strsplit(colnames(macrophages_all), paste0(c("Macrophages.M0", "Macrophages.M1", "Macrophages.M2"), collapse = "|")))
      colnames(macrophages_all) = paste0(methods_signatures, "Macrophages.cells")
      deconv = deconv[,-which(colnames(deconv)%in%c(colnames(m0), colnames(m1), colnames(m2)))]
      deconv = cbind(deconv, macrophages_all)
    }else if(cell_clusters[i]%in%"Macrophages_no_defined_M1"){
      deconv = deconvolution[,grep("Macrophages.M1", colnames(deconvolution))]
    }else if(cell_clusters[i]%in%"Macrophages_no_defined_M2"){
      deconv = deconvolution[,grep("Macrophages.M2", colnames(deconvolution))]
    }else if(cell_clusters[i]%in%"Macrophages_no_defined_M0"){
      deconv = deconvolution[,grep("Macrophages.M0", colnames(deconvolution))]
    }else if(cell_clusters[i]%in%"Macrophages_no_defined_Monocytes"){
      deconv = deconvolution[,grep("Monocytes", colnames(deconvolution))]
    }else if(cell_clusters[i]%in%"Macrophages_M2_M1"){
      deconv = deconvolution[,grep("Macrophages.M1", colnames(deconvolution))]
    }else if(cell_clusters[i]%in%"Macrophages_M2_M2"){
      deconv = deconvolution[,grep("Macrophages.M2", colnames(deconvolution))]
    }else if(cell_clusters[i]%in%"Macrophages_M2_M0"){
      deconv = deconvolution[,grep("Macrophages.M0", colnames(deconvolution))]
    }else if(cell_clusters[i]%in%"Macrophages_M2_Monocytes"){
      deconv = deconvolution[,grep("Monocytes", colnames(deconvolution))]
    }else{
      idx = grep(cell_clusters[i], colnames(deconvolution))
      if(length(idx)==0){
          cells_discard = c(cells_discard, cell_clusters[i])
      }
        
      deconv = deconvolution[,idx] 
    }
    
    ground = groundtruth[,cell_clusters[i],drop=F]
    
    ###Scatter plots
    if(scatter == T){
      if(ncol(deconv)!=0){
        pdf(paste0("Scatter_plots_", colnames(ground)))
        scatter_plots(deconv, ground, corr_type)
        dev.off()
      }
    }
    
    x = correlation(cbind(deconv, ground), corr_type)
    x = x[which(x$measure1==colnames(ground)),] #only taking corr against ground truth
      
    for (j in 1:ncol(corr_matrix)) {
      idx = grep(colnames(corr_matrix)[j], x$measure2)
      if(length(idx) == 0){
        corr_matrix[i,j] = NA
      }else{
        corr_matrix[i,j] = x$r[idx]
      }
    }
  }


  ###Benchmarking plot
  if(length(cells_discard)>0){
    corr_matrix = corr_matrix[-which(rownames(corr_matrix)%in%cells_discard),]
  }else{
    corr_matrix = corr_matrix
  }
  corr_matrix[nrow(corr_matrix)+1,] = colMeans(corr_matrix, na.rm = T)
  rownames(corr_matrix)[nrow(corr_matrix)] = "average"
  
  ##Order methods
  corr_matrix = t(corr_matrix) %>%
    data.frame() %>%
    arrange(average) %>%
    t() %>%
    data.frame()
  
  corr_df <- melt(corr_matrix)
  corr_df = corr_df %>%
    mutate(Cells = rep(rownames(corr_matrix), ncol(corr_matrix)))
  
  g <- corr_df %>%
    ggplot(aes(Cells, variable, fill=value, label=round(value,2))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = paste0(corr_type, "'s\nCorrelation"), title=file_name) +
    scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
    geom_text() +
    theme_classic() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggpubr::rotate_x_text(angle = 45) + theme(axis.text.x=ggtext::element_markdown()) + theme(axis.text.y=ggtext::element_markdown())
  
  pdf(paste0("Benchmark_plot_", file_name), width = 14)
  print(g)
  dev.off()
  
  return(corr_matrix)
  
}
