
### CellTFusion functions for Emilie DEG analysis

#' Run differential expression analysis with edgeR/limma-voom
#'
#' Filters low-expression genes, applies TMM normalization, runs voom
#' transformation, fits a linear model, and returns the top differentially
#' expressed genes via \code{limma::topTable}.
#'
#' @param counts A raw count matrix (genes x samples).
#' @param coldata A data frame of sample metadata whose row names match
#'   the column names of \code{counts}.
#' @param group_col Character. Name of the column in \code{coldata} used as
#'   the grouping factor for differential expression.
#' @param ref_level Character or \code{NULL}. Reference level for the group
#'   factor. If \code{NULL}, the default factor ordering is used.
#'
#' @return A data frame of differentially expressed genes (p.adj < 0.05) as
#'   returned by \code{limma::topTable}, with columns \code{logFC},
#'   \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val}, and \code{B}.
#'
#' @keywords internal
run_deg_analysis <- function(counts, coldata, group_col, ref_level = NULL, pval = 0.05) {
  # Prepare counts
  counts_mat <- as.matrix(counts)
  mode(counts_mat) <- "numeric"
  counts_mat <- counts_mat[, rownames(coldata)]

  # Create group factor
  group <- factor(coldata[[group_col]])

  # Set reference level if provided
  if (!is.null(ref_level)) {
    group <- stats::relevel(group, ref = ref_level)
  }

  # Create DGE object and filter
  dge <- edgeR::DGEList(counts = counts_mat, group = group)
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)

  # Design matrix
  design <- stats::model.matrix(~ group)

  # voom transformation
  v <- limma::voom(dge, design)

  # Fit model
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  # Extract coefficient name (second column of design)
  coef_name <- colnames(design)[2]

  # Get results
  res <- limma::topTable(fit, coef = coef_name, p.value = pval, number = Inf)

  return(res)
}


#' Compute Transcription Factor (TF) activity
#'
#' Infers transcription factor (TF) activity from a gene expression matrix using the VIPER algorithm (Alvarez et al., 2016). The function requires a TF-target gene regulatory network, which can be provided by the user or obtained from OmnipathR resources such as CollecTRI or Dorothea. ARACNE-inferred networks are also supported.
#'
#' @param RNA.counts A gene expression matrix with genes as rows and samples as columns. The matrix should be normalized (e.g., TPM, log2CPM, etc.).
#' @param TF.collection Character. The source of the TF-target network. Options are `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.
#' - `"CollecTRI"` and `"Dorothea"` use prebuilt collections from OmnipathR.
#' - `"ARACNE"` allows user input of a custom network file in a 3-column format: `regulator`, `target`, and `mutual information`.
#' @param min_targets_size Integer. Minimum number of target genes per regulon required for TF activity inference. Default is 5.
#' @param universe Optional. A user-specified data frame of TF-target interactions. If not provided, the function will fetch the relevant network based on the `TF.collection` argument.
#' @param cancer.type Optional character. Cancer type label used when caching the TF collection.
#' @param cores Integer. Number of cores used by VIPER inference. Default is 4.
#' @param scale Logical. If TRUE (default), z-score scales the TF activity matrix across samples.
#' @param return Logical; if TRUE, saves matrix in Results/ folder. Default is TRUE.
#' @param file.name Optional character suffix used when writing the TF activity matrix to disk.
#'
#' @return A data frame of inferred and scaled TF activity scores, with samples as rows and TFs as columns.
#'
#' @references
#' Alvarez, M. et al. (2016). Functional characterization of somatic mutations in cancer using network-based inference of protein activity. *Nature Genetics*, 48(8), 838-847. https://doi.org/10.1038/ng.3593
#'
#' Tuerei, D., Korcsmaros, T., & Saez-Rodriguez, J. (2016). OmniPath: guidelines and gateway for literature-curated signaling pathway resources. *Nature Methods*, 13(12), 966-967. https://doi.org/10.1038/nmeth.4077
#'
#' Garcia-Alonso, L. et al. (2019). Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Research*. https://doi.org/10.1101/gr.240663.118
#'
#' Lachmann, A. et al. (2016). ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information. *Bioinformatics*, 32(14), 2233-2235. https://doi.org/10.1093/bioinformatics/btw216
#'
#' Margolin, A.A. et al. (2006). ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. *BMC Bioinformatics*, 7(Suppl 1), S7. https://doi.org/10.1186/1471-2105-7-S1-S7
#'
#' @examples
#' data("counts.norm.tuto")
#' tfs_activity <- compute.TFs.activity(counts.norm.tuto, cores = 1)
#'
compute.TFs.activity <- function(RNA.counts, TF.collection = "CollecTRI", min_targets_size = 5, universe = NULL, cancer.type = NULL, cores = 3, return = TRUE, file.name = NULL){

  tf_cache_file <- "Results/TF_target_collection.csv"

  if(TF.collection == "ARACNE"){

    if(is.null(cancer.type)){
      # auto-discover when only one network exists
      candidates <- list.files("input/ARACNE", pattern = "^network\\.txt$",
                               recursive = TRUE, full.names = TRUE)
      if(length(candidates) == 0)
        stop("TF.collection = 'ARACNE' requires a 'cancer.type' or a network.txt under input/ARACNE/")
      if(length(candidates) > 1)
        stop("Multiple ARACNe networks found. Specify 'cancer.type' (e.g. cancer.type = 'skcm'):\n",
             paste(dirname(dirname(candidates)), collapse = "\n"))
      aracne.network <- candidates[1]
      cat("Auto-detected ARACNe network:", aracne.network, "\n")
    } else {
      aracne.network <- file.path("~/Documents/CellTFusion_paper/input/ARACNE", cancer.type, "network/network.txt")
      if(!file.exists(aracne.network))
        stop("ARACNe network not found for cancer type '", cancer.type, "': ", aracne.network)
    }

    cat("Loading ARACNe network from:", aracne.network, "\n")
    # Read network edges, filter to genes present in expression matrix
    aracne_net <- utils::read.table(aracne.network, header = TRUE, sep = "\t") %>%
      dplyr::select(source = Regulator, target = Target) %>%
      dplyr::filter(source %in% rownames(RNA.counts) & target %in% rownames(RNA.counts))

    # Compute Spearman correlation between every TF and its targets in one matrix op
    # (this is exactly what TFmode1 does internally)
    all_tfs     <- unique(aracne_net$source)
    all_targets <- unique(aracne_net$target)
    cor_mat <- suppressWarnings(
      stats::cor(t(RNA.counts[all_tfs, , drop = FALSE]),
                 t(RNA.counts[all_targets, , drop = FALSE]),
                 method = "spearman")
    )

    # mor = sign of Spearman correlation - +1 activation, -1 repression
    universe <- aracne_net %>%
      dplyr::mutate(mor = cor_mat[cbind(source, target)]) %>%
      dplyr::mutate(mor = sign(mor)) %>%
      dplyr::filter(!is.na(mor) & mor != 0)

    cat("Computing TF activities...\n")
    
    sample_acts <- decoupleR::decouple( mat     = RNA.counts,
                                        network = universe,
                                        .source = "source",
                                        .target = "target",
                                        minsize = min_targets_size
                                      ) %>%
      dplyr::filter(.data$statistic == "consensus") %>%
      decoupleR::pivot_wider_profile(id_cols     = source,
                                     names_from  = condition,
                                     values_from = score) %>%
      as.matrix() %>%
      t()

  } else {

    if(TF.collection == "CollecTRI"){
      if(is.null(universe)){
        if(file.exists(tf_cache_file)){
          universe = utils::read.csv(tf_cache_file, row.names = 1)
          cat("Using cached TF-target collection from", tf_cache_file, "\n")
        } else {
          universe = decoupleR::get_collectri(organism = 'human', split_complexes = F)
          utils::write.csv(universe, tf_cache_file)
        }
      }
    } else if(TF.collection == "Dorothea"){
      if(is.null(universe)){
        if(file.exists(tf_cache_file)){
          universe = utils::read.csv(tf_cache_file, row.names = 1)
        } else {
          universe = dplyr::filter(dorothea::dorothea_hs, .data$confidence %in% c("A", "B")) %>%
            dplyr::mutate(source = .data$tf) %>%
            dplyr::select(-tf)
          utils::write.csv(universe, tf_cache_file)
        }
      }
    }

    sample_acts <- decoupleR::decouple(mat     = RNA.counts,
                                       network = universe,
                                      .source = "source",
                                      .target = "target", 
                                       minsize = min_targets_size,
                                    ) %>%
                                      dplyr::filter(.data$statistic == "consensus") %>%
                                      decoupleR::pivot_wider_profile(id_cols     = source,
                                                                     names_from  = condition,
                                                                     values_from = score) %>%
                                      as.matrix() %>%
                                      t()

  }

  sample_acts <- sample_acts[colnames(RNA.counts), , drop = FALSE]

  if(return){
    utils::write.csv(sample_acts, paste0("Results/TF_matrix_", file.name, ".csv"))
  }

  result <- data.frame(sample_acts)
  colnames(result) <- make.names(colnames(result))
  return(result)

}

#' Plot top up/down transcription factors
#'
#' Ranks TFs by their activity/statistic value and plots the top and bottom
#' \code{n_top} as a horizontal bar chart colored by direction. Adapted from
#' \code{plot_top_features()} (kinase_tf_mini_tuto/code/utils.R) to work with
#' the single-row (sample x TF) output of \code{compute.TFs.activity()}, which
#' is transposed relative to what that function expects.
#'
#' @param tfs_deg A data frame as returned by \code{compute.TFs.activity()} with a single row (one sample/contrast) and TFs as columns.
#' @param n_top Integer. Number of top upregulated and downregulated TFs to plot.
#'
#' @return A ggplot object.
#'
plot_top_TFs <- function(tfs_deg, n_top = 10) {

  if (nrow(tfs_deg) != 1) {
    stop("tfs_deg must have a single row (one sample/contrast); got ", nrow(tfs_deg), ".")
  }

  data <- as.data.frame(t(tfs_deg))
  colnames(data) <- "value"

  arranged <- data %>%
    tibble::rownames_to_column(var = "id") %>%
    dplyr::arrange(dplyr::desc(value))

  top_up <- dplyr::slice_head(arranged, n = n_top)
  top_down <- dplyr::slice_tail(arranged, n = n_top)

  p <- dplyr::bind_rows(list(up = top_up, down = top_down), .id = "status") %>%
    dplyr::mutate(id = forcats::fct_inorder(id)) %>%
    ggplot2::ggplot(ggplot2::aes(x = value, y = id, fill = status)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
    ggplot2::theme_bw()

  return(p)
}