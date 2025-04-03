.vals <- new.env(parent = emptyenv())

#' Prepare App
#'
#' @export
init_shinyTreelabel <- function(data, treelabels = where(is_treelabel),
                                design = ~ 1,
                                pseudobulk_by = NULL,
                                metaanalysis_over = NULL,
                                gene_expr_by = NULL,
                                contrasts = vars(Intercept = cond()),
                                de_test = c("limma", "glmGamPoi"),
                                col_data = NULL,
                                count_assay = "counts",
                                seurat_assay = "RNA",
                                treelabel_threshold = 0.5,
                                reduced_dims = NULL){

  de_test <- rlang::arg_match(de_test)
  if(is.matrix(data) || is(data, "Matrix")){
    counts <- data
    if(is.null(col_data)){
      stop("If 'data' is a matrix, the 'col_data' must not be 'NULL'")
    }
    if("PCA" %in% reduced_dims){
      logcounts <- transformGamPoi::shifted_log_transform(counts)
      dim_reductions <- list(PCA = irlba::prcomp_irlba(t(logcounts), n = 2)$x)
    }else{
      dim_reductions <- list()
    }
  }else if(is(data, "SummarizedExperiment")){
    counts <- SummarizedExperiment::assay(data, count_assay)
    col_data <- cbind(SummarizedExperiment::colData(data), col_data)
    dim_reductions <- lapply(SingleCellExperiment::reducedDims(data), \(x) x[,1:2,drop=FALSE])
    if(is.null(reduced_dims)){
      reduced_dims <- SingleCellExperiment::reducedDimNames(data)
    }
  }else if(inherits(data, "Seurat")){
    counts <- data[[seurat_assay]]$counts
    col_data <- bind_cols(data[[]], col_data)
    red <- SeuratObject::Reductions(data)
    dim_reductions <- lapply(red, \(name) SeuratObject::Embeddings(data, reduction = name)[,1:2,drop=FALSE])
    names(dim_reductions) <- red
    if(is.null(reduced_dims)){
      reduced_dims <- red
    }
  }else{
    stop("Cannot handle 'data' of class ", toString(class(data), width = 60))
  }
  contrasts <- rlang::quos_auto_name(contrasts)
  design_variables <- if(rlang::is_function(design)){
    lemur:::design_variable_to_quosures(design(col_data), data = col_data)
  }else{
    lemur:::design_variable_to_quosures(design, data = col_data)
  }
  col_data <- as_tibble(col_data) |>
    mutate(across({{treelabels}}, \(x) tl_modify(x, .scores >= treelabel_threshold)))


  vec <- vctrs::vec_ptype_common(!!! dplyr::select(col_data, where(is_treelabel)))


  .vals$col_data <- col_data
  .vals$treelabel_names <- names(tidyselect::eval_select({{treelabels}}, data = col_data))
  .vals$tree <- tl_tree(vec)
  .vals$root <- tl_tree_root(vec)
  .vals$dim_reductions <- dim_reductions[intersect(names(dim_reductions), reduced_dims)]
  .vals$counts <- counts
  .vals$metaanalysis_over <- metaanalysis_over
  .vals$gene_expr_by <- gene_expr_by
  .vals$pseudobulk_by <- c(pseudobulk_by, design_variables)
  .vals$design <- design
  .vals$contrasts <- contrasts
  .vals$de_test <- de_test
  .vals$precalc_res <- list(psce = list(), full_de = list(), meta = list())
}


check_init <- function(){
  stopifnot(all(c("col_data","treelabel_names", "tree", "root",  "dim_reductions",  "counts", "metaanalysis_over",
                  "gene_expr_by", "pseudobulk_by", "design", "contrasts", "de_test", "precalc_res") %in% names(.vals)))
}



#' Speed-up app by precalculating the results
#'
#' @param spec the specification which results to precalculate. Either a
#'   string 'all' or 'nothing' or a list with two elements 'treelabels' and
#'   'nodes' which list the elements to precalculate.
#' @param verbose should a progress bar be printed.
#' @param results a list produced by `precalculate_results`.
#'
#' @returns a list with all the precalculated results that needs to be
#'   passed to `install_precalculated_results`.
#'
#' @export
precalculate_results <- function(spec = c("all", "nothing"), verbose = TRUE){
  check_init()
  valid_precalculate_spec <- (is.character(spec) && precalculateprecalculates[1] %in% c("all", "nothing")) ||
    (is.list(spec) && all(names(spec) %in% c("treelabels", "nodes")))
  stopifnot(valid_precalculate_spec)

  treelabelSelectors <- .vals$treelabel_names
  nodes <- igraph::V(.vals$tree)$name
  aggr_by <- dplyr::transmute(.vals$col_data, !!! .vals$pseudobulk_by)
  col_data_cp <- .vals$col_data
  for(new_name in colnames(aggr_by)){
    col_data_cp[[new_name]] <- aggr_by[[new_name]]
  }


  res <- list(psce = list(), full_de = list(), meta = list(),
              full_da = list())
  psce_key_lookup <- rlang::new_environment()
  # Calculate psce
  for(ts in treelabelSelectors){
    if(check_treelabel_precalculuate_spec(ts, spec)){
      i <- 1
      total_steps <- length(nodes)
      n <- .vals$root
      step <- "init"
      if(verbose){
        cli::cli_progress_step("Step {i}/{total_steps} | ETA: {cli::pb_eta}: '{ts}' calculating '{n}' ({step})",
                               msg_done = paste0(" Finished '", ts, "'"), total = total_steps)
      }
      for(n in nodes){
        if(check_nodes_precalculuate_spec(n, spec)){
          step <- "pseudobulk"
          if(verbose) tryCatch(cli::cli_progress_update(force = TRUE, set = i), error = \(err){})
          psce_val <-  make_pseudobulk(ts, n, return_selector = TRUE)
          key <- rlang::hash(psce_val$selector)
          step <- "differential abundance"
          if(verbose) tryCatch(cli::cli_progress_update(force = TRUE, inc = 0.25), error = \(err){})
          full_da <- make_differential_abundance_analysis(col_data_cp, treelabel = ts, node = n,
                                                          aggregate_by = all_of(colnames(aggr_by)))
          if(exists(key, envir = psce_key_lookup)){
            psce <- res$psce[[psce_key_lookup[[key]]]]
            full_de_res <- res$full_de[[rlang::hash(psce)]]
            meta_res <- res$meta[[rlang::hash(full_de_res)]]
          }else{
            psce_key_lookup[[key]] <- paste0(ts, "-", n)
            psce <- psce_val$psce
            step <- "differential expression"
            if(verbose) tryCatch(cli::cli_progress_update(force = TRUE, inc = 0.25), error = \(err){})
            full_de_res <- suppressWarnings({
              make_full_de_results(psce)
            })
            step <- "meta analysis"
            if(verbose) tryCatch(cli::cli_progress_update(force = TRUE, inc = 0.25), error = \(err){})
            meta_res <- suppressWarnings({
              make_meta_analysis(full_de_res)
            })
          }
          res$full_da[[paste0(ts, "-", n)]] <- full_da
          res$psce[[paste0(ts, "-", n)]] <- psce
          res$full_de[[rlang::hash(psce)]] <- full_de_res
          res$meta[[rlang::hash(full_de_res)]] <- meta_res
        }
        i <- i+1
      }
      if(verbose) {
        cli::cli_process_done()
      }
    }
  }
  res
}


#' @export
#' @rdname precalculate_results
install_precalculated_results <- function(results){
  stopifnot(is.list(results))
  stopifnot(all(c("psce", "full_de", "meta") %in% names(results)))
  .vals$precalc_res <- results
}

check_treelabel_precalculuate_spec <- function(treelabel_name, spec){
 if(is.character(spec) && spec[1] == "all"){
   TRUE
 }else if(is.character(spec) && spec[1] == "nothing"){
   FALSE
 }else{
   spec$treelabels[1] == "all" || treelabel_name %in% spec$treelabels
 }
}

check_nodes_precalculuate_spec <- function(node, spec){
  if(is.character(spec) && spec == "all"){
    TRUE
  }else if(is.character(spec) && spec == "nothing"){
    FALSE
  }else{
    spec$nodes[1] == "all" || node %in% spec$nodes
  }
}

