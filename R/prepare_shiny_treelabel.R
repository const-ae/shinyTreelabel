
dummy_precalc_result <- list(
  col_data = \() NULL,
  dim_reductions = \(name) NULL,
  counts = \(gene, cell) NULL,
  da_results = \(treelabel, top_selection, targets) NULL,
  da_meta_results = \(treelabel, top_selection, targets) NULL
)

run_shinyTreelabel <- function(spec, sce = NULL, col_data = NULL, precalc_results = NULL){
  stopifnot(! is.null(sce) || ! is.null(col_data) || ! is.null(precalc_results))

  if(is.null(precalc_results)){
    precalc_results <- dummy_precalc_result
  }else if(is.character(precalc_results)){
    precalc_results <- PrecalculatedResults$new(dbfile = precalc_results)
  }

  obj <- list(
    col_data = \(){
      cd <- col_data %||% precalc_results$col_data() %||% colData(sce)
      as_tibble(cd)
    },
    dim_reductions = \(name){
        precalc_results$dim_reductions(name) %||% reducedDim(sce, name)
    },
    counts = \(gene, cell){
      precalc_results$dim_reductions(name) %||% assay(sce, spec$count_assay_name)[gene,cell]
    },
    da_results = \(treelabel, top_selection, targets){
      precalc_results$da_results(treelabel, top_selection, targets) %||%
        get_differential_abundance_results(obj, spec, treelabel, top_selection, targets)
    },
    da_meta_results = \(treelabel, top_selection, targets){
      precalc_results$da_meta_results(treelabel, top_selection, targets) %||%
        get_meta_differential_abundance_results(obj, spec, treelabel, top_selection, targets)
    }
  )

  shiny::shinyApp(singlecell_treelabel_ui(spec), singlecell_treelabel_server2_gen(spec, obj))
}


#' Prepare App
#'
#' @export
init_shinyTreelabel <- function(sce, treelabels = where(is_treelabel),
                                design = ~ 1,
                                pseudobulk_by = NULL,
                                metaanalysis_over = NULL,
                                gene_expr_by = NULL,
                                contrasts = vars(Intercept = cond()),
                                de_test = c("limma", "glmGamPoi"),
                                dim_reductions = NULL,
                                count_assay = "counts",
                                treelabel_threshold = 0.5){

  if(! is(sce, "SummarizedExperiment")){
    stop("Cannot handle 'sce' of class ", toString(class(sce), width = 60))
  }

  col_data <- as_tibble(SummarizedExperiment::colData(sce))
  de_test <- rlang::arg_match(de_test)
  if(is.null(dim_reductions)){
    dim_reductions <- SingleCellExperiment::reducedDimNames(sce)
  }
  contrasts <- rlang::quos_auto_name(contrasts)
  design_variables <- if(rlang::is_function(design)){
    lemur:::design_variable_to_quosures(design(col_data), data = col_data)
  }else{
    lemur:::design_variable_to_quosures(design, data = col_data)
  }

  col_data <- col_data |>
    mutate(across({{treelabels}}, \(x) tl_modify(x, .scores >= treelabel_threshold)))
  vec <- vctrs::vec_ptype_common(!!! dplyr::select(col_data, where(is_treelabel)))
  treelabel_names <- names(tidyselect::eval_select({{treelabels}}, data = col_data))

  list(
    treelabel_names = treelabel_names,
    tree = tl_tree(vec),
    root = tl_tree_root(vec),
    metaanalysis_over = metaanalysis_over,
    gene_expr_by = gene_expr_by,
    pseudobulk_by = c(pseudobulk_by, design_variables),
    design = design,
    contrasts = contrasts,
    de_test = de_test,
    count_assay_name = count_assay,
    dim_reduction_names = dim_reductions %||% character(0L),
    gene_names = rownames(sce) %||% as.character(seq_len(1:nrow(sce)))
  )
}





check_init <- function(){
  stopifnot(all(c("col_data","treelabel_names", "tree", "root",  "dim_reductions",  "counts", "metaanalysis_over",
                  "gene_expr_by", "pseudobulk_by", "design", "contrasts", "de_test", "precalc_res") %in% names(.vals)))
}





precalculate_results2 <- function(spec, sce, output = c("dim_reductions", "da", "de", "de_meta", "pseudobulks"),
                                  treelabel_sel = NULL, dbfile = ":memory:", verbose = TRUE){
  storage <- PrecalculatedResults$new(dbfile)

  col_data <- as_tibble(colData(sce))
  storage$add_col_data(col_data)

  if(is.null(treelabel_sel)){
    treelabel_sel <- spec$treelabel_names
  }else{
    stopifnot(all(treelabel_sel %in% spec$treelabel_names))
  }
  nodes <- igraph::V(spec$tree)$name
  aggr_by <- dplyr::transmute(col_data, !!! spec$pseudobulk_by)
  col_data_cp <- col_data
  for(new_name in colnames(aggr_by)){
    col_data_cp[[new_name]] <- aggr_by[[new_name]]
  }
  if("dim_reductions" %in% output){
    all_redDims <- reducedDims(sce)[spec$dim_reduction_names] |>
      lapply(\(x) if(ncol(x) == 0){
        matrix(0, nrow = length(spec$gene_names), ncol = 2)
      }else if(ncol(x) == 1){
        cbind(x, 0)
      }else{
        x[,1:2,drop=FALSE]
      }) |>
      enframe(name = "name", value = "embedding")
    storage$add_reducedDimensions_rows(all_redDims)
  }

  if("da" %in% output){
    for(ts in treelabel_sel){
      i <- 1; total_steps <- length(nodes)
      n <- spec$root; step <- "init"
      if(verbose){
        cli::cli_progress_step("Step {i}/{total_steps} | ETA: {cli::pb_eta}: '{ts}' calculating '{n}' ({step})",
                               msg_done = paste0(" Finished '", ts, "'"), total = total_steps)
      }
      for(n in nodes){
        full_da <- make_differential_abundance_analysis(spec, col_data_cp, treelabel = ts, node = n,
                                                        aggregate_by = all_of(colnames(aggr_by)))
        if(nrow(full_da) > 0){
          full_da$top <- n
        }
        storage$add_da_rows(full_da)
        if(verbose) tryCatch(cli::cli_progress_update(force = TRUE, inc = 1), error = \(err){})
      }
      i <- i+1
    }
    if(verbose) cli::cli_process_done()
    meta <- tbl(storage$con, "da") |>
      dplyr::collect() |>
      calculate_meta_differential_abundance_results()
    storage$add_da_meta_rows(meta)
  }

  storage
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

