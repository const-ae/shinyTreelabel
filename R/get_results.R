
########## Differential Abundance Section ##########

get_differential_abundance_results <- function(obj, spec, treelabel, top_sel, da_targets){
  # The first bit is necessary to handle complex pseudobulk_by specifications
  col_data <- obj$col_data()
  aggr_by <- dplyr::transmute(col_data, !!! spec$pseudobulk_by)
  for(new_name in colnames(aggr_by)){
    col_data[[new_name]] <- aggr_by[[new_name]]
  }
  calculate_differential_abundance_results(spec, col_data, treelabel, node = top_sel,
                                       aggregate_by = all_of(colnames(aggr_by)),
                                       targets = vars(!!!rlang::syms(da_targets)))
}

calculate_differential_abundance_results <- function(spec, col_data, treelabel, node, aggregate_by,
                                                 targets = NULL){
  # Check if it makes sense to continue
  celltype_sym <- rlang::sym(node)
  treelabel_sym <- rlang::sym(treelabel)
  cell_sel <- col_data |>
    mutate(..sel = (tl_eval(!! treelabel_sym, !!celltype_sym) == 1) %|% FALSE) |>
    pull(..sel)
  if(! any(cell_sel)){
    return(NULL)
  }

  # Now do the actual work
  tidyr::expand_grid(..meta = unique(col_data[[spec$metaanalysis_over]]),
                     ..contrast = spec$contrasts) |>
    rowwise() |>
    group_map(\(dat, .){
      sel_dat <- filter(col_data, !!rlang::sym(spec$metaanalysis_over) == dat$..meta)
      design_act <- if(rlang::is_function(spec$design)){
        spec$design(sel_dat)
      }else{
        spec$design
      }
      res <- treelabel::test_abundance_changes(sel_dat, design = design_act,
                                               aggregate_by = {{aggregate_by}},
                                               contrast = !!dat$..contrast[[1]],
                                               treelabels = all_of(treelabel),
                                               targets = {{targets}},
                                               reference = !! rlang::sym(node))
      res |>
        mutate(..meta = dat$..meta, ..contrast = rlang::as_label(dat$..contrast[[1]]))
    }) |>
    dplyr::bind_rows()
}


calculate_meta_differential_abundance_results <- function(da_results){
  da_results |>
    summarize((\(yi, sei){
      if(mean(is.finite(sei)) >= 0.5){
        res <- quick_metafor(yi, sei)
        tibble(LFC = res$b, LFC_se = res$se, pval = pmax(1e-22, res$pval),
               tau = sqrt(res$tau2))
      }else{
        tibble(LFC = 0, LFC_se = Inf, pval = 1, tau = 0)
      }
    })(LFC, LFC_se),
    .by = c(treelabel, top, target, ..contrast))
}



########## Differential Expression Section ##########

get_differential_expression_results <- function(obj, spec, treelabel, celltype, gene){
  # The first bit is necessary to handle complex pseudobulk_by specifications
  col_data <- obj$col_data()
  aggr_by <- dplyr::transmute(col_data, !!! spec$pseudobulk_by)
  for(new_name in colnames(aggr_by)){
    col_data[[new_name]] <- aggr_by[[new_name]]
  }

  calculate_differential_expression_results(spec, col_data, obj$counts(), treelabel, celltype,
                                       aggregate_by = vars(!!! rlang::syms(colnames(aggr_by))))
}


calculate_differential_expression_results <- function(spec, col_data, counts, treelabel, celltype, aggregate_by){
  # Make pseudobulk
  psce <- make_pseudobulk(spec, col_data, counts, treelabel, celltype, aggregate_by)
  # Calculate DE
  res <- make_full_de_results(spec, psce)
  if(! is.null(res)){
    res$treelabel <- treelabel
  }
  res
}



make_pseudobulk <- function(spec, col_data, counts, treelabel, celltype, aggregate_by, return_selector=FALSE){
  celltype_sym <- rlang::sym(celltype)
  treelabel_sym <- rlang::sym(treelabel)
  cell_sel <- col_data |>
    mutate(..sel = (tl_eval(!! treelabel_sym, !!celltype_sym) == 1) %|% FALSE) |>
    pull(..sel)
  if(! any(cell_sel)){
    if(return_selector){
      list(psce = NULL, selector = cell_sel)
    }else{
      NULL
    }
  }else{
    extra_pseudobulk_vars <- if(! is.null(spec$metaanalysis_over)){
      vars(!! rlang::sym(spec$metaanalysis_over))
    }else{
      vars()
    }
    psce <- pseudobulk(counts, col_data = col_data, pseudobulk_by = c(aggregate_by, extra_pseudobulk_vars), filter = cell_sel)
    if(return_selector){
      list(psce = psce, selector = cell_sel)
    }else{
      psce
    }
  }
}

pseudobulk <- function(counts, col_data, pseudobulk_by, filter = TRUE){
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=counts), colData = col_data)
  sce <- sce[,filter]
  res <- glmGamPoi::pseudobulk(sce, group_by = pseudobulk_by, verbose = FALSE)
  SingleCellExperiment::logcounts(res) <- transformGamPoi::shifted_log_transform(res)
  res
}

make_full_de_results <- function(spec, psce){
  if(is.null(psce)){
    NULL
  }else if(! is.null(spec$metaanalysis_over)){
    meta_vals <- SummarizedExperiment::colData(psce)[[spec$metaanalysis_over]]
    levels <- unique(as.character(meta_vals))
    res <- lapply(levels, \(level){
      psce_subset <- psce[,meta_vals == level]
      if(ncol(psce_subset) > 0){
        design_act <- if(rlang::is_function(spec$design)){
          spec$design(SummarizedExperiment::colData(psce_subset))
        }else{
          spec$design
        }
        if(spec$de_test == "limma"){
          tryCatch({
            de_limma(SingleCellExperiment::logcounts(psce_subset), design_act, SummarizedExperiment::colData(psce_subset), spec$contrasts)
          }, error = function(err){
            NULL
          })
        }else if(spec$de_test == "glmGamPoi"){
          tryCatch({
            de_glmGamPoi(SingleCellExperiment::counts(psce_subset), design_act, SummarizedExperiment::colData(psce_subset), spec$contrasts)
          }, error = function(err){
            NULL
          })
        }else{
          stop("Invalid de_test argument: '", de_test, "'")
        }
      }
    })
    names(res) <- levels
    bind_rows(res, .id = spec$metaanalysis_over)
  }else{
    design_act <- if(rlang::is_function(spec$design)){
      spec$design(SummarizedExperiment::colData(psce))
    }else{
      spec$design
    }
    de_limma(SingleCellExperiment::logcounts(psce), design_act, SummarizedExperiment::colData(psce), spec$contrasts)
  }
}


de_limma <- function(values, design, col_data, quo_contrasts){
  sel <- rowSums(values) > 0
  values_subset <- values[sel,,drop=FALSE]
  extra_res_rows <- data.frame(name = rownames(values)[!sel],
                               pval = 1, adj_pval = 1)

  des <- lemur:::handle_design_parameter(design, data = values_subset, col_data = col_data)
  fit <- lemur:::limma_fit(values_subset, des$design_matrix, col_data)
  bind_rows(lapply(seq_along(quo_contrasts), \(idx){
    cntrst <- quo_contrasts[[idx]]
    res <- bind_rows(lemur:::limma_test_de(fit, !!cntrst, des$design_formula),
                     extra_res_rows)
    res$contrast <- names(quo_contrasts)[idx]
    res
  }))
}

de_glmGamPoi <- function(values, design, col_data, quo_contrasts){
  fit <- glmGamPoi::glm_gp(values, design, col_data = col_data, size_factors = "ratio",
                           ridge_penalty = 1e-5)
  bind_rows(lapply(seq_along(quo_contrasts), \(idx){
    cntrst <- quo_contrasts[[idx]]
    res <- glmGamPoi::test_de(fit, !!cntrst)
    res$contrast <- names(quo_contrasts)[idx]
    res
  })) |>
    mutate(t_statistic = sign(lfc) * sqrt(abs(f_statistic)))
}


calculate_meta_differential_expression_results <- function(de_results){
  de_results |> # we summarize over .vals$metaanalysis_over
    arrange(name, contrast) |>
    mutate(lfc_se = pmax(0.1, ifelse(lfc == 0 & t_statistic == 0, Inf, abs(lfc) / abs(t_statistic)))) |>
    summarize((\(y, se){
      if(sum(is.finite(se) & se <= 10) >= 3){
        y <- y[is.finite(se)]
        se <- se[is.finite(se)]
        res <- quick_metafor(y, se)
        tibble(lfc = res$b, lfc_se = res$se, pval = pmax(1e-308, res$pval), tau = sqrt(res$tau2))
      }else {
        tibble(lfc = NA, lfc_se = NA, pval = NA, tau = NA)
      }
    })(lfc ,lfc_se),
    .by = c(name, contrast, treelabel, celltype)) |>
    mutate(adj_pval = p.adjust(pval, method = "BH"), .by = contrast) |>
    transmute(treelabel, celltype, name, pval, adj_pval, t_statistic = lfc/lfc_se, lfc, contrast)
}




