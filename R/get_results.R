
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

  # Now do the actual work.
  meta_sym <- rlang::syms(spec$metaanalysis_over)

  col_data |>
    group_by(!!! meta_sym) |>
    group_modify(.keep = TRUE, \(sel_dat, key){
      design_act <- if(rlang::is_function(spec$design)) spec$design(sel_dat) else spec$design
      bind_rows(lapply(spec$contrasts, \(cntrst){
        treelabel::test_abundance_changes(sel_dat, design = design_act,
                                                 aggregate_by = {{aggregate_by}},
                                                 contrast = !!cntrst,
                                                 treelabels = all_of(treelabel),
                                                 targets = {{targets}},
                                                 reference = !! rlang::sym(node)) |>
          mutate(contrast = rlang::as_label(cntrst))
      }))
    }) |>
    group_by(contrast, target, treelabel, .add = TRUE) |>
    (\(dat){grouped_df(dat, vars = rev(group_vars(dat)))})() |>
    transmute(obj = tibble(LFC, LFC_se, conf_low = LFC - qnorm(0.975) * LFC_se,
                           conf_high = LFC + qnorm(0.975) * LFC_se, tau = 0)) |>
    recursive_summarize(obj = calc_meta_analysis(obj), .early_stop = 3) |> #  Reduce away all but the last metaanalysis step
    ungroup() |>
    unnest(obj) |>
    dplyr::rename(..meta = dplyr::last(spec$metaanalysis_over)) |>
    mutate(top = node)
}

calc_meta_analysis <- function(obj, mean = "LFC", se = "LFC_se"){
  x <- obj[[mean]]
  y <- obj[[se]]
  x <- x[is.finite(y)]
  y <- y[is.finite(y)]
  if(length(x) == 0){
    tibble({{mean}} := NA, {{se}} := NA, conf_low = NA, conf_high = NA, tau = NA)
  }else if(length(x) == 1){
    tibble({{mean}} := x, {{se}} := y, conf_low = x - qnorm(0.975) * y, conf_high = x + qnorm(0.975) * y, tau = 0)
  }else{
    tryCatch({
      met <- metafor::rma(x, sei = y)
      conf <- confint(met)
      tibble({{mean}} := drop(met$b), {{se}} := met$se, conf_low = met$ci.lb, conf_high = met$ci.ub, tau = sqrt(met$tau2))
    }, error = function(err){
      tibble({{mean}} := NA, {{se}} := NA, conf_low = NA, conf_high = NA, tau = NA)
    })
  }
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
    .by = c(treelabel, top, target, contrast))
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
    res$celltype <- celltype
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
      vars(!!! rlang::syms(spec$metaanalysis_over))
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
    meta_sym <- rlang::syms(spec$metaanalysis_over)
    res <- as_tibble(SummarizedExperiment::colData(psce)) |>
      rowid_to_column(var = "..row_id") |>
      group_by(!!! meta_sym) |>
      group_modify(.keep = TRUE, \(sel_dat, key){
        psce_subset <- psce[,sel_dat$..row_id]
        internal_res <- if(ncol(psce_subset) > 0){
          design_act <- if(rlang::is_function(spec$design)){
            spec$design(sel_dat)
          }else{
            spec$design
          }
          if(spec$de_test == "limma"){
            tryCatch({
              de_limma(SingleCellExperiment::logcounts(psce_subset), design_act, sel_dat, spec$contrasts)
            }, error = function(err){
              tibble()
            })
          }else if(spec$de_test == "glmGamPoi"){
            tryCatch({
              de_glmGamPoi(SingleCellExperiment::counts(psce_subset), design_act, sel_dat, spec$contrasts)
            }, error = function(err){
              tibble()
            })
          }else{
            stop("Invalid de_test argument: '", de_test, "'")
          }
        }
        internal_res
      }) |>
      mutate(lfc_se = pmax(0.1, ifelse(lfc == 0 & t_statistic == 0, Inf, abs(lfc) / abs(t_statistic)))) |>
      transmute(contrast, name, obj = tibble(lfc, lfc_se, conf_low = lfc - qnorm(0.975) * lfc_se,
                             conf_high = lfc + qnorm(0.975) * lfc_se, tau = 0)) |>
      group_by(contrast, name, .add = TRUE) |>
      (\(dat){grouped_df(dat, vars = rev(group_vars(dat)))})()

    # The res_b is very expensive to calculate if there are many rows
    res_a <- res |> filter(n() == 1) |> ungroup(!!! head(meta_sym, n = length(meta_sym) - 1))
    res_b <- res |> filter(n() > 1)  |> recursive_summarize(obj = calc_meta_analysis(obj, mean = "lfc", se = "lfc_se"), .early_stop = 2)

    bind_rows(res_a, res_b) |>
      ungroup() |>
      unnest(obj) |>
      mutate(pval = pnorm(abs(lfc) / lfc_se, lower.tail = FALSE) * 2) |>
      dplyr::rename(..meta = dplyr::last(spec$metaanalysis_over))
  }else{
    design_act <- if(rlang::is_function(spec$design)){
      spec$design(SummarizedExperiment::colData(psce))
    }else{
      spec$design
    }
    res <- if(spec$de_test == "limma"){
      tryCatch({
        de_limma(SingleCellExperiment::logcounts(psce), design_act, SummarizedExperiment::colData(psce), spec$contrasts)
      }, error = function(err){
        tibble()
      })
    }else if(spec$de_test == "glmGamPoi"){
      tryCatch({
        de_glmGamPoi(SingleCellExperiment::counts(psce), design_act, SummarizedExperiment::colData(psce), spec$contrasts)
      }, error = function(err){
        tibble()
      })
    }else{
      stop("Invalid de_test argument: '", de_test, "'")
    }
    res |>
      mutate(lfc_se = pmax(0.1, ifelse(lfc == 0 & t_statistic == 0, Inf, abs(lfc) / abs(t_statistic)))) |>
      transmute(contrast, name, lfc, lfc_se) |>
      mutate(conf_low = lfc - qnorm(0.975) * lfc_se, conf_high = lfc + qnorm(0.975) * lfc_se, tau = 0,
             pval = pnorm(abs(lfc) / lfc_se, lower.tail = FALSE) * 2) |>
      dplyr::rename(..meta = dplyr::last(spec$metaanalysis_over))
  }
}


de_limma <- function(values, design, col_data, quo_contrasts){
  sel <- rowSums(values) > 0
  values_subset <- values[sel,,drop=FALSE]
  extra_res_rows <- data.frame(name = rownames(values)[!sel],
                               pval = rep(1, sum(!sel)), adj_pval = rep(1, sum(!sel)))

  des <- lemur:::handle_design_parameter(design, data = values_subset, col_data = col_data)
  fit <- lemur:::limma_fit(values_subset, des$design_matrix, col_data)
  bind_rows(lapply(seq_along(quo_contrasts), \(idx){
    cntrst <- quo_contrasts[[idx]]
    res <- bind_rows(lemur:::limma_test_de(fit, !!cntrst, des$design_formula),
                     extra_res_rows)
    res$contrast <- names(quo_contrasts)[idx]
    res
  })) |>
    mutate(t_statistic = sign(lfc) * qnorm(pval / 2))
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
    mutate(t_statistic = sign(lfc) * qnorm(pval / 2, lower.tail = FALSE))
}


calculate_meta_differential_expression_results <- function(de_results){
  de_results |> # we summarize over .vals$metaanalysis_over
    arrange(name, contrast) |>
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
    mutate(conf_low = lfc - qnorm(0.975) * lfc_se, conf_high = lfc + qnorm(0.975) * lfc_se) |>
    mutate(adj_pval = p.adjust(pval, method = "BH"), .by = contrast) |>
    transmute(treelabel, celltype, name, pval, adj_pval, lfc, contrast)
}




