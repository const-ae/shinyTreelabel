
get_differential_abundance_results <- function(obj, spec, treelabel, top_sel, da_targets){
  # The first bit is necessary to handle complex pseudobulk_by specifications
  col_data <- obj$col_data()
  aggr_by <- dplyr::transmute(col_data, !!! spec$pseudobulk_by)
  for(new_name in colnames(aggr_by)){
    col_data[[new_name]] <- aggr_by[[new_name]]
  }

  make_differential_abundance_analysis(spec, col_data, treelabel, node = top_sel,
                                       aggregate_by = all_of(colnames(aggr_by)),
                                       targets = vars(!!!rlang::syms(da_targets)))
}

make_differential_abundance_analysis <- function(spec, col_data, treelabel, node, aggregate_by,
                                                 targets = NULL){
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


get_meta_differential_abundance_results <- function(obj, spec, treelabel, top_sel, da_targets){
  calculate_meta_differential_abundance_results(
    obj$get_differential_abundance_results(treelabel, top_sel, da_targets)
  )
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



