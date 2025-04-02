


#'
#'
#' @export
singlecell_treelabel_ui <- function(){
  de_res_ui <- shiny::div(
    conditionalPanel("output.deResultsAreNULL", shiny::p("No cells available that match selection.")),
    conditionalPanel("! output.deResultsAreNULL", shinycssloaders::withSpinner(plotOutput(outputId = "deVolcano"))),
    conditionalPanel("! output.deResultsAreNULL", shinycssloaders::withSpinner(DT::DTOutput(outputId = "deTable"))),
    shinychat::chat_ui("deLLMChat", placeholder = "Feel free to ask follow-up questions")
  )

  metaanalysis_selector <- if(! is.null(.vals$metaanalysis_over)){
    levels <- if(is.factor(.vals$col_data[[.vals$metaanalysis_over]])){
      levels(.vals$col_data[[.vals$metaanalysis_over]])
    }else{
      unique(as.character(.vals$col_data[[.vals$metaanalysis_over]]))
    }
    selectInput("metaanalysisSelector", "Display results for:",
                choices = c(if(length(levels >= 2)) "Meta-Analysis", levels), selected = levels[1])
  }


  navbarPage(title = "Explore your single cell data",
             tabPanel("Overview (UMAP etc.)",
                      selectInput("treelabelSelector0", "Treelabel selector", choices = .vals$treelabel_names,
                                  selected = .vals$treelabel_names[1], multiple = FALSE),
                      treeSelectorUI('celltype_selectorUMAP'),
                      selectInput("redDimSelector", "Dimension reduction", choices = names(.vals$dim_reductions)),
                      plotOutput(outputId = "redDimPlot")
             ),
             tabPanel("Differential Expression",
                      selectInput("treelabelSelector1", "Treelabel selector", choices = .vals$treelabel_names,
                                  selected = .vals$treelabel_names[1], multiple = FALSE),
                      treeSelectorUI('celltype_selectorDEView'),
                      metaanalysis_selector,
                      de_res_ui
             ),
             tabPanel("Differential Abundance",
                      selectInput("treelabelSelector2", "Treelabel selector", choices = .vals$treelabel_names,
                                  selected = .vals$treelabel_names[1], multiple = FALSE),
                      treeSelectorUI('celltype_selectorDiffAbundance'),
                      shinycssloaders::withSpinner(plotOutput(outputId = "diffAbundancePlotsOverview")),
                      shinycssloaders::withSpinner(uiOutput("diffAbundancePlotsFullPlaceholder"))),
             tabPanel("Gene-level analysis",
                      selectInput("treelabelSelector3", "Treelabel selector", choices = .vals$treelabel_names,
                                  selected = .vals$treelabel_names[1], multiple = FALSE),
                      treeSelectorUI('celltype_selectorGeneView'),
                      selectizeInput("geneSelector", "Gene", choices = rownames(.vals$counts)[1:10], selected = rownames(.vals$counts)[1],
                                     options = list(maxOptions = 50, placeholder = "select a gene")),
                      shinycssloaders::withSpinner(plotOutput(outputId = "selGeneExpressionStratified")),
                      shinycssloaders::withSpinner(plotOutput(outputId = "selGeneExpressionContrasted"))
             )
  )
}

#'
#'
#' @export
singlecell_treelabel_server <- function(input, output, session){
  check_init()

  ##### Treelabel selector

  treelabelSelector <- reactiveVal("")
  observeEvent(input$treelabelSelector0, {
    treelabelSelector(input$treelabelSelector0)
    updateSelectInput(session, inputId = "treelabelSelector1", selected = input$treelabelSelector0)
    updateSelectInput(session, inputId = "treelabelSelector2", selected = input$treelabelSelector0)
    updateSelectInput(session, inputId = "treelabelSelector3", selected = input$treelabelSelector0)
  })
  observeEvent(input$treelabelSelector1, {
    treelabelSelector(input$treelabelSelector1)
    updateSelectInput(session, inputId = "treelabelSelector0", selected = input$treelabelSelector1)
    updateSelectInput(session, inputId = "treelabelSelector2", selected = input$treelabelSelector1)
    updateSelectInput(session, inputId = "treelabelSelector3", selected = input$treelabelSelector1)
  })
  observeEvent(input$treelabelSelector2, {
    treelabelSelector(input$treelabelSelector2)
    updateSelectInput(session, inputId = "treelabelSelector0", selected = input$treelabelSelector2)
    updateSelectInput(session, inputId = "treelabelSelector1", selected = input$treelabelSelector2)
    updateSelectInput(session, inputId = "treelabelSelector3", selected = input$treelabelSelector2)
  })
  observeEvent(input$treelabelSelector3, {
    treelabelSelector(input$treelabelSelector3)
    updateSelectInput(session, inputId = "treelabelSelector0", selected = input$treelabelSelector3)
    updateSelectInput(session, inputId = "treelabelSelector1", selected = input$treelabelSelector3)
    updateSelectInput(session, inputId = "treelabelSelector2", selected = input$treelabelSelector3)
  })

  selectable_nodes <- reactiveVal(igraph::V(.vals$tree)$name)
  observeEvent(treelabelSelector(), {
    req(treelabelSelector())
    tl_vec <- .vals$col_data[[treelabelSelector()]]
    all_nodes <- igraph::V(.vals$tree)$name
    score_mat <- tl_score_matrix(tl_vec)
    selectable <- colSums(score_mat > 0, na.rm=TRUE) > 0
    selectable_nodes(colnames(score_mat)[selectable])
  })

  ########## Reduced Dimension View ##########
  cellTypeSelectorUMAPView <- treeSelectorServer("celltype_selectorUMAP", selection_mode = "multiple",
                                                 update_selectable_nodes = selectable_nodes)

  output$redDimPlot <- renderPlot({
    if(input$redDimSelector %in% names(.vals$dim_reductions)){
      sel_nodes <- as.character(cellTypeSelectorUMAPView$selected_nodes())
      .vals$col_data |>
        mutate(reddim = .vals$dim_reductions[[input$redDimSelector]]) |>
        mutate(coloring = tl_tree_filter(!! rlang::sym(.vals$treelabel_names[1]), \(x) sel_nodes) |>
                 tl_name()) |>
        mutate(coloring = factor(ifelse(coloring == .vals$root, NA, coloring), levels = sel_nodes)) |>
        ggplot(aes(x = reddim[,1], y = reddim[,2], color = coloring)) +
          geom_point(size = 0.3) +
          coord_fixed() +
          labs(x = paste0(input$redDimSelector, "1"), y = paste0(input$redDimSelector, "2")) +
          scale_color_discrete(drop = FALSE) +
          guides(color = guide_legend(override.aes = list(size = 2)))
    }else{
      warning("Cannot find '", input$redDimSelector, "' in the stored reduced dimensions ",
              toString(names(.vals$dim_reductions), width = 60))
    }
  })

  cellTypeSelector <- reactiveVal("Immune")
  cellTypeSelectorGeneView <- treeSelectorServer("celltype_selectorGeneView",
                                                 update_selectable_nodes = selectable_nodes,
                                                 update_selected_node = cellTypeSelector)
  observeEvent(cellTypeSelectorGeneView$selected_nodes(), {
    cellTypeSelector(cellTypeSelectorGeneView$selected_nodes())
  })
  cellTypeSelectorDEView <- treeSelectorServer("celltype_selectorDEView",
                                               update_selectable_nodes = selectable_nodes,
                                               update_selected_node = cellTypeSelector)
  observeEvent(cellTypeSelectorDEView$selected_nodes(), {
    cellTypeSelector(cellTypeSelectorDEView$selected_nodes())
  })


  ########## Differential Expression View ##########
  psce <- reactive({
    req(cellTypeSelector())
    req(treelabelSelector())
    print(paste0("Start psce"))
    key <- paste0(treelabelSelector(), "-", cellTypeSelector())
    psce <- if(is.null(.vals$precalc_res$psce[[key]])){
      make_pseudobulk(treelabelSelector(), cellTypeSelector())
    }else{
      print(paste0("Looking up ", key, " in .vals$precalc_res"))
      .vals$precalc_res$psce[[key]]
    }
    print("Done psce")
    psce
  }) |>
    bindCache(cellTypeSelector(), treelabelSelector()) |>
    bindEvent(cellTypeSelector(), treelabelSelector())

  full_de_results <- reactive({
    req(psce())
    print("Start full_de_results")
    key <- rlang::hash(psce())
    res <- if(is.null(.vals$precalc_res$full_de[[key]])){
      make_full_de_results(psce())
    }else{
      print(paste0("Looking up ", key, " in .vals$precalc_res"))
      .vals$precalc_res$full_de[[key]]
    }
    print("Done full_de_results")
    res
  }) |>
    bindCache(psce()) |>
    bindEvent(psce())

  de_result_meta <- reactive({
    print("Start de_result_meta")
    req(full_de_results())
    key <- rlang::hash(full_de_results())
    res <- if(is.null(.vals$precalc_res$meta[[key]])){
      make_meta_analysis(full_de_results())
    }else{
      print(paste0("Looking up ", key, " in .vals$precalc_res"))
      .vals$precalc_res$meta[[key]]
    }
    print("Done de_result_meta")
    res
  }) |>
    bindCache(full_de_results()) |>
    bindEvent(full_de_results())

  de_result_display <- reactive({
    if(is.null(full_de_results()) || nrow(full_de_results()) == 0){
      NULL
    }else if(is.null(.vals$metaanalysis_over)){
      full_de_results()
    }else if(input$metaanalysisSelector != "Meta-Analysis"){
      full_de_results() |>
        filter(!! rlang::sym(.vals$metaanalysis_over) == input$metaanalysisSelector)
    }else{
      de_result_meta()
    }
  })

  output$deTable <- DT::renderDT({
    print("Start renderDT deTable")
    if(is.null(de_result_display())){
      NULL
    }else{
      col_is_numeric <- vapply(de_result_display(), is.numeric, FUN.VALUE = logical(1L))
      slice_min(de_result_display(), pval, n = 50, with_ties = FALSE, by = contrast) |>
        DT::datatable() |>
        DT::formatSignif(which(col_is_numeric))
    }
  })

  output$deResultsAreNULL <- reactive(is.null(de_result_display()))
  outputOptions(output, 'deResultsAreNULL', suspendWhenHidden = FALSE)
  de_gene_message <- reactiveVal()
  # Handle Chat
  chat <- ellmer::chat_ollama(system_prompt = .system_prompt, model = "llama3.3")
  observeEvent(de_result_display(), {
    if(! is.null(de_result_display())){
      sel_genes <- de_result_display() |>
        filter(adj_pval < 0.1) |>
        filter(contrast == contrast[1])
      up_genes <- sel_genes |> filter(lfc > 0) |>  slice_min(pval, n = 5, with_ties = FALSE) |> pull(name)
      down_genes <- sel_genes |> filter(lfc < 0) |>  slice_min(pval, n = 5, with_ties = FALSE) |> pull(name)
      message <- if(length(up_genes) > 0 &  length(down_genes) > 0){
        paste0("I am comparing ", sel_genes$contrast[1], " in ", input$cellTypeSelector, " cells. ",
               "The following genes are significantly up regulated (", toString(up_genes), ") and ",
               "these genes are significantly down regulated (", toString(down_genes), "). ",
               "Do you see a common biological pattern?")
      }else if(length(up_genes) > 0 ){
        paste0("I am comparing ", sel_genes$contrast[1], " in ", input$cellTypeSelector, " cells. ",
               "The following genes are significantly up regulated (", toString(up_genes), "). ",
               "Do you see a common biological pattern?")
      }else if(length(down_genes) > 0 ){
        paste0("I am comparing ", sel_genes$contrast[1], " in ", input$cellTypeSelector, " cells. ",
               "The following genes are significantly down regulated (", toString(down_genes), "). ",
               "Do you see a common biological pattern?")
      }else{
        paste0("I am comparing ", sel_genes$contrast[1], " in ", input$cellTypeSelector,
               " cells, but there are no significantly differentially expressed genes.")
      }
      shinychat::chat_clear(id="deLLMChat")
      de_gene_message(message)
      shinychat::chat_append(id="deLLMChat", message, role = "user")
    }
  })
  observeEvent(input$deLLMChat_user_input, {
    response <- chat$chat_async(input$deLLMChat_user_input)
    shinychat::chat_append("deLLMChat", response)
    NULL
  })
  observeEvent(de_gene_message(), {
    response <- chat$chat_async(de_gene_message())
    shinychat::chat_append("deLLMChat", response)
    NULL # The extra NULL is important because otherwise the function waits for the LLM to complete.
  })


  # Handle Volcano plots
  output$deVolcano <- renderPlot({
    print("Start renderPlot deVolcano")
    de_res <- de_result_display()
    res <- if(is.null(de_res)){
      NULL
    }else{
      de_result_display() |>
        ggplot(aes(x = lfc, y = -log10(pval))) +
        geom_point() +
        geom_hline(data = \(x) x |> group_by(contrast) |> filter(adj_pval < 0.1) |> slice_max(pval, n = 1, with_ties = FALSE),
                   aes(yintercept = -log10(pval))) +
        facet_wrap(vars(contrast)) +
        scale_x_continuous(limits = c(-3.2, 3.2), expand = expansion(add = 0), oob = scales::oob_squish) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    }
    res
  })

  ############## DE of individual Genes ################
  updateSelectizeInput(session, 'geneSelector', choices = rownames(.vals$counts), server = TRUE)


  output$selGeneExpressionStratified <- renderPlot({
    print("Start renderPlot selGeneExpressionStratified")
    sel_gene <- input$geneSelector
    if(! is.null(psce()) && sel_gene %in% rownames(psce())){
      expr_by_sym <- rlang::sym(.vals$gene_expr_by)
      if(! is.null(.vals$metaanalysis_over)){
        metaanalysis_over_sym <- rlang::sym(.vals$metaanalysis_over)
      }
      pl <- SummarizedExperiment::colData(psce()) |>
        as_tibble() |>
        mutate(expr = SingleCellExperiment::logcounts(psce())[sel_gene,]) |>
        ggplot(aes(x = !! expr_by_sym, y = expr)) +
          ggbeeswarm::geom_quasirandom() +
          stat_summary(fun.data = mean_se, geom = "crossbar", color = "red")
      if(! is.null(.vals$metaanalysis_over)){
        pl <- pl + facet_wrap(vars(!! metaanalysis_over_sym))
      }
      print("Done renderPlot selGeneExpressionStratified")
      pl
    }
  })

  output$selGeneExpressionContrasted <- renderPlot({
    print("Start renderPlot selGeneExpressionContrasted")
    sel_gene <- input$geneSelector
    metaanalysis_over_sym <- if(! is.null(.vals$metaanalysis_over)){
      rlang::sym(.vals$metaanalysis_over)
    }else{
      "Data"
    }
    if(! is.null(psce()) && sel_gene %in% rownames(psce())){
      res <- full_de_results() |> filter(name == sel_gene)
      if(! is.null(.vals$metaanalysis_over)){
        meta_res <- de_result_meta() |> filter(name == sel_gene) |> mutate(!!metaanalysis_over_sym := "Meta")
        res <- bind_rows(res, meta_res) |>
          mutate(!!metaanalysis_over_sym := forcats::fct_relevel(as.factor(!!metaanalysis_over_sym), "Meta", after = 0L)) |>
          mutate(!!metaanalysis_over_sym := forcats::fct_rev(!! metaanalysis_over_sym))
      }
      pl <- res |>
        mutate(lfc_se = ifelse(lfc == 0 & t_statistic == 0, Inf, lfc / t_statistic)) |>
        mutate(conf.low = lfc - qnorm(0.975) * lfc_se,
               conf.high = lfc + qnorm(0.975) * lfc_se) |>
        ggplot(aes(x = lfc, y = !!metaanalysis_over_sym)) +
          geom_vline(xintercept  = 0) +
          geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = !!metaanalysis_over_sym == "Meta"), show.legend = FALSE) +
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
          scale_x_continuous(limits = c(-3.2, 3.2), expand = expansion(add = 0), oob = scales::oob_squish) +
          facet_wrap(vars(contrast))
      print("Done renderPlot selGeneExpressionContrasted")
      pl
    }
  })


  ############## Differential Abundance ################
  cellTypeSelectorDAView <- treeSelectorServer("celltype_selectorDiffAbundance", selection_mode = "hierarchical",
                                               update_selectable_nodes = selectable_nodes)

  da_results <- reactive({
    req(cellTypeSelectorDAView$selected_nodes())
    req(cellTypeSelectorDAView$top_selection())
    req(treelabelSelector())

    sel_nodes <- as.character(cellTypeSelectorUMAPView$selected_nodes())
    sel_treelabel <- as.character(treelabelSelector())
    aggr_by <- dplyr::transmute(.vals$col_data, !!! .vals$pseudobulk_by)
    col_data_cp <- .vals$col_data
    for(new_name in colnames(aggr_by)){
      col_data_cp[[new_name]] <- aggr_by[[new_name]]
    }

    top_sel <- cellTypeSelectorDAView$top_selection()
    da_targets <- setdiff(cellTypeSelectorDAView$selected_nodes(), top_sel)
    key <- paste0(sel_treelabel, "-", top_sel)
    res <- if(key %in% names(.vals$precalc_res$full_da)){
      print(paste0("Looking up ", key, " in .vals$precalc_res$full_da"))
      .vals$precalc_res$full_da[[key]] |>
        filter(target %in% da_targets)
    }else{
      make_differential_abundance_analysis(col_data_cp, sel_treelabel, node = top_sel,
                                           aggregate_by = all_of(colnames(aggr_by)), targets = vars(!!!rlang::syms(da_targets)))
    }
    res
  })

  da_meta_res <- reactive({
    req(da_results())
    da_results() |>
      filter(mean(is.finite(LFC_se)) >= 0.5, .by = c(target, ..contrast)) |>
      summarize((\(yi, sei){
        res <- quick_metafor(yi, sei)
        tibble(LFC = res$b, LFC_se = res$se, pval = pmax(1e-22, res$pval), tau = sqrt(res$tau2))
      })(LFC, LFC_se),
      .by = c(target, ..contrast))
  })

  output$diffAbundancePlotsOverview <- renderPlot({
    da_meta_res() |>
      mutate(conf.lower = LFC - qnorm(0.975) * LFC_se,
             conf.high = LFC + qnorm(0.975) * LFC_se) |>
      mutate(includes_zero = conf.lower < 0 & conf.high > 0) |>
      mutate(target = factor(target, levels = cellTypeSelectorDAView$selected_nodes())) |>
      ggplot(aes(x = LFC, y = target)) +
        geom_vline(xintercept = 0) +
        geom_pointrange(aes(xmin = conf.lower, xmax = conf.high, color = includes_zero)) +
        scale_color_manual(values = c("FALSE" = "red", "TRUE" = colorspace::lighten("red", 0.7))) +
        scale_x_continuous(limits = c(-3.2, 3.2), expand = expansion(add = 0), oob = scales::oob_squish) +
        facet_wrap(vars(..contrast))
  })

  output$diffAbundancePlotsFull <- renderPlot({
    req(da_results())
    req(da_meta_res())

    da_results() |>
      bind_rows(da_meta_res() |> mutate(..meta = "Meta")) |>
      mutate(..meta = fct_relevel(..meta, "Meta", after = 0)) |>
      mutate(..meta = fct_rev(..meta)) |>
      mutate(conf.lower = LFC - qnorm(0.975) * LFC_se,
             conf.high = LFC + qnorm(0.975) * LFC_se) |>
      mutate(includes_zero = conf.lower < 0 & conf.high > 0) |>
      mutate(target = fct_rev(factor(target, levels = cellTypeSelectorDAView$selected_nodes()))) |>
      ggplot(aes(x = LFC, y = ..meta)) +
        geom_vline(xintercept = 0) +
        geom_pointrange(aes(xmin = conf.lower, xmax = conf.high, color = includes_zero)) +
        scale_color_manual(values = c("FALSE" = "red", "TRUE" = colorspace::lighten("red", 0.7))) +
        facet_grid(vars(target), vars(..contrast)) +
        scale_x_continuous(limits = c(-3.2, 3.2), expand = expansion(add = 0), oob = scales::oob_squish)
  })

  output$diffAbundancePlotsFullPlaceholder <- renderUI({
    plotOutput(outputId = "diffAbundancePlotsFull",
               height = length(unique(da_meta_res()$target)) * 400)
  })
}


make_pseudobulk <- function(treelabel, celltype, return_selector=FALSE){
  celltype_sym <- rlang::sym(celltype)
  treelabel_sym <- rlang::sym(treelabel)
  cell_sel <- .vals$col_data |>
    mutate(..sel = (tl_eval(!! treelabel_sym, !!celltype_sym) == 1) %|% FALSE) |>
    pull(..sel)
  if(! any(cell_sel)){
    if(return_selector){
      list(psce = NULL, selector = cell_sel)
    }else{
      NULL
    }
  }else{
    extra_pseudobulk_vars <- if(! is.null(.vals$metaanalysis_over)){
      vars(!! rlang::sym(.vals$metaanalysis_over))
    }else{
      vars()
    }
    psce <- pseudobulk(.vals$counts, col_data = .vals$col_data, pseudobulk_by = c(.vals$pseudobulk_by, extra_pseudobulk_vars), filter = cell_sel)
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

make_full_de_results <- function(psce){
  if(is.null(psce)){
    NULL
  }else if(! is.null(.vals$metaanalysis_over)){
    meta_vals <- SummarizedExperiment::colData(psce)[[.vals$metaanalysis_over]]
    levels <- unique(as.character(meta_vals))
    res <- lapply(levels, \(level){
      psce_subset <- psce[,meta_vals == level]
      if(ncol(psce_subset) > 0){
        if(.vals$de_test == "limma"){
          tryCatch({
            de_limma(SingleCellExperiment::logcounts(psce_subset), .vals$design, SummarizedExperiment::colData(psce_subset), .vals$contrasts)
          }, error = function(err){
            NULL
          })
        }else if(.vals$de_test == "glmGamPoi"){
          tryCatch({
            de_glmGamPoi(SingleCellExperiment::counts(psce_subset), .vals$design, SummarizedExperiment::colData(psce_subset), .vals$contrasts)
          }, error = function(err){
            NULL
          })
        }else{
          stop("Invalid de_test argument: '", de_test, "'")
        }
      }
    })
    names(res) <- levels
    bind_rows(res, .id = .vals$metaanalysis_over)
  }else{
    de_limma(SingleCellExperiment::logcounts(psce), .vals$design, SummarizedExperiment::colData(psce), .vals$contrasts)
  }
}

de_limma <- function(values, design, col_data, quo_contrasts){
  des <- lemur:::handle_design_parameter(design, data = values, col_data = col_data)
  fit <- lemur:::limma_fit(values, des$design_matrix, col_data)
  bind_rows(lapply(seq_along(quo_contrasts), \(idx){
    cntrst <- quo_contrasts[[idx]]
    res <- lemur:::limma_test_de(fit, !!cntrst, des$design_formula)
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

make_meta_analysis <- function(full_de_results){
  if(is.null(full_de_results) || ncol(full_de_results) == 0 || nrow(full_de_results) == 0){
    NULL
  }else{
    full_de_results |> # we summarize over .vals$metaanalysis_over
      arrange(name, contrast) |>
      mutate(lfc_se = pmax(0.1, ifelse(lfc == 0 & t_statistic == 0, Inf, abs(lfc) / abs(t_statistic)))) |>
      filter(mean(is.finite(lfc_se)) >= 0.5, .by = c(name, contrast)) |>
      summarize((\(yi, sei){
        res <- quick_metafor(yi, sei)
        tibble(lfc = res$b, lfc_se = res$se, pval = pmax(1e-22, res$pval), tau = sqrt(res$tau2))
      })(lfc ,lfc_se),
      .by = c(name, contrast)) |>
      mutate(adj_pval = p.adjust(pval, method = "BH"), .by = contrast) |>
      transmute(name, pval, adj_pval, t_statistic = lfc/lfc_se, lfc, contrast)
  }
}


make_differential_abundance_analysis <- function(col_data, treelabel, node, aggregate_by,
                                                 targets = NULL){
  tidyr::expand_grid(..meta = unique(col_data[[.vals$metaanalysis_over]]),
                     ..contrast = .vals$contrasts) |>
    rowwise() |>
    group_map(\(dat, .){
      sel_dat <- filter(col_data, !!rlang::sym(.vals$metaanalysis_over) == dat$..meta)
      res <- treelabel::test_abundance_changes(sel_dat, design = .vals$design,
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
