


#'
#'
#' @export
singlecell_treelabel_ui <- function(spec){
  # de_res_ui <- shiny::div(
  #   conditionalPanel("output.deResultsAreNULL", shiny::p("No cells available that match selection.")),
  #   conditionalPanel("! output.deResultsAreNULL", shinycssloaders::withSpinner(plotOutput(outputId = "deVolcano"))),
  #   conditionalPanel("! output.deResultsAreNULL", shinycssloaders::withSpinner(DT::DTOutput(outputId = "deTable"))),
  #   shinychat::chat_ui("deLLMChat", placeholder = "Feel free to ask follow-up questions")
  # )

  metaanalysis_selector <- if(! is.null(spec$metaanalysis_over)){
    levels <- spec$metaanalysis_levels
    selectInput("metaanalysisSelector", "Display results for:",
                choices = c(if(length(levels >= 2)) "Meta-Analysis", levels), selected = levels[1])
  }


  navbarPage(title = "Explore your single cell data",
             tabPanel("Overview (UMAP etc.)",
                      selectInput("treelabelSelector0", "Treelabel selector", choices = spec$treelabel_names,
                                  selected = spec$treelabel_names[1], multiple = FALSE),
                      treeSelectorUI('celltype_selectorUMAP'),
                      selectInput("redDimSelector", "Dimension reduction", choices = spec$dim_reduction_names),
                      plotOutput(outputId = "redDimPlot")
             ),
             tabPanel("Differential Expression",
                      selectInput("treelabelSelector1", "Treelabel selector", choices = spec$treelabel_names,
                                  selected = spec$treelabel_names[1], multiple = FALSE),
                      treeSelectorUI('celltype_selectorDEView'),
                      metaanalysis_selector,
                      shinycssloaders::withSpinner(plotOutput(outputId = "deVolcano")),
                      shinycssloaders::withSpinner(DT::DTOutput(outputId = "deTable"))
             ),
             tabPanel("Differential Abundance",
                      selectInput("treelabelSelector2", "Treelabel selector", choices = spec$treelabel_names,
                                  selected = spec$treelabel_names[1], multiple = FALSE),
                      treeSelectorUI('celltype_selectorDiffAbundance'),
                      shinycssloaders::withSpinner(plotOutput(outputId = "diffAbundancePlotsOverview")),
                      shinycssloaders::withSpinner(uiOutput("diffAbundancePlotsFullPlaceholder"))
                      ),
             # tabPanel("Gene-level analysis",
             #          selectInput("treelabelSelector3", "Treelabel selector", choices = obj$treelabel_names,
             #                      selected = obj$treelabel_names[1], multiple = FALSE),
             #          treeSelectorUI('celltype_selectorGeneView'),
             #          selectizeInput("geneSelector", "Gene", choices = rownames(obj$counts)[1:10], selected = rownames(obj$counts)[1],
             #                         options = list(maxOptions = 50, placeholder = "select a gene")),
             #          shinycssloaders::withSpinner(plotOutput(outputId = "selGeneExpressionStratified")),
             #          shinycssloaders::withSpinner(plotOutput(outputId = "selGeneExpressionContrasted"))
             # )
  )
}


singlecell_treelabel_server2_gen <- function(spec, obj) { function(input, output, session){
  ##### Treelabel selector
  # A lot of boiler plate, but this just keeps the treelabelSelectors across the tabs in sync
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

  selectable_nodes <- reactiveVal(igraph::V(spec$tree)$name)
  observeEvent(treelabelSelector(), {
    req(treelabelSelector())
    tl_vec <- obj$col_data()[[treelabelSelector()]]
    all_nodes <- igraph::V(spec$tree)$name
    score_mat <- tl_score_matrix(tl_vec)
    selectable <- colSums(score_mat > 0, na.rm=TRUE) > 0
    selectable_nodes(colnames(score_mat)[selectable])
  })

  ########## Reduced Dimension View ##########
  cellTypeSelectorUMAPView <- treeSelectorServer("celltype_selectorUMAP", spec, selection_mode = "multiple",
                                                 update_selectable_nodes = selectable_nodes)
  color_choice_fnc_gen <- color_generator_fnc("multiple")

  reduced_dim <- reactive({
    obj$dim_reductions(input$redDimSelector)
  })

  output$redDimPlot <- renderPlot({
    req(reduced_dim())

    sel_nodes <- as.character(cellTypeSelectorUMAPView$selected_nodes())
    color_map <- structure(color_choice_fnc_gen(sel_nodes)(sel_nodes), names = sel_nodes)
    obj$col_data() |>
      mutate(reddim = reduced_dim()) |>
      mutate(coloring = tl_name(tl_tree_filter(!! rlang::sym(spec$treelabel_names[1]), \(x) setdiff(sel_nodes, spec$root)))) |>
      mutate(coloring = factor(ifelse(coloring == spec$root, NA, coloring), levels = sel_nodes)) |>
      ggplot(aes(x = .data$reddim[,1], y = .data$reddim[,2], color = .data$coloring)) +
        geom_point(size = 0.3) +
        coord_fixed() +
        labs(x = paste0(input$redDimSelector, "1"), y = paste0(input$redDimSelector, "2")) +
        scale_color_manual(values = color_map, drop = FALSE) +
        guides(color = guide_legend(override.aes = list(size = 2)))

  })

  ####### Differential abundance output #######
  cellTypeSelectorDAView <- treeSelectorServer("celltype_selectorDiffAbundance", spec,
                                               selection_mode = "hierarchical",
                                               update_selectable_nodes = selectable_nodes)

  da_results <- reactive({
    req(cellTypeSelectorDAView$selected_nodes(), cellTypeSelectorDAView$top_selection(),
        treelabelSelector())

    obj$da_results(
      treelabel = as.character(treelabelSelector()),
      top_selection = cellTypeSelectorDAView$top_selection(),
      targets = setdiff(cellTypeSelectorDAView$selected_nodes(), cellTypeSelectorDAView$top_selection())
    )
  })

  da_meta_res <- reactive({
    req(da_results(), cellTypeSelectorDAView$selected_nodes(),
        cellTypeSelectorDAView$top_selection(), treelabelSelector())

    obj$da_meta_results(
      treelabel = as.character(treelabelSelector()),
      top_selection = cellTypeSelectorDAView$top_selection(),
      targets = setdiff(cellTypeSelectorDAView$selected_nodes(), cellTypeSelectorDAView$top_selection())
    )

  })


  output$diffAbundancePlotsOverview <- renderPlot({
    req(da_meta_res())

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
      mutate(..meta = forcats::fct_relevel(..meta, "Meta", after = 0)) |>
      mutate(..meta = forcats::fct_rev(..meta)) |>
      mutate(conf.lower = LFC - qnorm(0.975) * LFC_se,
             conf.high = LFC + qnorm(0.975) * LFC_se) |>
      mutate(includes_zero = conf.lower < 0 & conf.high > 0) |>
      mutate(target = forcats::fct_rev(factor(target, levels = cellTypeSelectorDAView$selected_nodes()))) |>
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

  ####### Differential expression output #######
  cellTypeSelector <- reactiveVal(spec$root)
  cellTypeSelectorGeneView <- treeSelectorServer("celltype_selectorGeneView", spec,
                                update_selectable_nodes = selectable_nodes, update_selected_node = cellTypeSelector)
  observeEvent(cellTypeSelectorGeneView$selected_nodes(), cellTypeSelector(cellTypeSelectorGeneView$selected_nodes()))
  cellTypeSelectorDEView <- treeSelectorServer("celltype_selectorDEView", spec,
                                update_selectable_nodes = selectable_nodes, update_selected_node = cellTypeSelector)
  observeEvent(cellTypeSelectorDEView$selected_nodes(), cellTypeSelector(cellTypeSelectorDEView$selected_nodes()))

  de_results <- reactive({
    req(cellTypeSelector, treelabelSelector())
    obj$de_results(
      treelabel = as.character(treelabelSelector()),
      celltype = cellTypeSelector()
    )
  })

  de_meta_results <- reactive({
    req(cellTypeSelector, treelabelSelector())
    obj$de_meta_results(
      treelabel = as.character(treelabelSelector()),
      celltype = cellTypeSelector()
    )
  })

  de_result_display <- reactive({
    req(de_results())
    if(is.null(spec$metaanalysis_over)){
      de_results()
    }else if(input$metaanalysisSelector != "Meta-Analysis"){
      de_results() |>
        filter(!! rlang::sym(spec$metaanalysis_over) == input$metaanalysisSelector)
    }else{
      de_meta_results()
    }
  })

  output$deTable <- DT::renderDT({
    req(de_result_display())
    col_is_numeric <- vapply(de_result_display(), is.numeric, FUN.VALUE = logical(1L))
    slice_min(de_result_display(), pval, n = 50, with_ties = FALSE, by = contrast) |>
      DT::datatable() |>
      DT::formatSignif(which(col_is_numeric))
  })
  output$deVolcano <- renderPlot({
    req(de_result_display())
    de_result_display() |>
      ggplot(aes(x = lfc, y = -log10(pval))) +
      geom_point(aes(color = adj_pval < 0.1)) +
      geom_hline(data = \(x) x |> group_by(contrast) |> filter(adj_pval < 0.1) |> slice_max(pval, n = 1, with_ties = FALSE),
                 aes(yintercept = -log10(pval))) +
      ggrepel::geom_text_repel(data = \(x) slice_min(x, pval, n = 10, with_ties = FALSE, by = contrast),
                               aes(label = name), size = 4) +
      facet_wrap(vars(contrast)) +
      scale_x_continuous(limits = c(-3.2, 3.2), expand = expansion(add = 0), oob = scales::oob_squish) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey"))
  })

}}

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
  cellTypeSelectorUMAPView <- treeSelectorServer("celltype_selectorUMAP", obj, selection_mode = "multiple",
                                                 update_selectable_nodes = selectable_nodes)

  output$redDimPlot <- renderPlot({
    if(input$redDimSelector %in% names(.vals$dim_reductions)){
      sel_nodes <- as.character(cellTypeSelectorUMAPView$selected_nodes())
      .vals$col_data |>
        mutate(reddim = .vals$dim_reductions[[input$redDimSelector]]) |>
        mutate(coloring = tl_tree_filter(!! rlang::sym(.vals$treelabel_names[1]), \(x) sel_nodes) |>
                 tl_name()) |>
        mutate(coloring = factor(ifelse(coloring == .vals$root, NA, coloring), levels = sel_nodes)) |>
        ggplot(aes(x = .data$reddim[,1], y = .data$reddim[,2], color = .data$coloring)) +
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

  cellTypeSelector <- reactiveVal(obj$root)
  cellTypeSelectorGeneView <- treeSelectorServer("celltype_selectorGeneView", obj,
                                                 update_selectable_nodes = selectable_nodes,
                                                 update_selected_node = cellTypeSelector)
  observeEvent(cellTypeSelectorGeneView$selected_nodes(), {
    cellTypeSelector(cellTypeSelectorGeneView$selected_nodes())
  })
  cellTypeSelectorDEView <- treeSelectorServer("celltype_selectorDEView", obj,
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
    res <- if(is.null(de_res) || nrow(de_res) == 0){
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
  cellTypeSelectorDAView <- treeSelectorServer("celltype_selectorDiffAbundance", obj,
                                               selection_mode = "hierarchical", update_selectable_nodes = selectable_nodes)

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
      calculate_differential_abundance_results(col_data_cp, sel_treelabel, node = top_sel,
                                           aggregate_by = all_of(colnames(aggr_by)),
                                           targets = vars(!!!rlang::syms(da_targets)))
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



