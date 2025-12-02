


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
    levels <- dplyr::last(spec$metaanalysis_levels)
    selectInput("metaanalysisSelector", "Display results for:",
                choices = c(if(length(levels >= 2)) "Meta-Analysis", levels), selected = levels[1])
  }

  shiny_theme <- bslib::bs_theme(
    version = 5,
    bootswatch = "flatly"
    # base_font = bslib::font_google("Inter"),
    # heading_font = bslib::font_google("Inter")
  )


  navbarPage(title = "Explore your single cell data",
             theme = shiny_theme, #shinythemes::shinytheme("yeti"),
             tabPanel("Overview (UMAP etc.)",
                      selectInput("treelabelSelector0", "Treelabel selector", choices = spec$treelabel_names,
                                  selected = spec$treelabel_names[1], multiple = FALSE),
                      p("Click on any cell type to see it's location on the 2D plot."),
                      treeSelectorUI('celltype_selectorUMAP'),
                      selectInput("redDimSelector", "Dimension reduction", choices = spec$dim_reduction_names),
                      plotOutput(outputId = "redDimPlot")
             ),
             tabPanel("Differential Expression",
                      selectInput("treelabelSelector1", "Treelabel selector", choices = spec$treelabel_names,
                                  selected = spec$treelabel_names[1], multiple = FALSE),
                      p("Click on any cell type to see the top differentially expressed genes."),
                      treeSelectorUI('celltype_selectorDEView'),
                      metaanalysis_selector,
                      div(
                        style = "max-width: 1400px; margin: auto;",
                        shinycssloaders::withSpinner(plotOutput(outputId = "deVolcano")),
                        div(style = "width=100%; overflow-x: auto;",
                          shinycssloaders::withSpinner(DT::DTOutput(outputId = "deTable"))
                        )
                      )
             ),
             tabPanel("Differential Abundance",
                      selectInput("treelabelSelector2", "Treelabel selector", choices = spec$treelabel_names,
                                  selected = spec$treelabel_names[1], multiple = FALSE),
                      p("Click any two cell types to calculate their differential abundance."),
                      p("If the nodes are siblings, the best top level reference is determined. The reference is always colored yellow."),
                      treeSelectorUI('celltype_selectorDiffAbundance'),
                      div(
                        style = "max-width: 1400px; margin: auto;",
                        shinycssloaders::withSpinner(plotOutput(outputId = "diffAbundancePlotsOverview")),
                        shinycssloaders::withSpinner(uiOutput("diffAbundancePlotsFullPlaceholder"))
                      )
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

    umap_col_data <- obj$col_data() |>
      mutate(reddim = reduced_dim())
    if(length(setdiff(sel_nodes, spec$root)) == 0){
      umap_col_data <- umap_col_data |> mutate(coloring = spec$root)
    }else{
      umap_col_data <- umap_col_data |>
        mutate(coloring = tl_name(tl_tree_filter(!! rlang::sym(input$treelabelSelector0), \(x) setdiff(sel_nodes, spec$root))))
    }
    umap_col_data |>
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

    sel_nodes <- cellTypeSelectorDAView$selected_nodes()
    color_map <- structure(color_choice_fnc_gen(sel_nodes)(sel_nodes), names = sel_nodes)

    da_meta_res() |>
      mutate(conf.lower = LFC - qnorm(0.975) * LFC_se,
             conf.high = LFC + qnorm(0.975) * LFC_se) |>
      mutate(includes_zero = conf.lower < 0 & conf.high > 0) |>
      mutate(target = factor(target, levels = sel_nodes)) |>
      ggplot(aes(x = LFC, y = target)) +
        geom_vline(xintercept = 0) +
        geom_pointrange(aes(xmin = conf.lower, xmax = conf.high, color = target, alpha = includes_zero), linewidth = 1.5, fatten = 4, show.legend = FALSE) +
        scale_color_manual(values = color_map, drop = FALSE) +
        scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.5)) +
        scale_x_continuous(limits = c(-3.2, 3.2), expand = expansion(add = 0), oob = scales::oob_squish) +
        facet_wrap(vars(contrast))
  })

  output$diffAbundancePlotsFull <- renderPlot({
    req(da_results())
    req(da_meta_res())

    sel_nodes <- cellTypeSelectorDAView$selected_nodes()
    color_map <- structure(color_choice_fnc_gen(sel_nodes)(sel_nodes), names = sel_nodes)

    da_results() |>
      bind_rows(da_meta_res() |> mutate(..meta = "Meta")) |>
      mutate(..meta = forcats::fct_relevel(..meta, "Meta", after = 0)) |>
      mutate(..meta = forcats::fct_rev(..meta)) |>
      mutate(conf.lower = LFC - qnorm(0.975) * LFC_se,
             conf.high = LFC + qnorm(0.975) * LFC_se) |>
      mutate(includes_zero = conf.lower < 0 & conf.high > 0) |>
      mutate(target = forcats::fct_rev(factor(target, levels = sel_nodes))) |>
      ggplot(aes(x = LFC, y = ..meta)) +
        geom_vline(xintercept = 0) +
        geom_pointrange(aes(xmin = conf.lower, xmax = conf.high, color = target, alpha = includes_zero), linewidth = 1.5, fatten = 5, show.legend = FALSE) +
        scale_color_manual(values = color_map, drop = FALSE) +
        scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.5)) +
        facet_grid(vars(target), vars(contrast)) +
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
        filter(..meta == input$metaanalysisSelector)
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
    req(nrow(de_result_display()) > 0)

    de_result_display() |>
      mutate(adj_pval = p.adjust(pval, method = "BH")) |>
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




