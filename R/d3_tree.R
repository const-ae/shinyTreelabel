treeSelectorUI <- function(id){
  tagList(
    htmltools::htmlDependency(
      name = "d3_tree",
      package = "shinyTreelabel",
      version = "0.1",
      src = "../inst/www",
      script = list(src = "d3_tree.js", type =  "module")
    ),
    tags$script(HTML(r"(
                        import { D3TreeSelector } from './www/d3_tree.js';
                        new D3TreeSelector("$ID$");
                      )" |> stringr::str_replace("\\$ID\\$", id)
    ), type = "module"),
    tags$div(id = NS(id, "d3tree_holder"))
  )
}

treeSelectorServer <- function(id, selection_mode = c("single", "multiple", "hierarchical")){
  selection_mode <- match.arg(selection_mode)

  moduleServer(id, function(input, output, session){
    session$sendCustomMessage(NS(id, "firstTreeFullData"), igraph_tree_to_nested_list(.vals$tree, .vals$root, \(x) list(selected=FALSE)))

    selected_nodes <- new.env(parent = emptyenv())
    top_selection <- NULL

    observeEvent(input$d3TreeClick, {
      print(paste0("Clicked: ", input$d3TreeClick))
      if(selection_mode == "single"){
        rm(list = ls(envir = selected_nodes), envir = selected_nodes)
        selected_nodes[[input$d3TreeClick]] <- TRUE
        color_choice_fnc <- \(.) "orange"
      }else if(selection_mode %in% c("multiple", "hierarchical")){
        if(exists(input$d3TreeClick, envir = selected_nodes)){
          rm(list = input$d3TreeClick, envir = selected_nodes)
        }else{
          selected_nodes[[input$d3TreeClick]] <- TRUE
        }
        sel <- ls(envir = selected_nodes)
        color_choice_fnc <- if(length(sel) > 0){
          scales::col_factor(scales::pal_hue()(length(sel)), domain = sel)
        }
      }

      if(selection_mode == "hierarchical"){
        if(! is.null(top_selection) && input$d3TreeClick == top_selection && exists(input$d3TreeClick, envir = selected_nodes)){
          # This counts as deselecting the top node
          rm(list = input$d3TreeClick, envir = selected_nodes)
          sel <- ls(envir = selected_nodes)
          color_choice_fnc <- if(length(sel) > 0){
            scales::col_factor(scales::pal_hue()(length(sel)), domain = sel)
          }
          top_selection <<- NULL
        }else if(length(sel) > 1){
            top_selection <<- igraph_smallest_common_ancestor(sel, .vals$tree)
        }
      }

      nested_list <- igraph_tree_to_nested_list(.vals$tree, .vals$root, callback = \(x){
        if(! is.null(top_selection) && x == top_selection){
          list(selected = TRUE, selectionColor = "orange", top_selection = TRUE)
        }else if(exists(x, envir = selected_nodes)){
          list(selected = TRUE, selectionColor = color_choice_fnc(x))
        }else{
          NULL
        }
      })
      session$sendCustomMessage(NS(id, "treeFullData"), nested_list)
    })
  })
}

d3_tree_ui <- function(){
  resources <- system.file("www", package = "shinyTreelabel")
  navbarPage(title = "Explore your single cell data",
             tabPanel("Overview (UMAP etc.)",
                      treeSelectorUI('umap_selector')),
             tabPanel("Differential Expression",
                      selectInput("treelabelSelector", "Treelabel selector", choices = letters, selected = "A"),
                      treeSelectorUI('de_selector')
             ),
             tabPanel("Differential Abundance",
                      treeSelectorUI('da_selector')),
             tabPanel("Gene-level analysis")
  )
}

d3_tree_server <- function(input, output, session){
  treeSelectorServer("umap_selector", selection_mode = "multiple")
  treeSelectorServer("de_selector")
  treeSelectorServer("da_selector", selection_mode = "hierarchical")
}

igraph_tree_to_nested_list <- function(tree, root = "root", callback = \(x) NULL){
  all_names <- igraph::V(tree)$name
  hash_map <- new.env(parent = emptyenv(), size = length(all_names))
  for(n in all_names){
    hash_map[[n]] <- new.env(parent = emptyenv(), size = length(all_names))
  }
  result <- new.env(parent = emptyenv(), size = length(all_names))
  result[[root]] <- hash_map[[root]]
  for(v in igraph::V(tree)$name){
    node <- hash_map[[v]]
    for(c in igraph::neighbors(tree, v, mode = "out")$name){
      node[[c]] <- hash_map[[c]]
    }
  }
  nested_env_list <- function(x, name){
    if(is.environment(x)){
      out <- as.list(x)
      child <- lapply(seq_along(out), \(idx) nested_env_list(out[[idx]], names(out)[idx]))
      c(list(name = name, children = child), callback(name))
    }else{
      c(list(name = name, children = list()), callback(name))
    }
  }
  nested_env_list(result[[root]], root)
}

igraph_smallest_common_ancestor <- function(selected, graph){
  find_ancestors <- function(graph, node) {
    ancestors <- c()
    parent <- igraph::neighbors(graph, node, mode = "in")
    while (length(parent) > 0) {
      ancestors <- c(ancestors, parent)
      parent <- igraph::neighbors(graph, parent, mode = "in")
    }
    c(node, names(ancestors))
  }

  ancestor_list <- lapply(selected, \(n) find_ancestors(tree, n))
  Reduce(intersect, ancestor_list)[1]
}

