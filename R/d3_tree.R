treeSelectorUI <- function(id){
  tagList(
    htmltools::htmlDependency(
      name = "d3_tree",
      package = "shinyTreelabel",
      version = "0.1",
      src = "www",
      script = list(src = "d3_tree.js", type =  "module"),
      stylesheet = "d3_tree.css"
    ),
    tags$script(HTML(r"(
                        import { D3TreeSelector } from './www/d3_tree.js';
                        new D3TreeSelector("$ID$");
                      )" |> stringr::str_replace("\\$ID\\$", id)
    ), type = "module"),
    tags$div(id = NS(id, "d3tree_holder"),
             class = "d3-tree-container")
  )
}

color_generator_fnc <- function(selection_mode){
  function(selected_nodes){
    if(selection_mode == "single"){
      color_choice_fnc <- \(.) "orange"
    }else if(selection_mode %in% c("multiple", "hierarchical")){
      color_choice_fnc <- if(length(selected_nodes) > 0){
        scales::col_factor(scales::pal_hue()(length(selected_nodes)), domain = selected_nodes)
      }else{
        \(.) "orange"
      }
    }else{
      stop("Illegal selection_mode: ", selection_model)
    }
  }
}

treeSelectorServer <- function(id, spec, selection_mode = c("single", "multiple", "hierarchical"),
                               update_selectable_nodes = NULL, update_selected_node = NULL){
  selection_mode <- match.arg(selection_mode)
  if(! is.null(update_selectable_nodes)){
    stopifnot(is.reactive(update_selectable_nodes))
  }else{
    update_selectable_nodes <- reactiveVal(NULL)
  }
  if(! is.null(update_selected_node)){
    stopifnot(is.reactive(update_selected_node))
  }else{
    update_selected_node <- reactiveVal(NULL)
  }

  color_choice_fnc_gen <- color_generator_fnc(selection_mode)

  make_tree_for_javascript <- function(selected_nodes, top_selection, selectable_nodes = NULL){
    if(is.null(selectable_nodes)){
      selectable_nodes <- igraph::V(spec$tree)$name
    }

    color_choice_fnc <- color_choice_fnc_gen(selected_nodes)
    nested_list <- igraph_tree_to_nested_list(spec$tree, spec$root, callback = \(x){
      if(! is.null(top_selection) && x == top_selection){
        list(selected = TRUE, selectionColor = "orange", top_selection = TRUE, selectable = TRUE)
      }else if(x %in% selected_nodes){
        list(selected = TRUE, selectionColor = color_choice_fnc(x), selectable = TRUE)
      }else{
        list(selected = FALSE, selectable = x %in% selectable_nodes)
      }
    })
  }


  moduleServer(id, function(input, output, session){
    session$sendCustomMessage(NS(id, "firstTreeFullData"), igraph_tree_to_nested_list(spec$tree, spec$root, \(x) list(selected=FALSE)))

    selected_nodes <- character(0L)
    top_selection <- NULL

    selected_nodes_react <- reactiveVal(selected_nodes)
    top_selection_react <- reactiveVal(top_selection)

    observeEvent(input$d3TreeClick, {
      print(paste0("Clicked: ", input$d3TreeClick))
      if(selection_mode == "single"){
        selected_nodes <<- input$d3TreeClick
      }else if(selection_mode %in% c("multiple", "hierarchical")){
        if(input$d3TreeClick %in% selected_nodes){
          selected_nodes <<- remove_element_from_character_vector(selected_nodes, input$d3TreeClick)
        }else{
          selected_nodes <<- c(selected_nodes, input$d3TreeClick)
        }
      }

      if(selection_mode == "hierarchical"){
        if(! is.null(top_selection) && input$d3TreeClick == top_selection && input$d3TreeClick %in% selected_nodes){
          # This counts as deselecting the top node
          selected_nodes <<- remove_element_from_character_vector(selected_nodes, input$d3TreeClick)
          top_selection <<- NULL
        }else if(length(selected_nodes) > 1){
          top_selection <<- igraph_smallest_common_ancestor(selected_nodes, spec$tree)
        }
      }
      nested_list <- make_tree_for_javascript(selected_nodes, top_selection, update_selectable_nodes())
      session$sendCustomMessage(NS(id, "treeFullData"), nested_list)

      selected_nodes_react(selected_nodes)
      top_selection_react(top_selection)
    })

    observeEvent(update_selectable_nodes(), {
      if(! all(update_selectable_nodes() %in% igraph::V(spec$tree)$name)){
        stop("Illegal update to selectable nodes: ", toString(update_selectable_nodes()))
      }
      selected_nodes <<- intersect(selected_nodes, update_selectable_nodes())

      if(! is.null(top_selection) && ! top_selection %in% update_selectable_nodes()){
        top_selection <<- NULL
      }
      nested_list <- make_tree_for_javascript(selected_nodes, top_selection, update_selectable_nodes())
      session$sendCustomMessage(NS(id, "treeFullData"), nested_list)
      selected_nodes_react(selected_nodes)
      top_selection_react(top_selection)
    })

    observeEvent(update_selected_node(), {
      sel_nodes <- update_selected_node()
      if(is.list(sel_nodes)){
        top_selection <<- sel_nodes$top_selection
        selected_nodes <<- sel_nodes$selected_nodes
      }else{
        top_selection <<- NULL
        selected_nodes <<- sel_nodes
      }
      nested_list <- make_tree_for_javascript(selected_nodes, top_selection, update_selectable_nodes())
      session$sendCustomMessage(NS(id, "treeFullData"), nested_list)
      selected_nodes_react(selected_nodes)
      top_selection_react(top_selection)
    })

    list(selected_nodes = selected_nodes_react, top_selection = top_selection_react)
  })
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

  ancestor_list <- lapply(selected, \(n) find_ancestors(graph, n))
  Reduce(intersect, ancestor_list)[1]
}

remove_element_from_character_vector <- function(vec, element){
  stopifnot(length(element) == 1)
  pos <- which(vec == element)
  if(length(pos) == 0){
    warning("Trying to remove ", element, " but not found in vector")
  }else{
    vec[-pos]
  }
}
