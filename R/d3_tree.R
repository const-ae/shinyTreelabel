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

treeSelectorServer <- function(id){
  moduleServer(id, function(input, output, session){
    session$sendCustomMessage(NS(id, "firstTreeFullData"), igraph_tree_to_nested_list(.vals$tree, .vals$root, \(x) list(selected=FALSE)))
    observeEvent(input$d3TreeClick, {
      clicked_node <- input$d3TreeClick
      print(paste0("Clicked: ", clicked_node))
      nested_list <- igraph_tree_to_nested_list(.vals$tree, .vals$root, callback = \(x){
        if(x == clicked_node){
          list(selected = TRUE, selectionColor = "orange")
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
             tabPanel("Differential Abundance"),
             tabPanel("Gene-level analysis")
  )
}

d3_tree_server <- function(input, output, session){
  treeSelectorServer("umap_selector")
  treeSelectorServer("de_selector")
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

