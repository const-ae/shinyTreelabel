

PrecalculatedResults <- R6::R6Class("PrecalculatedResults",
  public = list(
    con = NULL,
    dbfile = NULL,

    initialize = function(dbfile = ":memory:"){
      self$dbfile <- dbfile
      self$con <- DBI::dbConnect(duckdb::duckdb(), dbdir = dbfile)
    },
    clear = function(){
      self$clear_reducedDimensions()
      self$clear_da()
      self$clear_da_meta()
      self$clear_de()
      self$clear_de_meta()
      self$clear_col_data()
    },

    add_reducedDimensions_rows = function(reducedDimensions){
      DBI::dbWriteTable(self$con, "reducedDimensions", reducedDimensions, append = TRUE)
    },
    add_da_rows = function(da_results){
      DBI::dbWriteTable(self$con, "da", da_results, append = TRUE)
    },
    add_da_meta_rows = function(da_results){
      DBI::dbWriteTable(self$con, "da_meta", da_results, append = TRUE)
    },
    add_de_rows = function(da_results){
      DBI::dbWriteTable(self$con, "de", da_results, append = TRUE)
    },
    add_de_meta_rows = function(da_results){
      DBI::dbWriteTable(self$con, "de_meta", da_results, append = TRUE)
    },
    add_col_data = function(col_data){
      treelabel_columns <- colnames(dplyr::select(col_data, where(is_treelabel)))
      expanded_col_data <- col_data |>
        mutate(across(all_of(treelabel_columns), \(x) x |> tl_score_matrix() |> as_tibble())) |>
        unpack(treelabel_columns, names_sep = "_--_")
      treelabel_df <- tibble(name = treelabel_columns) |>
        mutate(treeroot = map_chr(name, \(n) tl_tree_root(col_data[[n]])))
      igraph_df <- tibble(name = treelabel_columns) |>
        mutate(edges = map(name, \(n){
          mat <- igraph::as_edgelist(tl_tree(col_data[[n]]))
          colnames(mat) <- c("left", "right")
          as_tibble(mat)
        })) |>
        unnest(edges, names_sep = "_")

      DBI::dbWriteTable(self$con, "col_data_expanded", expanded_col_data)
      DBI::dbWriteTable(self$con, "tree_defs", igraph_df)
      DBI::dbWriteTable(self$con, "treelabel_defs", treelabel_df)
      private$col_data_cache <- NULL
    },
    clear_reducedDimensions = function(){
      DBI::dbRemoveTable(self$con, "reducedDimensions", fail_if_missing = FALSE)
    },
    clear_da = function(){
      DBI::dbRemoveTable(self$con, "da", fail_if_missing = FALSE)
    },
    clear_da_meta = function(){
      DBI::dbRemoveTable(self$con, "da_meta", fail_if_missing = FALSE)
    },
    clear_de = function(){
      DBI::dbRemoveTable(self$con, "de", fail_if_missing = FALSE)
    },
    clear_de_meta = function(){
      DBI::dbRemoveTable(self$con, "de_meta", fail_if_missing = FALSE)
    },
    clear_col_data = function(){
      DBI::dbRemoveTable(self$con, "col_data_expanded", fail_if_missing = FALSE)
      DBI::dbRemoveTable(self$con, "treelabel_defs", fail_if_missing = FALSE)
      DBI::dbRemoveTable(self$con, "tree_defs", fail_if_missing = FALSE)
      private$col_data_cache <- NULL
    },


    # Implement the obj interface
    col_data = \(){
      if(! is.null(private$col_data_cache)){
        private$col_data_cache
      }else if(all(c("col_data_expanded", "tree_defs", "treelabel_defs") %in% DBI::dbListTables(self$con))){
        treelabel_root <- dplyr::tbl(self$con, "treelabel_defs") |> dplyr::collect() |> deframe()
        igraph_df <- dplyr::tbl(self$con, "tree_defs") |> dplyr::collect()
        expanded_col_data <- dplyr::tbl(self$con, "col_data_expanded") |> dplyr::collect()
        trees <- igraph_df |>
          summarize(tree = list(igraph::graph_from_edgelist(cbind(edges_left, edges_right))), .by = name) |>
          tibble::deframe()
        for(n in names(treelabel_root)){
          sel_col <- paste0(n, "_--_")
          scores <- as.matrix(dplyr::select(expanded_col_data, starts_with(sel_col)))
          colnames(scores) <- stringr::str_remove(colnames(scores), sel_col)
          col_data[[n]] <- treelabel::treelabel(scores, tree = trees[[n]], tree_root = treelabel_root[n], propagate_up = "none")
        }
        cd <- col_data |> dplyr::select(- contains("_--_"))
        private$col_data_cache <- cd
        cd
      }
    },
    dim_reductions = \(name) {
      if("reducedDimensions" %in% DBI::dbListTables(self$con)){
        res <- dplyr::tbl(self$con, "reducedDimensions") |>
          filter(.data$name == .env$name) |>
          collect(n = 1)
        matrix(res$embedding[[1]], ncol = 2)
      }
    },
    counts = \(gene, cell) NULL,
    da_results = \(treelabel, top_selection, targets, lazy = FALSE){
      if("da" %in% DBI::dbListTables(self$con)){
        tl_test <- if(! rlang::is_missing(treelabel)){
          rlang::quo(.data$treelabel %in% .env$treelabel)
        }else{
          TRUE
        }
        top_sel_test <- if(! rlang::is_missing(top_selection)){
          rlang::quo(.data$top %in% .env$top_selection)
        }else{
          TRUE
        }
        target_test <- if(! rlang::is_missing(targets)){
          rlang::quo(.data$target %in% .env$targets)
        }else{
          TRUE
        }
        res <- dplyr::tbl(self$con, "da") |>
          dplyr::filter(!! tl_test, !! top_sel_test, !! target_test)
        if(lazy){
          res
        }else{
          dplyr::collect(res)
        }
      }
    },
    da_meta_results = \(treelabel, top_selection, targets, lazy = FALSE){
      if("da_meta" %in% DBI::dbListTables(self$con)){
        tl_test <- if(! rlang::is_missing(treelabel)){
          rlang::quo(.data$treelabel %in% .env$treelabel)
        }else{
          TRUE
        }
        top_sel_test <- if(! rlang::is_missing(top_selection)){
          rlang::quo(.data$top %in% .env$top_selection)
        }else{
          TRUE
        }
        target_test <- if(! rlang::is_missing(targets)){
          rlang::quo(.data$target %in% .env$targets)
        }else{
          TRUE
        }
        res <- dplyr::tbl(self$con, "da_meta") |>
          dplyr::filter(!! tl_test, !! top_sel_test, !! target_test)
        if(lazy){
          res
        }else{
          dplyr::collect(res)
        }
      }
    },
    de_results = \(treelabel, celltype, gene, lazy = FALSE){
      if("de" %in% DBI::dbListTables(self$con)){
        tl_test <- if(! rlang::is_missing(treelabel)){
          rlang::quo(.data$treelabel %in% .env$treelabel)
        }else{
          TRUE
        }
        celltype_test <- if(! rlang::is_missing(celltype)){
          rlang::quo(.data$celltype %in% .env$celltype)
        }else{
          TRUE
        }
        gene_test <- if(! rlang::is_missing(gene)){
          rlang::quo(.data$gene %in% .env$gene)
        }else{
          TRUE
        }
        res <- dplyr::tbl(self$con, "de") |>
          dplyr::filter(!! tl_test, !! celltype_test, !! gene_test)
        if(lazy){
          res
        }else{
          dplyr::collect(res)
        }
      }
    },
    de_meta_results = \(treelabel, celltype, gene, lazy = FALSE){
      if("de_meta" %in% DBI::dbListTables(self$con)){
        tl_test <- if(! rlang::is_missing(treelabel)){
          rlang::quo(.data$treelabel %in% .env$treelabel)
        }else{
          TRUE
        }
        celltype_test <- if(! rlang::is_missing(celltype)){
          rlang::quo(.data$celltype %in% .env$celltype)
        }else{
          TRUE
        }
        gene_test <- if(! rlang::is_missing(gene)){
          rlang::quo(.data$gene %in% .env$gene)
        }else{
          TRUE
        }

        res <- dplyr::tbl(self$con, "de_meta") |>
          dplyr::filter(!! tl_test, !! celltype_test, !! gene_test)
        if(lazy){
          res
        }else{
          dplyr::collect(res)
        }
      }
    }
  ),
  private = list(
    col_data_cache = NULL
  )
)



object_to_string <- function(obj) {
  qs::base91_encode(qs::qserialize(obj), quote_character =  "'")
}

string_to_object <- function(str) {
  qs::qdeserialize(qs::base91_decode(str))
}
