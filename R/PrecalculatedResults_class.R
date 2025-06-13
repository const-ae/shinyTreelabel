

PrecalculatedResults <- R6::R6Class("PrecalculatedResults",
  public = list(
    con = NULL,
    dbfile = NULL,

    initialize = function(dbfile = ":memory:"){
      self$dbfile <- dbfile
      self$con <- DBI::dbConnect(duckdb::duckdb(), dbdir = dbfile)
    },
    clear = function(){
      DBI::dbRemoveTable(self$con, "da", fail_if_missing = FALSE)
      DBI::dbRemoveTable(self$con, "da_meta", fail_if_missing = FALSE)
      DBI::dbRemoveTable(self$con, "de", fail_if_missing = FALSE)
      DBI::dbRemoveTable(self$con, "de_meta", fail_if_missing = FALSE)
      DBI::dbRemoveTable(self$con, "col_data", fail_if_missing = FALSE)
      DBI::dbRemoveTable(self$con, "reducedDimensions", fail_if_missing = FALSE)
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
      DBI::dbWriteTable(self$con, "col_data", data.frame(x = object_to_string(col_data)),
                        overwrite = TRUE)
      private$col_data_cache <- NULL
    },


    # Implement the obj interface
    col_data = \(){
      if(! is.null(private$col_data_cache)){
        private$col_data_cache
      }else if("col_data" %in% DBI::dbListTables(self$con)){
        cd <- tbl(self$con, "col_data") |>
          dplyr::collect(n = 1)
        cd <- string_to_object(cd$x)
        private$col_data_cache <- cd
        cd
      }else{
        NULL
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
  rawToChar(serialize(obj, NULL, ascii = TRUE))
}

string_to_object <- function(str) {
  unserialize(charToRaw(str))
}
