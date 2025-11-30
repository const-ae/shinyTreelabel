
recursive_summarize <- function(data, ..., .early_stop = 0){
  stopifnot(dplyr::is.grouped_df(data))
  stopifnot(.early_stop >= 0)
  while(length(dplyr::groups(data)) > .early_stop){
    data <- dplyr::summarize(data, ..., .groups = "drop_last")
  }
  if(.early_stop == 0){
    # Do it once more
    data <- dplyr::summarize(data, ...)
  }
  data
}
