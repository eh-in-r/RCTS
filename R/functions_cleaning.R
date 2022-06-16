#' Function with as input a dataframe. (this will be "Y" or "to_divide") It filters out rows with NA.
#'
#' @param df input
#' @importFrom stats na.omit
#' @return list with a dataframe where the rows with NA are filtered out, and a dataframe with only those rows
handleNA <- function(df) {
  rows_with_NA <- which(apply(df, 1, function(x) sum(is.na(x)) != 0)) # rownumbers of rows containing NA
  # message(paste("There are ",length(rows_with_NA),"of",nrow(df), "rows with NA's."))

  return(list(na.omit(df), rows_with_NA))
}

#' Removes NA's in LG (in function calculate_virtual_factor_and_lambda_group() )
#' @param df input
#' @importFrom dplyr mutate_if
#' @importFrom tidyr replace_na
#' @return matrix
handleNA_LG <- function(df) {
  result <- (as.matrix(data.frame(df) %>% mutate_if(is.numeric, replace_na, replace = 0)))
  return(result)
}

#' Function to evade floating point errors.
#'
#' Sets values that should be zero but are >0 (e.g. 1e-13) on zero.
#' @param LIMIT limit under which value is set to 0
#' @param x numeric input
#' @return numeric
evade_floating_point_errors <- function(x, LIMIT = 1e-13) {
  for (i in 1:length(x)) {
    if (abs(x[i]) < LIMIT) x[i] <- 0
  }
  return(x)
}
