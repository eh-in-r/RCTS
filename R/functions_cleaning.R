#' Function with as input a dataframe. (this will be "Y" or "to_divide") It filters out rows with NA.
#'
#' @param df input
#' @importFrom stats na.omit
#' @export
handleNA <- function(df) {
  rows_with_NA <- which(apply(df, 1, function(x) sum(is.na(x)) != 0)) # rownumbers of rows containing NA
  # message(paste("There are ",length(rows_with_NA),"of",nrow(df), "rows with NA's."))

  return(list(na.omit(df), rows_with_NA))
}

#' Removes NA's in LG (in function calculate_virtual_factor_and_lambda_group() )
#' @param df input
#' @importFrom dplyr mutate_if
#' @importFrom tidyr replace_na
#' @export
handleNA_LG <- function(df) {
  result <- (as.matrix(data.frame(df) %>% mutate_if(is.numeric, replace_na, replace = 0)))
  return(result)
}

#' Function to evade floating point errors.
#'
#' Sets values that should be zero but are >0 (e.g. 1e-13) on zero.
#' @param LIMIT limit under which value is set to 0
#' @param A input
#' @export
evade_floating_point_errors <- function(A, LIMIT = 1e-13) {
  for (i in 1:length(A)) {
    if (abs(A[i]) < LIMIT) A[i] <- 0
  }
  return(A)
}
