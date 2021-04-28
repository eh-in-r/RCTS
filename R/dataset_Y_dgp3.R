


#' Y_dgp3 contains a simulated dataset for DGP 3.
#'
#' Y = XB + LgFg. It has 3 groups and each group has 3 groupfactors.
#' At last there were 3 observable variables generated into it.
#' @format 300 x 30 matrix. Each row is one time series.
#' @examples
#' plot(Y_dgp3[,1:2], col = g_true_dgp3, xlab = "First column of Y",  ylab = "Second column of Y",
#' main = "Plot of the first two columns of the dataset Y. \nColors are the true groups.")
"Y_dgp3"
