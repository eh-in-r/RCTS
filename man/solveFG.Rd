% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_IFE_algorithm_functions.R
\name{solveFG}
\alias{solveFG}
\title{Helpfunction in update_g(), to calculate solve(FG x t(FG)) x FG}
\usage{
solveFG(TT, S, kg, factor_group, testing = FALSE)
}
\arguments{
\item{TT}{length of time series}

\item{S}{number of groups}

\item{kg}{vector with the estimated number of group specific factors for each group}

\item{factor_group}{estimated group specific factors}

\item{testing}{variable that determines if we are in 'testing phase'; defaults to FALSE (requires Matrix-package if set to TRUE)}
}
\value{
list: the number of elements in this list is equal to S (the number of groups).
Each of the elements in this list has a number rows equal to the number of group specific factors, and TT columns.
}
\description{
Helpfunction in update_g(), to calculate solve(FG x t(FG)) x FG
}
