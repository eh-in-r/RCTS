% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_IFE_algorithm_functions.R
\name{create_data_dgp2}
\alias{create_data_dgp2}
\title{Creates an instance of DGP 2, as defined in \insertCite{BoudtHeyndels2021;textual}{RCTS}.}
\usage{
create_data_dgp2(N, TT, S_true = 3, vars = 3, k_true = 0, kg_true = c(3, 3, 3))
}
\arguments{
\item{N}{number of time series}

\item{TT}{length of time series}

\item{S_true}{true number of groups}

\item{vars}{number of available observed variables}

\item{k_true}{true number of common_factors}

\item{kg_true}{vector with the true number of group factors for each group}
}
\value{
list
}
\description{
The default has 3 groups with each 3 group specific factors. Further it contains 0 common factors and 3 observed variables.
The output is a list where the first element is the simulated panel dataset (a dataframe with N (amount of time series) rows and T (length of time series) columns).
The second element contains the NxTxp array with the p observed variables. The third element contains the true group membership.
The fourth element contains the true beta's (this has p+1 rows and one column for each group).
The fifth element contains a list with the true group specific factors.
The sixth element contains a dataframe with N rows where each row contains the group specific factor loadings that corresponds to the group specific factors.
Further it contains the true group membership and an index (this corresponds to the rownumber in Y and X).
The seventh and eighth elements contain the true common factor(s) and its loadings respectively.
}
\examples{
create_data_dgp2(30, 10)
}
