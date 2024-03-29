% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_IFE_algorithm_functions.R
\name{handle_macropca_errors}
\alias{handle_macropca_errors}
\title{Helpfunction in robustpca().}
\usage{
handle_macropca_errors(
  object,
  temp,
  KMAX,
  number_eigenvectors,
  verbose = FALSE
)
}
\arguments{
\item{object}{input}

\item{temp}{this is the result of the trycatch block of using macropca on object}

\item{KMAX}{parameter kmax in MacroPCA}

\item{number_eigenvectors}{number of principal components that are needed}

\item{verbose}{when TRUE, it prints messages}
}
\value{
matrix of which the columns contain the chosen amount of eigenvectors of object
}
\description{
It handles possible thrown errors in MacroPCA.
}
