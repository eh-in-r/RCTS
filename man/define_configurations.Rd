% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_IFE_algorithm_functions.R
\name{define_configurations}
\alias{define_configurations}
\title{Constructs dataframe where the rows contains all configurations that are included and for which the estimators will be estimated.}
\usage{
define_configurations(S_cand, k_cand, kg_cand)
}
\arguments{
\item{S_cand}{candidates for S (number of groups)}

\item{k_cand}{candidates for k (number of common factors)}

\item{kg_cand}{candidates for kg (number of group specific factors)}
}
\value{
data.frame
}
\description{
Constructs dataframe where the rows contains all configurations that are included and for which the estimators will be estimated.
}
\examples{
define_configurations(2:4, 0, 2:3)
}
