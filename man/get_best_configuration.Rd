% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_IFE_algorithm_functions.R
\name{get_best_configuration}
\alias{get_best_configuration}
\title{Finds the first stable interval after the first unstable point. It then defines the value for C for the begin, middle and end of this interval.}
\usage{
get_best_configuration(
  list_vc,
  list_rc,
  list_rcj,
  C_candidates,
  S_cand,
  return_short = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{list_vc}{list with resulting expression(VC^2) for each run}

\item{list_rc}{list with resulting rc for each run}

\item{list_rcj}{list with resulting rcj for each run}

\item{C_candidates}{candidates for C}

\item{S_cand}{candidates for S (number of groups)}

\item{return_short}{if TRUE, the function returns the dataframe filtered for several specified potential candidates for C}

\item{verbose}{when TRUE, it prints messages}
}
\value{
data.frame with the optimized configuration for each candidate C (if return_short is FALSE) and
for each of the selected C's in the chosen stable interval (if return_short is TRUE).
}
\description{
Finds the first stable interval after the first unstable point. It then defines the value for C for the begin, middle and end of this interval.
}
\examples{
set.seed(1)
all_best_values <- calculate_best_config(add_configuration(initialise_df_results(TRUE),
  3, 0, c(3, 3, 3, rep(NA, 17))),
  data.frame(t(1:5)), 1:5)
rc <- fill_rc(initialise_rc(0:1, 1:5), all_best_values, 0)
rc <- fill_rc(rc, all_best_values, 1)
rcj <- fill_rcj(initialise_rcj(0:1, 1:5) , all_best_values, 0, 2:4, 2:4)
rcj <- fill_rcj(rcj, all_best_values, 1, 2:4, 2:4)
get_best_configuration(sort(runif(5)), rc, rcj, 1:5, 2:4, return_short = FALSE)
}
