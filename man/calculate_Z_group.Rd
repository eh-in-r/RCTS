% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_IFE_algorithm_functions.R
\name{calculate_Z_group}
\alias{calculate_Z_group}
\title{Calculates Z = Y - X*beta_est - LF. It is used to estimate the groupfactorstructure.}
\usage{
calculate_Z_group(
  Y,
  X,
  beta_est,
  g,
  lambda,
  comfactor,
  group,
  k,
  method_estimate_beta,
  initialise,
  vars_est
)
}
\arguments{
\item{Y}{Y: NxT dataframe with the panel data of interest}

\item{X}{X: NxTxp array containing the observable variables}

\item{beta_est}{estimated values of beta}

\item{g}{Vector with group membership for all individuals}

\item{lambda}{loadings of the estimated common factors}

\item{comfactor}{estimated common factors}

\item{group}{indexnumber of the group}

\item{k}{number of common factors to be estimated}

\item{method_estimate_beta}{defines how beta is estimated. Default case is an estimated beta for each individual. Default value is "individual." Possible values are "homogeneous", "group" or "individual".}

\item{initialise}{indicator of being in the initialisation phase}

\item{vars_est}{number of variables that will be included in the algorithm and have their coefficient estimated. This is usually equal to the number of observable variables.}
}
\value{
NxT matrix
}
\description{
Calculates Z = Y - X*beta_est - LF. It is used to estimate the groupfactorstructure.
}
