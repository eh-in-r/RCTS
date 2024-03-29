% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_IFE_algorithm_functions.R
\name{scaling_X}
\alias{scaling_X}
\title{Scaling of X.}
\usage{
scaling_X(X, firsttime, robust, vars)
}
\arguments{
\item{X}{input}

\item{firsttime}{Scaling before generating Y and before adding outliers: this is always with mean and sd. If this is FALSE, it indicates that
we are using the function for a second time, after adding the outliers. In the robust case it uses median and MAD, otherwise again mean and sd.}

\item{robust}{logical, scaling with robust metrics instead of with non-robust measures}

\item{vars}{number of observable variables}
}
\value{
3D-array with the same dimensions as X
}
\description{
Scaling of X.
}
