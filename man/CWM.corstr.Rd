% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Funs.R
\name{CWM.corstr}
\alias{CWM.corstr}
\title{Construct the correlation working matrix}
\usage{
CWM.corstr(rho, ndim, corstr)
}
\arguments{
\item{rho}{the common parameter needed to construct the correlation working matrix}

\item{ndim}{the dimension of the matrix that we will produce}

\item{corstr}{a character string specifying the correlation structure.
The following are permitted: \code{independence}, \code{exchangeable} and \code{ar1}.}
}
\description{
Construct the correlation working matrix
}
