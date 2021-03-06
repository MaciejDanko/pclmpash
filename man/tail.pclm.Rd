% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pclm_pash.r
\name{tail.pclm}
\alias{tail.pclm}
\title{Tail function for PCLM object}
\usage{
tail.pclm(object, n = 6L, type = c("lt", "aggregated", "nonaggregated"))
}
\arguments{
\item{object}{PCLM object.}

\item{n}{A single integer. If positive, size for the resulting object: number of rows for a life-table. If negative, all but the n last/first number of elements of x.}

\item{type}{which life-table  should be returned. One of \code{"lt"}, \code{"aggregated"} or \code{"nonaggregated"}.}
}
\description{
\emph{bold{Generic function}}
}
\author{
Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
}

