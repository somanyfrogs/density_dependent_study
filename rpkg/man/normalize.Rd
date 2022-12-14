% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_common.R
\name{normalize}
\alias{normalize}
\title{Normalization function}
\usage{
normalize(x, method = "0-1")
}
\arguments{
\item{x}{Numeric vector, to be normalized}

\item{method}{Character, which switches the normalization method
if method == "0-1", the x is normalized to range between 0-1
otherwise, x is normalized with mean 0 and unit variance}
}
\description{
Normalization function
}
