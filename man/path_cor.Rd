\name{path_cor}
\alias{path_cor}
\title{ Measures the tendency of edge weights to increase.}
\description{
  Returns the (Kendalls tau) correlation of the edge weights with the vector 1.. (number of weights).
}
\usage{
path_cor(edgew, method = "kendall")
}
\arguments{
  \item{edgew}{ A vector of edge weights. }
  \item{method}{ passed to \code{\link{cor}}}
}

