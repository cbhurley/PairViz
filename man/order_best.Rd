\name{order_best}
\alias{order_best}
\title{Uses brute-force enumeration  to find the best hamiltonian on the complete graph on 1..n.  }
\description{
Returns the best hamiltonian
}
\usage{
order_best(d, maxexact=9,nsamples=50000,path_weight=sum,
cycle=FALSE,path_dir = path_cor,...)
}
\arguments{
  \item{d}{A \code{dist}, used to provide edge weights.}
   \item{maxexact}{If the sequence length is  \code{<= maxexact}, find the overall best hamiltonian, otherwise compares \code{nsamples}
randomly generated permutations.}
    \item{nsamples}{If the sequence length is  \code{<= maxexact}, finds the  best of \code{nsamples} randomly generated permutations .}
  \item{cycle}{If \code{TRUE}, finds the shortest cycle, otherwise the shortest open path.}
  \item{path_weight}{Combines edge weights into a single path/cycle weight.}
  \item{path_dir}{If a function is provided, used to re-orient the cycle/path. Default function is \code{\link{path_cor}}.}
  \item{...}{Ignored.}
}
\value{A vector containing a permutation of 1..n}
\details{Requires package gtools. 
Currently it is possible to find the best hamiltonian by complete enumeration for up to 10 nodes.
When  \code{path_dir}  is non NULL, the returned hamiltonian is also optimally oriented using \code{best_orientation}, which compares orientations via \code{path_dir}.}
\references{see \link{overview}}
\author{ C.B. Hurley and R.W. Oldford }

\seealso{\code{\link{order_tsp}}.}
\examples{
require(PairViz)
order_best(eurodist)



}

\keyword{graphs}
\keyword{optimize}
