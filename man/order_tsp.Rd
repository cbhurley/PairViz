\name{order_tsp}
\alias{order_tsp}  
\title{Uses tsp to find the best hamiltonian on the complete graph on 1..n}
\description{
Returns shortest cycle or path via tsp solver from package TSP
}
\usage{
order_tsp(d, method = "nearest", cycle=FALSE,improve=FALSE,path_dir = path_cor,...)
}
\arguments{
  \item{d}{A \code{dist}, used to provide edge weights.}
  \item{method}{Options are \code{nearest_insertion}, \code{farthest_insertion}, \code{cheapest_insertion},         
	\code{arbitrary_insertion}, \code{nn}, \code{repetitive_nn},
	\code{2-opt} and if concorde package is loaded, \code{concorde}. See \code{\link[TSP]{solve_TSP}} for details.
 }
  \item{improve}{if \code{TRUE}, attempts to improve the solution using "2-opt".}
  \item{cycle}{If \code{TRUE}, finds the shortest cycle, otherwise the shortest open path.}
  \item{path_dir}{If a function is provided, used to re-orient the cycle/path. Default function is \code{\link{path_cor}}.}
  \item{...}{passed to \code{solve_tsp}}
}
\value{A vector containing a permutation of 1..n}
\details{Requires package TSP.
When  \code{path_dir}  is non NULL, the returned hamiltonian is also optimally oriented using \code{best_orientation}, which compares orientations via \code{path_dir}.}
\references{See package TSP.}
\author{ C.B. Hurley and R.W. Oldford }

\seealso{\code{\link{order_best}}, \code{\link[TSP]{solve_TSP}} in \pkg{TSP}.}
\examples{
require(PairViz)


rdist <- function(n) {
	d <- matrix(0,n,n)
	d[lower.tri(d)] <- runif(n*(n-1)/2)
	return(as.dist(d))
	}


order_tsp(rdist(7))


%order_tsp(eurodist) fails, because tour_length(TSP(eurodist)) fails, mysteriously!

edist <- as.dist(as.matrix(eurodist))
order_tsp(edist)

}

\keyword{graphs}
\keyword{optimize}
