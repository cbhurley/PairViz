\name{weighted_hpaths}
\alias{weighted_hpaths}  
\title{Constructs weight decreasing hamiltonian paths }
\description{
Returns a modified version of \code{paths}, where 
component paths/cycles are re-oriented so low weight edges occur first, and the component paths/cycles are then permuted so low-weight paths are first. 
}
\usage{
weighted_hpaths(d, path1 = NULL, paths=NULL, matrix=TRUE, cycle=NULL, 
path_weight=sum, path_dir = path_cor,...)
}
\arguments{
  \item{d}{A \code{dist}, used to provide edge weights.}
  \item{path1}{A vector giving a hamiltonian. 
 This will be the first path of the returned hamiltonian. The default is obtained from \code{order_tsp}.}
  \item{paths}{A matrix where each row is a hamiltonian. Default comes from \code{hpaths}.}
  \item{matrix}{if \code{TRUE}, returns a matrix where each row is a hamiltonian 
path, otherwise concatenates the rows into a vector. For odd \code{n}, the starting node is appended to close the eulerian.}
  \item{cycle}{If \code{TRUE}, the \code{weighted_hpaths} algorithm evaluates 
\code{path_weight} on  hamiltonian cycles, if \code{FALSE}, on open hamiltonian paths.
 Default  is \code{TRUE} for odd n and \code{FALSE} for even \code{n}.}
  \item{path_weight}{A function used  combine path weights into a single value. Default function is \code{\link{path_cor}}.}

    \item{path_dir}{A function used  to evaluate a path start and orientation.}
   \item{...}{passed to \code{path_weight}}
 }
\details{
If \code{path} is not provided, find the hamiltonian (path for even n, cycle for odd n) with the smallest total 
weight. Applying \code{path_dir} to edge weights, pick the  starting and point orientation for path1 giving the largest \code{path_dir} value. (For open paths, there are only two possible starts, for cycles there are n). 
Apply this node labelling to the  hamiltonians in the  rows of paths. 
Use criterion \code{path_dir}  again to find the best orientation for each of rows 2... of paths and
permute these rows in order of increasing \code{path_weight}.
}
\references{see \link{overview}}
\author{ C.B. Hurley and R.W. Oldford }

\seealso{\code{\link{hpaths}}, \code{\link{eulerian}}.}
\examples{

require(PairViz)

weighted_hpaths(dist(rnorm(6)))

weighted_hpaths(dist(rnorm(7)))



}

\keyword{graphs}
\keyword{optimize}
