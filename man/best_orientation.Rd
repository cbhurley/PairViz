\name{best_orientation}
\alias{best_orientation}       
\title{Re-orients a path to be weight-decreasing }
\description{

Re-orients a path/cycle, preserving adjacencies so that weights tend to decrease. From specifies the starting point, for cycles only.
}
\usage{
best_orientation(path, d, cycle=FALSE, path_dir= path_cor, from=NULL)
}
\arguments{
  \item{path}{A vector giving a hamiltonian. } 	
    \item{d}{A \code{dist}, used to provide edge weights.}

  \item{cycle}{If \code{TRUE}, the \code{path} is interpreted as a closed path.}
     \item{path_dir}{A function used  to evaluate a path start and orientation}
     \item{from}{Sepcifies the starting point, for cycles only.}
 }
\references{see \link{overview}}

\author{ C.B. Hurley and R.W. Oldford }

\seealso{\code{\link{hpaths}}, \code{\link{eulerian}}.}
\examples{

require(PairViz)

rdist <- function(n) {
	d <- matrix(0,n,n)
	d[lower.tri(d)] <- runif(n*(n-1)/2)
	return(as.dist(d))
	}
r <- rdist(7)
best_orientation(1:7,r) 
best_orientation(1:7,r,cycle=TRUE) 

}

\keyword{graphs}
\keyword{optimize}
