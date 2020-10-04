\name{hpaths}
\alias{zigzag}
\alias{hpaths}   
\alias{permute_hpaths}   
\title{Hamiltonian paths on the complete graph on 1..n, using Lucas-Walecki constructions.}
\description{

\code{zigzag} - Constructs hamiltonian paths where each pair (i,j) 
appears in at least one of the hamiltonians.

\code{hpaths} - Returns a hamiltonian decomposition on the complete graph with n nodes. See Details.
 
 \code{permute_hpaths} - Returns a modified version of \code{paths}, where vertices 
are re-labelled so that the first hamiltonian is \code{path1}.


}
\usage{
zigzag(n)
hpaths(n, matrix=TRUE,cycle=NULL,...)
permute_hpaths(path1,paths= hpaths(length(path1)), matrix=TRUE,...)
}
\arguments{
  \item{n}{a positive integer. For \code{hpaths}, \code{n} may also be a vector specifying the first hamiltonian.}
  \item{matrix}{if \code{TRUE}, returns a matrix where each row is a hamiltonian 
path, otherwise concatenates the rows into a vector.}
  \item{cycle}{If \code{TRUE}, returns hamiltonian cycles, i.e. every hamiltonian starts at the same node.
If \code{FALSE}, returned paths are open. Defaults to \code{TRUE} for odd n,\code{FALSE} for even n. }
  \item{path1}{A vector- This  will be the first hamiltonian of the returned hamiltonian decomposition.}
  \item{paths}{A matrix where each row is a hamiltonian.}
  \item{...}{Ignored.}
 }
\value{A numeric matrix where each row  contains a permutation of 1..n, or these rows concatenated into a vector if \code{matrix=FALSE}.}
\details{

\code{hpaths} - 
From graph theory we know that for odd n, the complete graph decomposes into (n-1)/2 edge distinct hamiltonian
cycles, while for even n the graph decomposes into n/2 edge distinct hamiltonian paths.
The default behaviour of the function  \code{hpaths} is to produce the cycle decomposition for odd n and the path
decomposition for even n.

However, if a \code{TRUE} value is supplied as argument cycle, the returned paths are cycles, and the result is a true  decomposition for odd n, but for even n the last hamiltonian has some duplicate edges.
 If a \code{FALSE} value is supplied as argument cycle, the returned paths are open, and the result is a true decomposition 
  for even n, but for odd n the last hamiltonian has some duplicate edges. 
}

\references{
	D.E. Lucas (1892), Recreations Matematiques, Vol II. Gauthier Villars, Paris.
	
Also	see \link{overview}}

\author{ C.B. Hurley and R.W. Oldford }

\seealso{\code{\link{weighted_hpaths}}, \code{\link{eseq}}.}
\examples{
require(PairViz)

zigzag(7)
hpaths(7) # the rows form a decomp. into hamiltonian cycles

# Now concatenate the rows and close the path
hpaths(7,matrix=FALSE)

# Form a decomposition into hamiltonian cycles- 
# this decomposition is not exact, as the last row duplicates edges
hpaths(7,cycle=FALSE)

# For even n, the default is a decomposition into hamiltonian paths, not cycles.
hpaths(6)

# If cycles are required for even n, 
# the decomposition will not be exact and the last row duplicates edges

hpaths(6,cycle=TRUE)

# If you want to specify the first hamiltonian of the decomposition, use
hpaths(1:7)


}

