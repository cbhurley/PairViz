\name{path_weights}
\alias{path_weights}
\alias{path_cis}
\alias{edge2dist}
\alias{dist2edge}
\alias{edge_index}
\alias{nnodes}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Utility functions to manipulate pairwise information.
}
\description{
  These  functions perform calculations on edge matrices containing pairwise information.}
\usage{
path_weights(edgew, path,  symmetric = TRUE,edge.index=edge_index(edgew),...)
path_cis(edgew, path,edge.index=edge_index(edgew),ci.pos=FALSE)
edge2dist(edgew, edge.index=edge_index(edgew))
dist2edge(d)
edge_index(x, order="default")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edgew}{A Matrix (or vector) whose ith row (or element) has weights for pair indexed by pair in row i of edge.index.
  For \code{edge2dist}, \code{edgew} should be a vector.}
  \item{path}{Vector of indices into rows of \code{edgew}.}
  \item{symmetric}{If \code{TRUE} edge weights are interpreted as symmetric. }
   \item{edge.index}{A 2-column matrix with each row giving indices for
 corresponding weight in \code{edgew}.}
  \item{ci.pos}{If TRUE, all CIs are mu(max) - mu(min), otherwise mu(right) - mu(left). }
   \item{d}{A \code{dist} or matrx of distances.}
    \item{order}{If "low.order.first" or "scagnostics", lists lowest index pairs first, otherwise lists pairs starting with 1, then 2 etc.}
  \item{x}{An edgew matrix or vector, or a positive integer.}
  \item{...}{Ignored}

}
\details{
 \code{path_weights} - Returns matrix of path weights so that the ith row of result contains weights for indices path[i], path[i+1]
 
  \code{path_cis} - Returns matrix of path confidence intervals so that the ith row of result contains intervals for mean-path[i] -  mean-path[i+1]
  
 \code{edge2dist} - Returns a \code{dist},
 containing elements of \code{edgew}. 
 
  \code{dist2edge} - Returns a vector of edge weights.
 
 \code{edge_index} -A generic function. Returns a 2-column matrix with one row for
 each edge.  Each row contains an index pair i,j. If \code{order} is "low.order.first" or "scagnostics", lists lowest index pairs first - this is the default ordering for class \code{scagdf}, otherwise lists pairs
  starting with 1, then 2 etc
  
  \code{nnodes} - Here \code{edgew} contains edge weights for a complete graph;  returns the number of nodes in this complete graph.

 }
\references{see \link{overview}}

\author{ C.B. Hurley and R.W. Oldford }

\examples{
require(PairViz)

s <- matrix(1:40,nrow=10,ncol=4)

edge2dist(s[,1])

path_weights(s,1:4)
path_weights(s,eseq(5))

fm1 <- aov(breaks ~ wool + tension, data = warpbreaks)
tuk <- TukeyHSD(fm1, "tension")[[1]]

# Here the first argument (weight matrix) can have number of columns

path_weights(tuk,c(1:3,1))  



# Here the first argument (weight matrix) should have an odd number of columns-
# the first is the mean difference, other column pairs are endpoints of CIs


path_cis(tuk[,-4],c(1:3,1))
 
 
 }
 
\keyword{graphs}
