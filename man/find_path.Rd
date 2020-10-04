\name{find_path}
\alias{find_path}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Constructs a path from a matrix of edge weights.
}
\description{
 Returns a path,  constructed by applying the function in \code{path}
  to the edge weights. If each edge has many weights, i.e if \code{edgew} is a matrix, these weights are first reduced by the function \code{combine} applied to the rows. If \code{path} is NULL, the returned path defaults to 1..\code{nnodes(edgew)}}
  
\usage{
find_path(edgew, path=NULL,  combine=sum, edge.index=edge_index(edgew),...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edgew}{Matrix (or vector) whose ith row (or element) has weights for pair indexed by pair in row i of edge.index.}
  \item{path}{a function
  	used to construct the index path.}
      \item{combine}{A function that combines the row of weights for an edge into a single numeric value.}
\item{edge.index}{A 2-column matrix with each row giving indices for
 corresponding weight in \code{edgew}.}

  \item{...}{passed to path construction function.}

}

\author{ C.B. Hurley and R.W. Oldford }

 
\keyword{graphs}
\keyword{optimize}
