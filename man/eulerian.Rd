\name{eulerian}
\docType{methods}
\alias{eulerian}
\alias{eulerian-methods}
\alias{eulerian,ANY-method}
\alias{eulerian,graphNEL-method}
\alias{eulerian,even_graph-method}
\alias{eulerian,matrix-method}
\alias{eulerian,numeric-method}

\title{ ~~ Methods for Function eulerian  ~~}
\description{
  A generic function that returns an eulerian (or nearly eulerian)
 path based on \code{self}.
}

\usage{
eulerian(self, start=NULL,weighted=TRUE)
}

\arguments{
  \item{self}{-- see below}
  \item{start}{-- see below}
  \item{weighted}{-- see below}
}

\value{A vector representing the eulerian- a character vector of node names for a graph, otherwise a numeric vector. If the graph is not connected, the result is a list of eulerians for each connected component.}
	
\section{Methods}{
\describe{
\item{self = "even_graph"}{Uses \code{etour} to construct the eulerian. If \code{weighted} is TRUE a weighted eulerian is constructed, otherwise weights are ignored.  A non-null \code{start} is the eulerian starting point.}

\item{self = "graphNEL"}{ Augments the graph using \code{mk_euler_graph}, then invokes eulerian again on the augmented verion. If self  is not connected, (approximate) eulerians are formed for each connected component, which are returned as a list. }

\item{self = "matrix"}{ Builds a graph using \code{mk_euler_graph}, then invokes eulerian again on the result.}

\item{self = "numeric"}{ Builds a graph with self nodes using \code{mk_euler_graph}, then invokes eulerian again on the result.}

\item{self = "ANY"}{ Builds a graph using \code{mk_euler_graph}, then invokes eulerian again on the result.}
}}

\references{C. Hierholzer (1873). Uber die Moglichkeit, einen 
	Linienzug ohne Wiederholung und ohne Unterbrechung zu umfahren. Math. Annalen VI, pp. 30-32.
	
	
Also, see \link{overview}}

\author{ C.B. Hurley and R.W. Oldford }

\keyword{methods}



\examples{

require(PairViz)

d <- as.matrix(eurodist)[1:8,1:8]   # pick the first 8 cities

eulerian(d)
eulerian(d, weighted=FALSE) # In this case, starts at city 1 and ends at city 8


}


