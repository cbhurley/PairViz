\name{etour}
\alias{etour}  

\title{Constructs eulerian tours on a graph. }
\description{
\code{etour}-- Constructs an eulerian tour on a graph using Hierholzer's algorithm. Returns a vector of node labels. If \code{weighted} is \code{TRUE} constructs a weight-decreasing eulerian using the modified Hierholzer's algorithm. Usually \code{etour} is not called directly, rather the generic function \code{eulerian} is used. 


}
\usage{
etour(g, start=NULL,weighted=TRUE)

}
\arguments{
  \item{g}{a graph satisfying \code{is_even_graph}}
  \item{start}{an optional starting node for the tour.}
  \item{weighted}{whether tour uses weights}
  
  	 }
     



\details{The supplied graph should satisfy\code{is_even_graph}. If \code{weighted} is \code{TRUE}
the lowest weight edge is found, and the tour starts at the
one of its nodes, picking the node with the bigger second-smallest edge weight.
 After that the tour follows weight-increasing edges.
If \code{weighted} is \code{FALSE}
weights are ignored. The returned tour is typically a closed path. However, if the last edge is a duplicated edge added to make the graph even, this edge is omitted and the result is an open path.	

	}
	
\references{see \link{overview}}
\author{ C.B. Hurley and R.W. Oldford }
\keyword{graphs}
\examples{
require(PairViz)

g <- mk_even_graph(5)

etour(g)
g <- mk_even_graph(6) # adds 3 extra edges to g, so all nodes are even
etour(g)
etour(g, start= "4") # modifies the starting node

eulerian(6) # The eulerian wrapper looks after making even graph,
#also returns numbers rather than nodes

# On a general graph.
v <- LETTERS[1:4]
g <- new("graphNEL",nodes=v)
g <- addEdge(v[1],v[3:4],g,1:2)
g <- addEdge(v[2],v[3:4],g,3:4)
etour(g) 

eulerian(g) # Equivalently, use eulerian wrapper 


n <- LETTERS[1:5]
g <- new("graphNEL",nodes=n)
g <- addEdge(n[1],n[2:3],g)
g <-addEdge(n[2],n[3:5],g)
g <-addEdge(n[4],n[3],g)
is_even_graph(g)
etour(mk_even_graph(g))

eulerian(g) # Equivalently, use eulerian wrapper 


}
