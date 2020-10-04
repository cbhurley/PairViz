\name{mk_complete_graph}\alias{mk_complete_graph}  \title{Constructs a complete graph.}\description{ Constructs a complete graph, actually an instance of \code{graph-NEL}}\usage{mk_complete_graph(d)}\arguments{  \item{d}{an integer vector of length 1 which specified the number of nodes, a character vector of nodes names,  a \code{dist}, or a symmetric matrix, either of which specify the nodes and edge weights.}}\value{- a \code{graph-NEL}}\author{ C.B. Hurley and R.W. Oldford }\keyword{graphs}\examples{
	
require(PairViz)
d <- dist(rnorm(5))
g <- mk_complete_graph(d)}