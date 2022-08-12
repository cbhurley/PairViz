\name{knn_graph}
\Rdversion{1.1}
\alias{knn_graph}
\alias{dn_graph}
\alias{mk_hypercube_graph}
\alias{mk_binary_graph}

\alias{mk_line_graph}
\alias{kspace_graph}

\alias{graph_product}
\alias{graph_compose}
\alias{graph_sum}

\alias{bipartite_graph}
\alias{iterated_line_graph}




%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to construct graphs.}
\description{
Functions to construct graphs- see details below.}
\usage{
knn_graph(g, k = 2)
dn_graph(g, d = 1, test=`<=`)
mk_binary_graph(n,sep="",delta=1,test=`==`)
mk_hypercube_graph(n,sep="")
mk_line_graph(g,sep="-")
kspace_graph(n,m, link=NULL,sep="-")
graph_product(g,h, type="cartesian",sep="-")
graph_compose(g,h,sep="-")
graph_sum(g,h,combineWeight=`+`)

bipartite_graph(n1,n2)
iterated_line_graph(g,sep="-")
}
\arguments{
  \item{g}{a graph}
  \item{h}{a graph}
  \item{n}{a positive integer, or a character vector.}
  \item{k}{a positive integer}
  \item{d}{an edge weight}
  \item{test}{used to select edges.}
  \item{sep}{used to form node names of new graph.}
  \item{m}{subsets of size m are nodes of kneser graph. }
  \item{link}{A positive number or NULL. If NULL, the returned graph is complete. Otherwise edges for subsets sharing link elements. }
  \item{type}{the type of graph product, one of "cartesian", "strong" or "tensor" }
    \item{n1}{a character vector.}
    \item{n2}{a character vector.}
\item{delta}{used to select edges.}
\item{combineWeight}{used to combine weights.}
  }
\details{
knn_graph- returns a symmetric k nearest neighbour graph

dn_graph- returns a  graph formed from g where edges of satisfy test(weight, d). The default retains edges whose weight are 1 are less. Nodes with no edges are also removed.

mk_hypercube_graph- returns a hypercube graph with $2^n$ nodes

mk_binary_graph(n,sep="",delta=1,test=`==`) - returns a  graph with \eqn{2^n} nodes. Undirected edges join nodes A and B whose binary vectors satisfy \eqn{a_i \leq b_i,i=1, \ldots, n} and test(\eqn{\sum(b_i - a_i)}, delta) is true.

mk_line_graph- returns the line graph of g

kspace_graph- returns a  graph where nodes  are subsets of size m from n. Edges are connect nodes whose subsets share link elements. The standard kneser graph has link=0. When link is NULL, returned graph is complete.

graph_product(g,h, type="cartesian",sep="-")- returns the graph product of g and h.

graph_compose(g,h,sep="-")- returns the graph composition of g and h.

bipartite_graph(n1,n2)- returns the complete bipartite graph with node sets n1 and n2.


graph_sum(g,h,combineWeight=`+`)- returns a graph whose nodes and edges are the union of those in g and h. Weights of common edges are combined using combineWeight.

iterated_line_graph- returns the iterated line graph of g, with compression of nodes as described in the reference Hurley and Oldford(2008) given below.


}

\references{
See any Graph Theory text. Also
C.B. Hurley and R.W. Oldford,
Graphs as navigational infrastructure for high dimensional data spaces. 2008 submitted.
}
\author{ C.B. Hurley and R.W. Oldford }

\examples{
	
# See the demo file nav.R	
}