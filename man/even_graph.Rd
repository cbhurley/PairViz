\name{even_graph}
\docType{class}
\alias{even_graph-class}
\alias{graph_NEL-class}
\alias{is_even_graph}
\alias{is_even_graph,even_graph-method}
\alias{is_even_graph,graphNEL-method}

\title{Class of graphs where all nodes have even degree}
\description{This class is an extension of \code{graphNEL-class}. For graphs of this class, euler tours may always be constructed. Objects of this class should  be created by  \code{mk_even_graph}}


\section{Slots}{
	This class has all slots from  \code{\link[graph]{graphNEL-class}} plus:
	 \describe{
    \item{\code{dummy_node}:}{Object of class \code{"character"}  }
    \item{\code{extra_edges}:}{Object of class \code{"character"}  }
    \item{\code{weighted}:}{Object of class \code{"logical"}  }
   % \item{\code{nodes}:}{Object of class \code{"vector"} ~~ }
   %  \item{\code{edgeL}:}{Object of class \code{"list"} ~~ }
   %  \item{\code{edgemode}:}{Object of class \code{"character"} ~~ }
   %  \item{\code{edgeData}:}{Object of class \code{"attrData"} ~~ }
   %  \item{\code{nodeData}:}{Object of class \code{"attrData"} ~~ }
   %  \item{\code{renderInfo}:}{Object of class \code{"renderInfo"} ~~ }
   %  \item{\code{graphData}:}{Object of class \code{"list"} ~~ }
  }
}

\section{Extends}{
Class \code{\link[graph]{graphNEL-class}}, directly.
Class \code{\link[graph]{graph-class}}, by class "graphNEL", distance 2.
}


\section{Methods}{
  \describe{
\item{is_even_graph}{\code{signature(g = "graphNEL")}: checks whether a graph has all nodes of even degree.}
      
\item{is_even_graph}{\code{signature(g = "even_graph")}: always TRUE.}
}}

\author{ C.B. Hurley and R.W. Oldford }

\references{see \link{overview}}



%\seealso{}
\examples{
showClass("even_graph")

}
\keyword{classes}

\keyword{graphs}
