\name{guided_pcp}
\alias{guided_pcp}

\title{ Guided parallel coordinate plot. }
\description{
 Draws a parallel coordinate plot, with an accompanying barchart showing an index (eg correlation, scagnostics) levels for each panel. An index legend is optional.
}
\usage{
guided_pcp(data, edgew=NULL, path = NULL, pathw=NULL,zoom=NULL,pcpfn=pcp,
     pcp.col = 1,lwd=0.5, panel.colors=NULL, pcp.mar=c(1.5,2,2,2), pcp.scale=TRUE,
     bar.col=1:9,bar.axes=FALSE, bar.mar=NULL,bar.ylim=NULL, reorder.weights=TRUE,
    layout.heights=NULL, layout.widths=c(10,1),
     main=NULL,legend=FALSE,cex.legend = 1,legend.mar=c(1,4,1,1),...)
     }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data frame or matrix. }
  \item{edgew}{Matrix (or vector) whose rows give index values for each pair of variables.}
  \item{path}{an index vector specifying variable order, or a function. If a function, \code{find_path(edgew,path,...)} constructs the index vector. }
   \item{pathw}{Matrix (or vector) whose rows give index values for each adjacent pair of variables in path. Usually this argument is NULL and \code{pathw} is computed from the \code{path} and \code{edgew}.}
\item{zoom}{If provided, a numeric vector specifying a subsequence of path to display.}
\item{pcpfn}{Function to draw the parallel coordinates.}
  \item{pcp.col}{ Line colors.  }
  \item{lwd}{ Line widths. }
  \item{panel.colors}{Background panel colors, passed to the\code{pcpfn} }
  \item{pcp.mar}{Controls PCP margin size. }
   \item{pcp.scale}{If TRUE, the variables will be scaled to 0-1 range, otherwise the data is not scaled. }
 \item{bar.col}{Bar colors. }
\item{bar.axes}{Draw barplot axes, if TRUE. }
  \item{bar.mar}{Controls barplot margin size. }
   \item{bar.ylim}{Vertical limits of bar plot. }
   \item{reorder.weights}{If TRUE, reorder barplot indices so large values are drawn at the bottom. }
\item{layout.heights}{Controls the layout.}
\item{layout.widths}{Controls the layout.}
\item{main}{Main title for PCP.}
\item{legend}{If TRUE, draws the barplot index legend.}
\item{cex.legend}{Controls legend text size.}
 \item{legend.mar}{Legend margin size. }
    \item{...}{ Optional arguments}
    }
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%  ~Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
%}
\references{see \link{overview}}

\author{ C.B. Hurley and R.W. Oldford }
\seealso{\code{\link{pcp}},\code{\link{catpcp}}}

\examples{
	 			
require(PairViz) 			


data <- mtcars[,c(1,3:6)]
cols <- c("red","green")[mtcars[,9]+1 ]    # transmission type, red=automatic

# add a correlation guide and find "better" hamiltonians...

# add a correlation guide...


corw <- dist2edge(as.dist(cor(data)))
edgew <- cbind(corw*(corw>0), corw*(corw<0))
         
# add a correlation guide to a PCP, positive cors shown in blue, negative in purple...

\dontrun{
dev.new(width=3,height=3)

par(cex.axis=.65)

guided_pcp(data,edgew, pcp.col=cols,
         main="Correlation guided PCP",bar.col = c("blue","purple"))

dev.new(width=7,height=3)
par(cex.axis=.65)

guided_pcp(data,edgew, path=eulerian, pcp.col=cols,lwd=2,
         main="Correlation guided Eulerian PCP",bar.col = c("blue","purple"),bar.axes=TRUE)



}

# Scagnostic guides are useful here- see the demos for more examples.
}

\keyword{hplot}
