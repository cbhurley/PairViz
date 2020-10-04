\name{table_plot}
\alias{table_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots rectangles on a grid }
\description{
   Plots rectangles on a grid-  a barchart/mosaic variant which facilitates pairwise comparisons.}
\usage{
table_plot(rectw, recth, col="grey50", gapx = NULL, gapy = NULL, 
spacex = 0.03, spacey = 0.03, xjust = "center", yjust = "center", 
xruler = NULL, yruler = NULL, color.ruler = "grey32",
pch0=1,xlab=NULL,ylab=NULL, plot=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rectw}{ An n*m matrix of  rectangle widths, or a vector of m column widths. }
  \item{recth}{ An n*m  matrix of  rectangle heights, or a vector of n row heights.}
  \item{col}{ Rectangle fill colours. }
  \item{gapx}{ Gaps in the x-direction. If provided should be a vector of length m-1. }
  \item{gapy}{ Gaps in the x-direction. If provided should be a vector of length n-1. }
  \item{spacex}{ A single value- extra space between columns as a fraction of maximum row total of rectw .}
  \item{spacey}{ A single value- extra space between rows as a fraction of maximum column total of recth . }
  \item{xjust}{Horizontal justification of rectangles- "center", "left", or "right".}
  \item{yjust}{ Vertical justification of rectangles- "center", "bottom", or "top".}
  \item{xruler}{Specifies position of rulers drawn parallel to x-axis. Values are a subset of ("top","center","bottom")}
  \item{yruler}{ Specifies position of rulers drawn parallel to y-axis. Values are a subset of ("left","center","right") }
  \item{color.ruler}{ Color for the rulers. }
  \item{pch0}{ Symbol for zero cell size. May be NULL. }
   \item{xlab}{ X label }
 \item{ylab}{ Y label }
  \item{plot}{ If TRUE, draw tge plot. Otherwise returns a matrix where each row is the coordinates of a the calculated rectangle. }

  \item{\dots}{ Passed to plot. }
}
\references{See \link{overview}}
\author{ Catherine Hurley}

\seealso{ See also \code{\link{barplot}}, \code{\link{mosaicplot}} }
\examples{
\dontrun{
require(PairViz) 	
		
tab <- apply(HairEyeColor, c(1, 2), sum)


dev.new()
par(mar=c(3,3,1,1))
par(cex=.6,mgp=c(2, -.5, 0))
table_plot(sqrt(tab),sqrt(tab))
# this table plot has cells with widths and heights proportional to the square root of cell counts.

tabp <- prop.table(tab,2)

table_plot(apply(tab,2,sum),tabp) # make cell widths proportional to 
#margin totals, heights to conditional prob

cols <- 2:5
table_plot(apply(tab,2,sum),tabp, yjust="bottom",col=cols,yruler=c("left","right")) 
# add colours, rulers and  bottom-justify

# The result is similar to the mosaic, without the mosaic effect of equalizing gaps. 
#In the table version the rectangles line up across rows, 
#so comparing heights, ie. conditional probs is easier.

o <- hpaths(1:4)[2,]
table_plot(apply(tab,2,sum)[o],tabp[,o], yjust="bottom",col=cols,yruler=c("left","right"))
# Permutes the columns so all pairs of columns can be compared. 
#In the second permutation can easily see that   
#p(black|blue eyes)> p(black|green eyes)


dev.new()
par(mar=c(3,3,1,1))
par(mgp=c(2, -.5, 0))
mosaicplot(t(tab)[,nrow(tab):1],col=rev(cols),main="")
# mosaic- good for seeing deviations from independence. hard to compare conditional probs,
# except for those in the bottom and top rows. 
}
	
}
\keyword{hplot}
