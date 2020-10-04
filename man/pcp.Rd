\name{pcp}
\alias{pcp}
\alias{catpcp}
\title{Enhanced parallel coordinate plots. }

\description{
pcp draws a parallel coordinate plot.It is a modified
version of \code{parcoord {MASS}}.
 Variables
may be reordered and panels colored in the display. 

catpcp draws a parallel coordinate plot variant for categorical data.}
\usage{
pcp(data, order = NULL, panel.colors = NULL, col = 1, lty = 1, 
horizontal = TRUE, mar = NULL, scale = TRUE, axis.width = 0,
axis.grid.col="grey70",connect=TRUE, ...)

catpcp(data, order = NULL, pcpbars, barvars = 1:ncol(data), 
pcpbars.border = "black", pcpbars.col = NULL, pcpbars.labels = FALSE, 
pcpbars.axis.at = NULL, pcpbars.axis.labels = NULL, 
axis.width = 0.2,connect=TRUE, ...)


}
\arguments{
\item{data}{ A data frame or matrix. }
\item{order}{an index vector specifying variable order. If NULL, all variables are used.
}
  \item{panel.colors}{either a vector or a matrix of panel colors. If a 
  vector is supplied, the ith color is used for the ith panel. If a matrix, dimensions 
  should match those of the variables. Diagonal entries are ignored. }
 \item{col}{ a vector of colours, recycled as necessary for each observation.  }
  \item{lty}{ a vector of line types, 
recycled as necessary for each observation.  }
  \item{horizontal}{ If TRUE, orientation is horizontal. }
  \item{mar}{ margin parameters, passed to \code{par}. }
   \item{scale}{
If TRUE, the variables are scaled to the unit interval.}
  \item{axis.width}{Width of each of the parallel axes.
}
 \item{axis.grid.col}{Color of variable axes. Use NULL for no axes.}
\item{connect}{If FALSE, line segments are not connected. Points are drawn if axis.width=0.
}
  \item{pcpbars}{A list, with one component per barvar. Component i is a matrix with the bottom and top of the bars for that variable.}
  \item{barvars}{Categorical variables where overlayed bars show the level frequency.}
  \item{pcpbars.border}{Border colour of the bars.}
  \item{pcpbars.col}{Colour of the bars.
}
  \item{pcpbars.labels}{Labels for the bars.
}
  \item{pcpbars.axis.at}{Axis label positions for the bars.
}
  \item{pcpbars.axis.labels}{Axis label text for the bars.}
   \item{\dots}{ other parameters, passed to pcp by catpcp}

}
\examples{
require(PairViz) 		
y <- as.data.frame(as.table(HairEyeColor))

colvar <- 3 # any of 1:3 will do
y <- y[order(y[,colvar]),] # ensures that cases are ordered by colour within each factor level
ylong <- apply(y[,-4],2, function(x) rep(x,times=y[,4]))

cols <- desaturate_color(rainbow(4,alpha=0.3),.5)
cols <- cols[as.numeric(as.factor(ylong[,colvar]))]




ds <- factor_spreadout(ylong)




dev.new(width=5,height=2.5)
par(mar=c(2,1,2,1))
par( cex.axis=.8,cex=.8)



catpcp(ds$data,col=cols,lwd=2,pcpbars=ds$bars,pcpbars.labels=TRUE,main="Hair Eye data")

	
	
	
	
	}
