\name{mc_plot}
\alias{mc_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multiple comparison plot.  }
\description{
 For grouped data. Draws boxplots for each group and overlays with confidence intervals for pairwise comparison of means.
}
\usage{
mc_plot(data, fit, path = eulerian, col = rainbow(length(data), s = 0.4), 
levels = c(0.9, 0.95, 0.99), varwidth = TRUE, frame.plot = FALSE, 
boxwex = 0.3, cex=0.75, zoom=NULL,ci.yusr=NULL,ci.pos=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A list of vectors, such as that returned by \code{split}. }
  \item{fit}{ Either an \code{\link{aov}} fit, or else  a matrix with columns estimate, followed by confidence intervals. If fit is not an \code{\link{aov}} fit, the path argument should be a vector.}
  \item{path}{an index vector or a function. If a vector, groups are plotted in order \code{data[path]}. By default, it is the  function \code{\link{eulerian}}, and produces an ordering where each pair of groups appears adjacently, with p-values roughly increasing as the sequence progresses.  }
  \item{col}{ A vector of colours, one per group.  }
  \item{levels}{ Vector of increasing confidence levels. }
  \item{varwidth}{Passed to  \code{\link{boxplot}}.  }
  \item{frame.plot}{ Passed to  \code{\link{boxplot}}.}
  \item{boxwex}{ Passed to  \code{\link{boxplot}}. }
   \item{cex}{ Passed to  \code{\link{boxplot}}. }
  \item{zoom}{If provided, a numeric vector specifying a subsequence of path to display.}
  \item{ci.yusr}{Specifies the vertical \code{par(usr) } for the confidence intervals. Defaults to max and min. }
   \item{ci.pos}{If TRUE, all CIs are mu(max) - mu(min), otherwise mu(right) - mu(left). }

  \item{...}{ Optional arguments, passed to  \code{boxplot} and \code{\link{overlayCI}}.}
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

\seealso{ See also  \code{\link{overlayCI}} }
\examples{
	 			
require(PairViz) 			

data(cancer)
bx <- with(cancer, split(sqrt(Survival),Organ))



a <-  aov(sqrt(Survival) ~ Organ,data=cancer)
\dontrun{
dev.new(height=4.5, width=9.5)
op <- par(no.readonly = TRUE)

par(cex.axis=.75, cex.main = 1.0, cex.lab=1)
par(mar=c(3,5,3,5))

mc_plot(bx,a,main="Pairwise comparisons of cancer types", ylab="Sqrt Survival")

par(op)
}
}

\keyword{hplot}
