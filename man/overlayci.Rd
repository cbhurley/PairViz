\name{overlayCI}
\alias{overlayCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to overlay confidence intervals on the current plot. }
\description{
  Overlays confidence intervals on the current plot. Also draws a right hand axis, a horizontal broken line at zero, and marks the significant comparisons with an arrow, i.e. the CIs that do not intersect zero.
}
\usage{
overlayCI(cis, xpos=NULL,ci.cols = NULL, ci.ex = 2, ci.ocol = "grey40", 
p.col = "grey40", pch = 1, sig.col = "red", sig.lwd = 1, yusr = NULL,
ci.label="Differences",ci.cex=0.5,arrow.length=0.1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cis}{A matrix containing the confidence intervals. Each row corresponds to a
  	different comparison, the first column is the estimated mean, and succesive pairs of columns give the lower and upper limits for different confidence levels. }
  \item{ci.cols}{ A vector of colours, one colour for each confidence level. Defaults to shades of grey. }
  \item{ci.ex}{Controls confidence interval line width. }
  \item{xpos}{Horizonal positions where CIs are drawn. Defaults to 1.5,2.5,3.5,..}
  \item{ci.ocol}{Colour of zero line.}
  \item{p.col}{Colour of point used for CI centre. }
  \item{pch}{ Symbol used for CI centre. }
  \item{sig.col}{ Colour of arrow marking significant comparisons. }
  \item{sig.lwd}{Width of arrow marking significant comparisons. }
  \item{yusr}{Specifies the vertical \code{par(usr) }.Defaults to max and min. }
  \item{ci.label}{Label drawn on right margin. }
  \item{ci.cex}{Controls size of CI mean point symbol. }
  \item{arrow.length}{Controls size arrow at right hand axis. }
  \item{...}{Ignored}
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
\note{ This function is called by \code{mc_plot}}
\seealso{ See Also as \code{\link{mc_plot}} }

\keyword{aplot}

