\name{spreadout}
\alias{rater_spreadout}
\alias{factor_spreadout}

\title{
Functions to prepare for categorical parallel coordinates, drawn by catpcp.
}
\description{
factor_spreadout spreads out the data at each factor level.
rater_spreadout spreads out the data at each rating level.
The rater version is appropriate when the variables (factors) have all the same levels.
}
\usage{
factor_spreadout(d)
rater_spreadout(d, levs, minspace = NULL,scale=FALSE)
}

\arguments{
  \item{d}{A data frame where each variable can be interpreted as a factor.
}
  \item{levs}{The rating levels. Specifying this controls the order of rating levels on each axis.}
  \item{minspace}{The minimum amount of space between the bars.}
  \item{scale}{If scale=FALSE, the ith rater values are spreadout about the value i. If scale=TRUE, all values are scaled to 0-1.}

}

\details{
factor_spreadout spreads out the data at each factor level. It returns a list with two components. The first is data, containing the spreadout data, scaled to 0-1. The second is bars, which is a list whose ith component gives the bottom and top of the bars for the ith variable of d.

rater_spreadout spreads out the data at each rater level. It returns a list with two components. The first is data, containing the spreadout data. If scale=FALSE, the ith rater values are spreadout about the value i. If scale=TRUE, all values are scaled to 0-1. The second component is bars, which is a list whose ith component gives the bottom and top of the bars for the ith variable of d.
}