## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("PairViz")

## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("graph", quietly = TRUE)){
#      install.packages("BiocManager")
#      BiocManager::install("graph")
#  }

## -----------------------------------------------------------------------------
requireNamespace("igraph")
if (!requireNamespace("igraph", quietly = TRUE)){
    install.packages("igraph")
}


igplot <- function(g,weights=FALSE,layout=igraph::layout_in_circle, 
                   vertex.size=60, vertex.color="lightblue",...){
    g <- igraph::graph_from_graphnel(as(g, "graphNEL"))
    op <- par(mar=c(1,1,1,1))
    if (weights){
      ew <- round(igraph::get.edge.attribute(g,"weight"),3)  
      igraph::plot.igraph(g,layout=layout,edge.label=ew,vertex.size=vertex.size,vertex.color=vertex.color,...)
    }
    else
    igraph::plot.igraph(g,layout=layout,vertex.size=vertex.size,vertex.color=vertex.color,...)
    par(op)
}

## ----eval=FALSE---------------------------------------------------------------
#  
#  if (!requireNamespace("Rgraphviz", quietly = TRUE)){
#      install.packages("BiocManager")
#      BiocManager::install("Rgraphviz")
#  }
#  
#  # For a graph g use
#  plot(g)

## ----fig.show='hold'----------------------------------------------------------
suppressPackageStartupMessages(library(PairViz))
k4 <- mk_complete_graph(4)
k5 <- mk_complete_graph(5)
igplot(k4)
igplot(k5)

## -----------------------------------------------------------------------------
eulerian(k5)

## ----fig.show='hold', fig.align='center'--------------------------------------
k5a <- removeEdge("1", "2", k5)
igplot(k5a)
eulerian(k5a)

## -----------------------------------------------------------------------------
eulerian(k4)

## -----------------------------------------------------------------------------
n <- LETTERS[1:5]
g <- new("graphNEL",nodes=n)
efrom <- n[c(1,1,2,2,2,4)]
eto <- n[c(2:3,3:5,5)]
ew <- c(8,9,5:7,1)
g <- addEdge(efrom, eto, g, ew) 

## ----fig.align="center"-------------------------------------------------------
igplot(g, weights=TRUE,edge.label.color="black")

## -----------------------------------------------------------------------------
eulerian(g)

## -----------------------------------------------------------------------------
eulerian(g, weighted=FALSE)

## ----fig.show='hold'----------------------------------------------------------
ew <- rep(1, length(edgeNames(k5)))
s1 <- c(1,5,9,10)
ew[s1]<- 5
ec <- rep("grey40", length(ew))
ec[s1]<- "cyan"
igplot(k5, edge.width=ew, edge.color=ec)

ew[]<- 1
s2 <- c(2,6:8)
ew[s2]<- 5
ec[] <- "grey40"
ec[s2]<- "magenta"
igplot(k5, edge.width=ew, edge.color=ec)

## ----fig.align='center'-------------------------------------------------------
ew[]<- 5
ec[s1]<- "cyan"
s3 <- 3:4
ec[s3]<- "rosybrown1"
igplot(k5, edge.width=ew, edge.color=ec)
igplot(k5, edge.width=ew, edge.color=ec)

## -----------------------------------------------------------------------------
hpaths(5)

## -----------------------------------------------------------------------------
hpaths(5, matrix=FALSE)

## -----------------------------------------------------------------------------
hpaths(6)
hpaths(6, matrix=FALSE)

## -----------------------------------------------------------------------------
hpaths(1:5)

## -----------------------------------------------------------------------------
set.seed(123)
k7 <- mk_complete_graph(7)
ew <- sample(numEdges(k7),numEdges(k7)) # a vector of edgeweights

## -----------------------------------------------------------------------------
d7 <- matrix(0,7,7)
d7[lower.tri(d7)] <- ew
d7[upper.tri(d7)]<-  t(d7)[upper.tri(d7)]
d7 
# or using the shortcut function edge2dist from PairViz
#d7 <- as.matrix(edge2dist(ew))

## ----fig.align="center",fig.width=6, fig.height=6-----------------------------
k7 <- mk_complete_graph(d7)
igplot(k7, weights=TRUE,edge.label.color="black", vertex.label.cex=2,vertex.size=30)

# Unfortunately, plot.igraph does not show graph edge weights  automatically, you have to 
# input them as above. You might want to check that the igraph
# matches that of ew.
igraph::E(igraph::graph_from_graphnel(k7))


## -----------------------------------------------------------------------------
weighted_hpaths(d7)
# this version returns the eulerian
weighted_hpaths(d7, matrix=FALSE)

## -----------------------------------------------------------------------------
o <- weighted_hpaths(d7, matrix=FALSE)
o1 <- o[1:8] # include the 8th to form the tour
d7e <- dist2edge(d7) 
# d7e is a vector giving edge weights in order (1,1)... (1,7), (2,3),.. (2,7) etc
h1weights <- path_weights(d7e, o1) # the edge weights for o1
# the same as
d7[cbind(head(o1,-1), o1[-1])]

h1weights
sum(h1weights)

## -----------------------------------------------------------------------------
o2 <- o[8:15]
sum(path_weights(d7e, o2))
o3 <- o[15:22]
sum(path_weights(d7e, o3))

## -----------------------------------------------------------------------------
order_best(d7,cycle=TRUE)
order_tsp(d7,cycle=TRUE)

## -----------------------------------------------------------------------------
e1 <- eseq(7)
e2 <- eseqa(7)
e3 <- eulerian(7) # same path as eulerian(k7, weighted=FALSE)
h1 <- hpaths(7, matrix=FALSE)

## ----fig.width=7, fig.height=6, echo=FALSE------------------------------------
par(mfrow=c(2,2))
par(mar=c(2,2,3,1))

plot(e1, type="n",main="eseq(7)")
lines(e1[1:4], col="tan1", lwd=1.5)
lines(4:11,e1[4:11], col="cyan", lwd=1.5)
lines(11:22,e1[11:22], col="magenta", lwd=1.5)

points(e1, pch=20); grid()

plot(e2, type="n",main="eseqa(7)")
e2x <- rep(NA, length(e2))
gap1 <- seq(1, length(e2), by=3)
gap1 <- sort(c(gap1, gap1+1))

e2x[gap1]<- e2[gap1]
lines(e2x, col="tan1", lwd=1.5)

e2x <- rep(NA, length(e2))
gap2 <- gap1+1
e2x[gap2]<- e2[gap2]
lines(e2x, col="cyan", lwd=1.5)

e3x <- rep(NA, length(e2))
gap3 <- gap2+1
e3x[gap3]<- e2[gap3]
lines(e3x, col="magenta", lwd=1.5)
points(e2, pch=20);grid()


plot(e3, type="n", main="eulerian(7)")


e3x <- rep(NA, length(e3))
#indy <- which(head(e3,-1) <= 1 | e3[-1] <= 1)
indy <- which(pmin(head(e3,-1), e3[-1]) <= 1)
e3x[indy] <- e3[indy]
e3x[indy+1] <- e3[indy+1]
lines(e3x, col="tan1", lwd=1.5)


e3x <- rep(NA, length(e3))
#indy <- which(head(e3,-1) %in% c(2,3) | e3[-1] %in% c(2,3))
indy <- which(pmin(head(e3,-1), e3[-1]) %in% c(2,3))

e3x[indy] <- e3[indy]
e3x[indy+1] <- e3[indy+1]
lines(e3x, col="cyan", lwd=1.5)



e3x <- rep(NA, length(e3))
#indy <- which(head(e3,-1) >=4 & e3[-1] >=4)
indy <- which(pmin(head(e3,-1), e3[-1]) >=4)
e3x[indy] <- e3[indy]
e3x[indy+1] <- e3[indy+1]
lines(e3x, col="magenta", lwd=1.5)

points(e3, pch=20);grid()

plot(h1, type="n", main="hpaths(7)"); 
lines(1:8,h1[1:8], col="cyan", lwd=1.5); lines(8:15,h1[8:15], col="magenta", lwd=1.5)
lines(15:22,h1[15:22], col="tan1", lwd=1.5)
points(h1, pch=20);grid()

## -----------------------------------------------------------------------------
e4 <- eulerian(d7) # same path as eulerian(k7)
h2 <- weighted_hpaths(d7, matrix=FALSE)

## ----fig.width=7, fig.height=3, echo=FALSE------------------------------------
par(mfrow=c(1,2))
par(mar=c(2,2,3,1))



plot(e4, type="n", main="eulerian(d7)")

e4x <- rep(NA, length(e4))
indy <- which(pmin(head(e4,-1), e4[-1]) %in% c(2,3))
e4x[indy] <- e4[indy]
e4x[indy+1] <- e4[indy+1]
lines(e4x, col="cyan", lwd=1.5)

e4x <- rep(NA, length(e4))
indy <- which(pmin(head(e4,-1), e4[-1]) <= 1)

e4x[indy] <- e4[indy]
e4x[indy+1] <- e4[indy+1]
lines(e4x, col="tan1", lwd=1.5)

e4x <- rep(NA, length(e4))
indy <-  which(pmin(head(e4,-1), e4[-1]) >=4)

e4x[indy] <- e4[indy]
e4x[indy+1] <- e4[indy+1]
lines(e4x, col="magenta", lwd=1.5)



points(e4, pch=20);grid()

plot(h2, type="n", main="weighted_hpaths(d7)"); 
lines(1:8,h2[1:8], col="cyan",lwd=1.5); lines(8:15,h2[8:15], lwd=1.5,col="magenta")
lines(15:22,h2[15:22], col="tan1", lwd=1.5)
points(h2, pch=20);grid()

## -----------------------------------------------------------------------------
d7e <- dist2edge(d7) 
path_weights(d7e, e4) # the edge weights for e4
path_weights(d7e, h2) # the edge weights for h2

## ----fig.width=7, fig.height=3, echo=FALSE------------------------------------
par(mfrow=c(1,2))
par(mar=c(2,2,3,1))

e4w <- path_weights(d7e, e4)
plot(e4w, type="l", col="grey80", main="edge weights of eulerian(d7)", ylab="weight")
points(e4w, pch=21, bg="cyan");grid()

h2w <- path_weights(d7e, h2)
plot(h2w, type="l", col="grey80", main="edge weights of hpaths(d7)", ylab="weight")

points(1:7,h2w[1:7], pch=21, bg="cyan", col="grey50")
points(8:14,h2w[8:14], pch=21, bg="magenta", col="grey50")
points(15:21,h2w[15:21], pch=21, bg="tan1",col="grey50")
grid()


