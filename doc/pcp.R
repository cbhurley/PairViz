## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(digits=3)

## ---- fig.align="center",fig.width=4, fig.height=2.5--------------------------
suppressPackageStartupMessages(library(PairViz))
data <- mtcars[,c(1,3:7)]
pcp(data,horizontal=TRUE,lwd=2, col="grey50",
main = "Standard PCP")

## ---- fig.align="center",fig.width=7, fig.height=2.5--------------------------
o <- hpaths(1:ncol(data),matrix=FALSE)
par(cex.axis=.7)
pcp(data,order=o,horizontal=TRUE,lwd=2, col="grey50",
main = "Hamiltonian decomposition")

## -----------------------------------------------------------------------------
corw <- as.dist(cor(data))

## -----------------------------------------------------------------------------
o <- eulerian(-corw)
o

## ---- fig.align="center",fig.width=7, fig.height=2.5--------------------------
par(cex.axis=.7)
pcp(data,order=o,horizontal=TRUE,lwd=2, col="grey50",
main = "Weighted eulerian")

## ---- fig.align="center",fig.width=7, fig.height=3----------------------------
corw <- dist2edge(corw)
edgew <- cbind(corw*(corw>0), corw*(corw<0))
par(cex.axis=.7)
guided_pcp(data,edgew, path=o,pcp.col="grey50" ,lwd=2,
         main="Weighted eulerian with correlation guide",
         bar.col = c("lightpink1","lightcyan1"),
         bar.ylim=c(-1,1),bar.axes=TRUE)


## ---- fig.align="center",fig.width=7, fig.height=2.5--------------------------
pathw <-  path_weights(corw,o)
corcols <- ifelse(pathw>0, "lightpink1", "lightcyan")
par(cex.axis=.7)
pcp(data,order=o,col="grey50" ,lwd=2,
         main="Weighted eulerian with correlation guide",
         panel.colors = corcols)

## -----------------------------------------------------------------------------
if (!requireNamespace("alr4", quietly = TRUE)){
    install.packages("alr4")
}

suppressPackageStartupMessages(library(alr4))
data(sleep1)
data <- na.omit(sleep1)
# these vars changed to factors in alr4, change from alr3
data$D <- as.numeric(data$D)
data$P <- as.numeric(data$P)
data$SE <- as.numeric(data$SE)

# logging the brain and body weights
data[,4:5] <- log(data[,4:5])

# short variable names
colnames(data) <- c("SW","PS" ,"TS" ,"Bd", "Br","L","GP","P" ,"SE" , "D"  )

# colours for cases, split Life values into 3 equal sized groups
cols1 <- scales::alpha(c("red","navy","lightblue3"   ),.6)
cols <- cols1[cut(rank(data$L),3,labels=FALSE)] 

## -----------------------------------------------------------------------------
library(scagnostics)
library(RColorBrewer)
sc <- scagnostics(data)
scags <- rownames(sc)
scags

## -----------------------------------------------------------------------------
scag_cols <- rev(brewer.pal(9, "Pastel1"))
names(scag_cols) <- scags

## -----------------------------------------------------------------------------
select_scagnostics <- function(sc,names){
	sc1 <- sc[names,]
	class(sc1) <- class(sc)
	return(sc1)
	}

## -----------------------------------------------------------------------------
scOut <- select_scagnostics(sc,"Outlying")
dOut <- edge2dist(scOut) # dOut is a dist
dOut <- as.matrix(dOut)
rownames(dOut) <- colnames(dOut)<- names(data)

## ----eval=F-------------------------------------------------------------------
#  o <- order_best(-dOut, maxexact=10)
#  o <-order_best(-dOut)
#  o <-order_tsp(-dOut)

## -----------------------------------------------------------------------------
o <-c( 2 , 4 , 1,  5 , 6,  7 , 3 , 8,  9, 10)

## ---- fig.align="center",fig.width=7, fig.height=3----------------------------
par(tcl = -.2, cex.axis=.8,mgp=c(3,.3,0))
guided_pcp(data,scOut, path=o,pcp.col=cols,lwd=1.4,
         main= "Best Hamiltonian for outliers",bar.col = scag_cols["Outlying"],legend=FALSE,bar.axes=TRUE,bar.ylim=c(0,max(scOut)))

## -----------------------------------------------------------------------------
outliers <- order(data$L, decreasing=T)[1:2]
rownames(data)[outliers]

## -----------------------------------------------------------------------------
colOut <- rep("grey50", nrow(data))
colOut[outliers[1]] <- "red" # Human
colOut[outliers[2]] <- "blue"

## -----------------------------------------------------------------------------
scSS <- t(select_scagnostics(sc,c("Striated", "Sparse")))

## ----eval=F-------------------------------------------------------------------
#  dSS <- edge2dist(scSS[,1])  + edge2dist(scSS[,2])
#  # You might think edge2dist(scSS[,1]+ scSS[,2]) would work, but as scSS[,1]+ scSS[,2] is
#  # not of class scagnostics, edge2dist will not fill the dist in the correct order
#  dSS <- as.matrix(dSS)
#  rownames(dSS) <- colnames(dSS)<- names(data)
#  order_best(-dSS,maxexact=10)

## ----eval=F-------------------------------------------------------------------
#  find_path(-scSS,   order_best,maxexact=10) # for the best path
#  # or
#  find_path(-scSS,   order_best) # for a nearly "best" path

## -----------------------------------------------------------------------------
o <-  c(4, 10 , 2 , 9 , 1,  7 , 8,  6 , 5 ,3) 

## ---- fig.align="center",fig.width=7, fig.height=3----------------------------
par(tcl = -.2, cex.axis=.8,mgp=c(3,.3,0))
guided_pcp(data,scSS, path=o,pcp.col=cols,lwd=1.4,
         main= "Best Hamiltonian for Striated + Sparse",
         bar.col = scag_cols[c("Striated", "Sparse")],
         legend=FALSE,bar.axes=TRUE,bar.ylim=c(0,.6))

## -----------------------------------------------------------------------------
o <- find_path(-scOut, eulerian)
o

## -----------------------------------------------------------------------------
pathw <-  path_weights(scOut,o)
head(pathw)

## ---- fig.align="center",fig.width=7, fig.height=3,fig.show='hold'------------
par(tcl = -.2, cex.axis=.6,mgp=c(3,.3,0))
guided_pcp(data,scOut, path=o[1:25],pcp.col=cols,lwd=1.4,
         main= "First 25: Eulerian for Outlying",
         bar.col = scag_cols["Outlying"],
         legend=FALSE,bar.axes=TRUE,bar.ylim=c(0,max(scOut)))

guided_pcp(data,scOut, path=o[25:50],pcp.col=cols,lwd=1.4,
         main= "Last 26: Eulerian for Outlying",
         bar.col = scag_cols["Outlying"],
         legend=FALSE,bar.axes=TRUE,bar.ylim=c(0,max(scOut)))

## -----------------------------------------------------------------------------
g <- mk_complete_graph(dOut)

## -----------------------------------------------------------------------------
g1 <- dn_graph(g,.2 , test=`>=`)

## -----------------------------------------------------------------------------
requireNamespace("igraph")
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

## ----fig.align="center", fig.width=3.5, fig.height=3.5------------------------
igplot(g1, weights=TRUE, layout=igraph::layout_as_tree)

## ----eval=FALSE---------------------------------------------------------------
#  e <- edgeMatrix(g1,duplicates=FALSE)
#  ew <- eWV(g1,e)
#  e <- matrix(nodes(g1)[e],ncol=2,byrow=TRUE)
#  g2 <- ftM2graphNEL(e,-ew,edgemode="undirected")

## -----------------------------------------------------------------------------
g2 <- dn_graph(mk_complete_graph(-dOut),-.2 , test=`<=`)

## ----fig.align="center",fig.width=5, fig.height=3-----------------------------
o <- eulerian(g2)
o
o <- match(o, names(data))
par(tcl = -.2, cex.axis=.86,mgp=c(3,.3,0))
guided_pcp(data,scOut, path=o,pcp.col=colOut,lwd=1.4,
         main= "Eulerian for high Outlying graph",
         bar.col = scag_cols["Outlying"],
         legend=FALSE,bar.axes=TRUE,bar.ylim=c(0,max(scOut)))

