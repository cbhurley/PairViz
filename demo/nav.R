# This file constructs some of the graphs used in "Graphs as navigational infrastructure for high dimensional data spaces."

library(PairViz)

# one of igraph (CRAN) or Rgraphviz (bioconductor) is required to draw the graphs

# install igraph with 
#install.packages("igraph")


# install Rgraphviz with the two lines below
#source("https://bioconductor.org/biocLite.R")
# biocLite("Rgraphviz")


# In this demo, we will also show plots of the graph structure, which requires package
# igraph from CRAN. The following utility function will be helpful:

requireNamespace("igraph")
igplot <- function(g,...){
    g <- igraph::graph_from_graphnel(as(g, "graphNEL"))
    op <- par(mar=c(1,1,1,1))
    igraph::plot.igraph(g,layout=igraph::layout.circle,vertex.size=60,vertex.color="lightblue",...)
    par(op)
  }
 



#----------------
# Versions of Fig 1,2,3
n <- colnames(iris)[1:4] <- c("SL","SW","PL","PW")
g <- mk_complete_graph(n)
igplot(g)

space2 <- kspace_graph(n,2)
igplot(space2)

lg <- mk_line_graph(g)
igplot(lg)



lgc <- complement(lg)
igplot(lgc)

space21 <- kspace_graph(n,2,1) # same as lg

space20 <- kspace_graph(n,2,0) # same as lgc

space32 <- kspace_graph(n,3,2)
igplot(space32)



#----------------
# Versions of Fig 7,8

n <- LETTERS[1:5]
g <- mk_complete_graph(n)
lg <- mk_line_graph(g,sep="")
igplot(lg)

petersen <- complement(lg)
igplot(petersen)

space31 <- kspace_graph(n,2,1) # same as lg
space30 <-  kspace_graph(n,2,0) # same as petersen

igplot(space32)

space31 <- kspace_graph(n,3,1)
igplot(space31)



#----------------
# Versions of Fig 9,10
require(gclus)
data(ozone)
n <- colnames(ozone)

g <- new("graphNEL", nodes=n)
igplot(g)
g <- addEdge(n[1], n[2:5],g)
igplot(g)
 
lg <- mk_line_graph(g)
igplot(lg)

    
d <- abs(cor(ozone, method="spearman"))
g <- mk_complete_graph(d)
q <- quantile(as.dist(d),.9)
g1 <- dn_graph(g,q, `>=`)
igplot(g1)
lg1 <- mk_line_graph(g1)
igplot(lg1)
igplot(complement(lg1))


#----------------
# Versions of Fig 14

n <- LETTERS[1:5]
space32 <-  kspace_graph(n,3,2) 


space31 <- kspace_graph(n,3,1) 

igplot(space31)

space32c <- complement(space32) # same as space31


#----------------
# Versions of Fig 15-19

xn <-paste("X",1:3,sep="")
yn <- paste("Y",1:2,sep="")


g <- PairViz::bipartite_graph(xn,yn) # same function name appears in igraph
igplot(g)
lg <- mk_line_graph(g,sep="")
igplot(lg)
igplot(complement(lg))

space3d <- mk_line_graph(lg)
igplot(space3d)

space3d <- iterated_line_graph(g,sep=".")
igplot(space3d) # this version does compression, not needed here


u <-  new("graphNEL", nodes=paste("U",1:3,sep=""))
u <- addEdge("U2",c("U1","U3"),u)
igplot(u)
v <-  new("graphNEL", nodes=paste("V",1:2,sep=""))
v <- addEdge("V2","V1",v)
igplot(v)

g1 <- graph_product(u,v,"cartesian")
igplot(g1)


g2 <- graph_product(u,v,"tensor")
igplot(g2)

g3<- graph_product(u,v,"strong")
igplot(g3)


 
g4<- graph_compose(u,v)
igplot(g4)

g5<- graph_compose(v,u)
igplot(g5)