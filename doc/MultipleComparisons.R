## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(digits=4)

## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(PairViz))
data(cancer)  # Need this step to load the data
str(cancer)   # Summary of structure of the data
# We can separate the survival times by which organ is affected
organs <- with(cancer, split(Survival, Organ))
# And record their names for use later!
organNames <- names(organs)
# the structure of the organs data
str(organs)

## ---- fig.align="center", fig.width=6, fig.height=5---------------------------

library(colorspace)
cols <- rainbow_hcl(5, c = 50)  # choose chromaticity of 50 to dull colours
boxplot(organs, col=cols, 
        ylab="Survival time", 
        main="Cancer treated by vitamin C")

## ---- fig.align="center", fig.width=6, fig.height=5---------------------------
# Split the data
sqrtOrgans <- with(cancer, split(sqrt(Survival), Organ))
boxplot(sqrtOrgans, col=cols, 
        ylab=expression(sqrt("Survival time")), 
        main="Cancer treated by vitamin C")

## -----------------------------------------------------------------------------
ord <- eulerian(5)
ord

## ---- fig.align="center", fig.width=7.5, fig.height=5-------------------------
boxplot(sqrtOrgans[ord], col=cols[ord], 
        ylab=expression(sqrt("Survival time")), 
        main="Cancer treated by vitamin C", cex.axis=.6)

## -----------------------------------------------------------------------------
ordHam <-  hpaths(5, matrix = FALSE)
ordHam

## -----------------------------------------------------------------------------
# Get the test results
test <- with(cancer,
             pairwise.t.test(sqrt(Survival), Organ))
pvals <- test$p.value
pvals

## -----------------------------------------------------------------------------
# First construct a vector, removing NAs.
weights <- pvals[!is.na(pvals)]
weights <-edge2dist(weights)

## -----------------------------------------------------------------------------
weights <- as.matrix(weights)
rownames(weights) <- organNames
colnames(weights)<- rownames(weights)
weights

## -----------------------------------------------------------------------------
g <- mk_complete_graph(weights)

## ----fig.width=5, fig.height=5, fig.align='center'----------------------------
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
igplot(g,weights=TRUE,edge.label.color="black")

## ----fig.width=5, fig.height=5,fig.align='center', eval=FALSE-----------------
#  library(Rgraphviz)
#  ew <- round(unlist(edgeWeights(g)),3)
#  ew <- ew[setdiff(seq(along=ew), removedEdges(g))]
#  names(ew) <- edgeNames(g)
#  plot(g,  "circo",edgeAttrs=list(label=ew))

## -----------------------------------------------------------------------------
low2highEulord <- eulerian(weights); colnames(weights)[low2highEulord]
## or equivalently
eulerian(g)

## -----------------------------------------------------------------------------
bestHam <- order_best(weights)
colnames(weights)[bestHam]

## ---- fig.align="center", fig.width=7.5, fig.height=5-------------------------
boxplot(sqrtOrgans[low2highEulord], col=cols[ord], 
        ylab=expression(sqrt("Survival time")), 
        main="Cancer treated by vitamin C", cex.axis=.6)

## ----fig.width=5, fig.height=5, align='center'--------------------------------
aovOrgans <-  aov(sqrt(Survival) ~ Organ,data=cancer)
TukeyHSD(aovOrgans,conf.level = 0.95)

## -----------------------------------------------------------------------------
tuk <-TukeyHSD(aovOrgans,conf.level = 0.95)
ptuk <- tuk$Organ[,"p adj"]
dtuk <- as.matrix(edge2dist(ptuk))
rownames(dtuk)<- colnames(dtuk)<- organNames
g <- mk_complete_graph(weights)
eulerian(dtuk)

## ----fig.width=5, fig.height=5, align='center'--------------------------------
par(mar=c(3,8,3,3))
plot(TukeyHSD(aovOrgans,conf.level = 0.95),las=1,tcl = -.3)

## ----fig.align="center", fig.width=7.5, fig.height=5--------------------------

mc_plot(sqrtOrgans,aovOrgans,main="Pairwise comparisons of cancer types", 
        ylab="Sqrt Survival",col=cols,cex.axis=.6)

## ----fig.align="center", fig.width=7.5, fig.height=5--------------------------

suppressPackageStartupMessages(library(multcomp))
fitVitC <- glht(aovOrgans, linfct = mcp(Organ= "Tukey"))

# this gives confidence intervals without a family correction
confint(fitVitC,level=.99,calpha = univariate_calpha())

# for short, define
cifunction<- function(f, lev) 
  confint(f,level=lev,calpha = univariate_calpha())$confint

# calculate the confidence intervals
conf <- cbind(cifunction(fitVitC,.9),cifunction(fitVitC,.95)[,-1],cifunction(fitVitC,.99)[,-1])


mc_plot(sqrtOrgans,conf,path=low2highEulord,
        main="Pairwise comparisons of cancer types", 
        ylab="Sqrt Survival",col=cols,cex.axis=.6)




## -----------------------------------------------------------------------------
if (!requireNamespace("Sleuth3", quietly = TRUE)){
    install.packages("Sleuth3")
}
library(Sleuth3)
mice <- case0501
str(mice)
levels(mice$Diet)
# get rid of "/"
levels(mice$Diet) <-  c("NN85", "NR40", "NR50", "NP" ,   "RR50" ,"lopro")

## ---- fig.align="center", fig.width=6, fig.height=5---------------------------

life <- with(mice, split(Lifetime ,Diet))
cols <- rainbow_hcl(6, c = 50) 
boxplot(life, col=cols, 
        ylab="Lifetime", 
        main="Diet Restriction and Longevity")

## -----------------------------------------------------------------------------
aovMice   <- aov(Lifetime ~ Diet-1, data=mice)
fitMice <- glht(aovMice,
          linfct=c("DietNR50 - DietNN85 = 0", 
          "DietRR50  - DietNR50 = 0",
          "DietNR40  - DietNR50 = 0",
          "Dietlopro - DietNR50 = 0",
          "DietNN85  - DietNP   = 0")) 
  summary(fitMice,test=adjusted("none")) # No multiple comparison adjust.
  confint(fitMice, calpha = univariate_calpha()) # No adjustment

## -----------------------------------------------------------------------------
g <- new("graphNEL", nodes=names(life))

## ----fig.align='center', fig.width=5, fig.height=5----------------------------
fitMiceSum <- summary(fitMice,test=adjusted("none"))
pvalues <- fitMiceSum$test$pvalues
pvalues
# Extract  labels from the p-values for the edges
edgeLabs <- unlist(strsplit(names(pvalues), " - "))
edgeLabs <- matrix(substring(edgeLabs,5), nrow=2)
g <- addEdge(edgeLabs[1,], edgeLabs[2,], g,pvalues)

pos <- rbind(c(-1,0), c(0,-1), c(0,0), c(-2,0),c(1,0), c(0,1))
igplot(g, weights=TRUE, layout=pos,vertex.size=32)

## -----------------------------------------------------------------------------
eulerian(g)

## ----fig.align='center', fig.width=5, fig.height=5----------------------------
g1 <- addEdge("NR40","NP",g,1)
g1 <- addEdge("lopro","RR50",g1,1)
igplot(g1, weights=TRUE, layout=pos,vertex.size=32)
eulerian(g1)

## ---- fig.align="center", fig.width=6, fig.height=5---------------------------

eul <- eulerian(g1)
# make eul numeric
eul <- match(eul, names(life))

fitMice1 <- glht(aovMice, linfct = mcp(Diet= "Tukey"), calpha = univariate_calpha())

# need to construct the confidence intervals for all pairs
conf <- cbind(cifunction(fitMice1,.9),cifunction(fitMice1,.95)[,-1],cifunction(fitMice1,.99)[,-1])
# these comparisons are not relevant
conf[c(1,4,5,7,8,9,13,14,15),]<- NA

mc_plot(life,conf,path=eul,
        main="Diet Restriction and Longevity", 
        ylab="Lifetime",col=cols,cex.axis=.6)


