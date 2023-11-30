########
## Co-phylogeny plot
## Gutierrez Garcia-K et al, 2024
## The genetic basis of Lactobacillus-host specificity for the commensal niche in Drosophila melanogaster revealed through live imaging of colonization dynamics
########


# Installing packages 

install.packages(c("ape", "ggtree", "phytools"))
library(ape)
library(ggtree)
library(phytools)

# Changing the path

getwd( )
setwd('/your/path/')
getwd( )

# Loading the phylogenies: Core tree and aSec Tree both trees are collapsed and Midpoint rooted

CoreTree <- read.tree("41OrgCoreGenome-midpoint-collapsed.tree")
CoreTree

aSecTree <- read.tree("aSecSystem-Cornell-midpoint-collapsed.tree")
aSecTree

# Plotting the trees 

ggtree(CoreTree, ladderize = TRUE, right = TRUE) + geom_tiplab(size= 2) + geom_nodelab(color='red', size =5)
ggtree(aSecTree,ladderize = TRUE, right = TRUE) + geom_tiplab(size= 2) + geom_nodelab(color='red', size =5)

# Checking the bootstrap
bootstrap_values <- CoreTree$node.label
bootstrap_valuesB <- aSecTree$node.label

# Showing the first 41 values 
head(bootstrap_values, 41)
head(bootstrap_valuesB,41)


# Converting the collapsed trees to a phylo object 

CoreTree_phylo <- as.phylo(CoreTree)
plotTree(CoreTree_phylo, fsize=0.6)

aSecTree_phylo <- as.phylo(aSecTree)
plotTree(aSecTree_phylo, fsize=0.6)

# Tip labels

CoreTree_phylo$tip.label
aSecTree_phylo$tip.label

# Assocciation matrix between CoreTree and aSec tree

assoc.rot<-read.csv("CoreGenome-aSec-assoc.csv")
assoc.rot

# Rotate both trees

obj.bothrorate<-cophylo(CoreTree_phylo,aSecTree_phylo,assoc=assoc.rot,
                        rotate=TRUE)
obj.bothrorate

plot(obj.bothrorate,link.type="curved",link.lwd=3,link.lty="solid",
     link.col=make.transparent("grey",0.25),fsize=0.8, mar = c(1.1,1.1,4.1,1.1))

title(main = "rotate both", font.main =3)
#nodelabels.cophylo(which="left",frame="none",cex=1) # frame="circle", "none"
#nodelabels.cophylo(which="right",frame="none",cex=0.8)

# Numbering the nodes 1-41 

plot(obj.bothrorate,link.type="curved",link.lwd=3,link.lty="solid",
     link.col=make.transparent("grey",0.25),fsize=0.8, mar = c(1.1,1.1,4.1,1.1))

title(main = "rotate both", font.main =3)

textCoreTree <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 41)
nodelabels.cophylo(which = "left", frame="none", cex=1, text = textCoreTree, node = 1:41 + Ntip(CoreTree_phylo),
                   adj = c(-1, 0.5))
nodelabels.cophylo(which = "right", frame="none", cex=1, text = textArbolA, node = 1:41 + Ntip(CoreTree_phylo), adj = c(1, 0.5))


# Adding the number of nodes manually CoreGenome Tree

plot(obj.bothrorate,link.type="curved",link.lwd=3,link.lty="solid",
     link.col=make.transparent("grey",0.25),fsize=0.8, mar = c(1.1,1.1,4.1,1.1))

nodelabels.cophylo(
  which = "left",
  frame = "none",
  cex = 0.8,
  text = c("100", "100", "100", "100", "100", "100", "100", "100", "100","100", "100", "100", "100", "100", "100", "100"), # Adding the bootstrap to the nodes 
  node = c(1, 2, 20, 27, 25, 27, 28, 26, 30, 31, 29, 32, 33, 35, 36, 34) + Ntip(CoreTree_phylo),
  adj = c(-0.1, 0.5)  # Adjust the text position 
)

# Adding the number of nodes manually aSec Tree

plot(obj.bothrorate,link.type="curved",link.lwd=3,link.lty="solid",
     link.col=make.transparent("grey",0.9),fsize=0.8, mar = c(1.1,1.1,4.1,1.1))

nodelabels.cophylo(
  which = "right",
  frame = "none",
  cex = 0.9,
  text = c("87", "100", "100", "75", "100", "83", "100", "100", "72", "83", "100"), # Adding the bootstrap to the nodes 
  node = c(30, 29, 28, 17, 27, 26, 25, 23, 18, 11, 24) + Ntip(aSecTree_phylo),
  adj = c(-0.1, 0.1)  # Adjust the text position 
)

# Final tree

# Abre un dispositivo gráfico para guardar en formato PNG
png("Co-phylogeny-Collapsed.png", width = 800, height = 600)  # Adjust the width and height 
pdf("Co-phylogeny-Collapsed.pdf", width = 8, height = 6)

plot(obj.bothrorate,link.type="curved",link.lwd=3,link.lty="solid",
     link.col=make.transparent("grey",0.9),fsize=0.8, mar = c(1.1,1.1,4.1,1.1))

nodelabels.cophylo(
  which = "left",
  frame = "none",
  cex = 0.8,
  text = c("100", "100", "100", "100", "100", "100", "100", "100", "100","100", "100", "100", "100", "100", "100", "100"), # Adding the bootstrap to the nodes 
  node = c(1, 2, 20, 27, 25, 27, 28, 26, 30, 31, 29, 32, 33, 35, 36, 34) + Ntip(CoreTree_phylo),
  adj = c(-0.1, 0.5) # Adjust the text position 
)

nodelabels.cophylo(
  which = "right",
  frame = "none",
  cex = 0.9,
  text = c("87", "100", "100", "75", "100", "83", "100", "100", "72", "83", "100"), # Adding the bootstrap to the nodes 
  node = c(30, 29, 28, 17, 27, 26, 25, 23, 18, 11, 24) + Ntip(aSecTree_phylo),
  adj = c(-0.1, 0.1)  # Adjust the text position 
)

# Close the graphic
dev.off()



#  Patristic distances

patristicdistA.original <- cophenetic.phylo(CoreTree)
patristicdistA.original
patristicdistB.original <- cophenetic.phylo(aSecTree)
patristicdistB.original

matrixAB <- read.csv("MatrixParaFit.csv", header= TRUE, row.names = 1)

parafit(patristicdistA.original, patristicdistB.original, matrixAB,  nperm = 999, test.links = TRUE, seed =NULL, correction = "cailliez", silent = FALSE)
#Global test:  ParaFitGlobal = 0.7466382 , p-value = 0.001 ( 999 permutations) 
