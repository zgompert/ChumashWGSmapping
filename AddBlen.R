## adds ml branch lengths to the tree
library(phangorn)
library(ape)

packageVersion("phangorn")
#[1] ‘2.12.1’
packageVersion("ape")
#[1] ‘5.8’

## read tree, just topology, from caster
tre<-read.tree("topo_cout_tree")

## read the alignment and convert to phyDat object
aln<-read.dna("clean.mfa",format="fasta")
alnPh<-phyDat(aln,type="DNA")

## add arbitrary branch lengths as a starting point
treb<-compute.brlen(tre,1)

## compute the initial likelihood and then optimize branch lenghts under JC
## given the topoloty
fit<-pml(tree=treb,data=alnPh)
fit_optimized <- optim.pml(fit, model = "JC", optEdge = TRUE, optGamma = FALSE, optInv = FALSE)

## I might want to increase these slightly so they are overestimates
fit_optimized$tree$edge.length<-fit_optimized$tree$edge.length * 1.5

## write the tree
write.tree(phy=fit_optimized$tree,file="topoB_cout_tree")

## root on refugio and write again
refugio<-fit_optimized$tree$tip.label[9:16]
trr<-root(fit_optimized$tree,outgroup=refugio)
write.tree(phy=trr,file="topoRoot_cout_tree")

