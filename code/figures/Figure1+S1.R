# Figure 1
# generate panels of Figure 1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plot tree

library(ape)
library(ggtree)

tree <- read.tree("data/tree/tree.nw")
names <- read.table("data/tree/tree_names_short.txt")
tax <- read.table("data/taxonomy/tree_taxonomy.txt", header=T, sep='\t')
tree$tip.label <- as.character(names[,1])

groupInfo <- tax$kingdom
group <- list()

u <- unique(groupInfo)
for (i in 1:length(u)){	
	c <- which(groupInfo == u[i])
	group[[i]] <- c
}
tree <- groupOTU(tree,group ) 



postscript("figures/Figure1/Fig1D.ps", width=5, height=5, paper="special", horizontal=T, onefile=F)

#ggtree(tree, layout='circular', aes(color=group), branch.length="none" )
ggtree(tree, layout='circular', aes(color=group), branch.length="none" ) + geom_tiplab2(size=1.8, family="Arial", aes(angle=angle)) + ggplot2::xlim(0, 40)

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot bechmark barplots

postscript("figures/Figure1/Fig1B.ps",  height=5, width=3, paper="special", horizontal=F, onefile=F)

par(mfrow=c(1,1))
par(mar=c(7,6,2,1))

# 472+2 total, add naked mole rat and killifish to the counts
pourcentage = c(90/474*100, 121/474*100, 298/474*100)

b <- barplot(pourcentage, ylim = c(0,75) , axes=F, main="", xlab="", ylab="", col=c("#888888", "#555555", "#222222"), border=F )
 
axis(1, at=b, labels=c( "Blast","RNAmmer", "Pipeline"), cex.axis=1.2, tick=F, las=2)    
axis(2, at=c(0, 25, 50, 75), cex.axis=1.3)
mtext("% identified 18S rRNA", side=2, line=3, cex=1.3)
 
dev.off()



postscript("figures/Figure1/Fig1C.ps", height=5, width=3, paper="special", horizontal=F, onefile=F)
 
sequences_rrna = c(474,298,216)

par(mfrow=c(1,1))
par(mar=c(7,6,2,1))

b <- barplot(sequences_rrna, ylim = c(0,500) , axes=F, main="", xlab="", ylab="", col=c("#888888", "#555555", "#222222"), border=F)

axis(1, at=b, labels=c("all", "with 18S", "final"), cex.axis=1.2, tick=F, las=2)    
axis(2, at=c(0, 250, 500), cex.axis=1.3)

mtext("Number of genomes", side=2, line=3, cex=1.3)

dev.off()

         




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Analyses of chaperone networks on tree
# Figure S1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



matrice_gaps <- as.matrix( read.table("data/tree/gaps_matrix_final.txt"), header=T )


postscript("figures/Supplement/FigS1A.ps", width=5, height=4, paper="special", horizontal=T, onefile=F)

par(mfrow=c(1,1))
par(mar=c(6,6,2,1))

plot(colMeans(matrice_gaps), type="p", pch=3, axes=F, xlab="", ylab="", main="", ylim=c(0, 0.75))

#axis(1)
axis(2, at=c(0, 0.25, 0.50, 0.75), labels=c(0, 25, 50, 75))
mtext("Species", side=1, line=1, cex=1.2)
mtext("Average % gaps", side=2, line=3, cex=1.2)

abline(h = 0.15 , col = "red")

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Convergence of MCMC PhyloBayes analysis

tree_mcmc_r1 <- read.table("data/tree/model_CATGTR_15_10000chains_r1.trace", header=T)
tree_mcmc_r2 <- read.table("data/tree/model_CATGTR_15_10000chains_r2.trace", header=T)
tree_mcmc_r3 <- read.table("data/tree/model_CATGTR_15_10000chains_r3.trace", header=T)


postscript("figures/Supplement/FigS1B.ps", height=4, width=8, paper="special", horizontal=F, onefile=T)

par(mfrow=c(1,2))
par(mar=c(6,6,2,1))

plot( tree_mcmc_r1$loglik[1:200] , type = "l" , col="grey30" , ylab = "" , xlab = "" ,axes = F, ylim=c(-110000 , -20000), lwd=3)
lines(tree_mcmc_r2$loglik[1:200] , col= "grey60", lwd=3)
lines(tree_mcmc_r3$loglik[1:200] , type = "l" , col="grey90" ,lwd=3)

axis(1, cex.axis = 1.2)
axis(2,at=c(-110000 ,-60000 , -20000),  cex.axis = 1.2)
mtext("log likelihood", side=2, line=2.5, cex=1.2)
mtext("Monte Carlo steps", side=1, line=2.5, cex=1.2)
legend(2000,-90000, c("Filtred (chain 1)", "Filtred (chain 2)", "Not filtred (chain 1)", "Not filtred (chain 2)"), lty=c(1,3), lwd=c(2.5,2.5),col=c("green4","springgreen" , "royalblue4" , "blue"), cex=.6, bty='n')



plot( tree_mcmc_r1$loglik , type = "l" , col="grey30" , ylab = "" , xlab = "" ,axes = F, ylim=c(-110000 , -20000), lwd=3)
lines(tree_mcmc_r2$loglik , col= "grey60", lwd=3)
lines(tree_mcmc_r3$loglik , type = "l" , col="grey90" ,lwd=3)

axis(1, cex.axis = 1.2)
axis(2,at=c(-110000 ,-60000 , -20000),  cex.axis = 1.2)
mtext("log likelihood", side=2, line=2.5, cex=1.2)
mtext("Monte Carlo steps", side=1, line=2.5, cex=1.2)


dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(proxy)

# chaperone networks
hsp <- read.csv("data/chaperones/hsp.txt", stringsAsFactors = F, sep='\t')
names <- read.table("data/tree/tree_list_nodes.txt")
chaperones <- c("Hsp20","Hsp40", "Hsp60", "Hsp70", "Hsp90", "Hsp100")
colnames(hsp) <- c("species", chaperones)

evol <- as.matrix(read.table("data/tree/tree_distancematrix.txt", sep='\t', header=T, row.names=1))



library(ape)

tree <- read.tree("data/tree/tree.nw")

c1 <- tax$kingdom == "Fungi"
c2 <- tax$kingdom == "Metazoa"
c3 <- tax$kingdom == "Viridiplantae"
c4 <- tax$kingdom == "none"


postscript("figures/Supplement/FigS1C.ps", width=6, height=4, paper="special", horizontal=T, onefile=F)

par(mfrow=c(1,1))
par(mar=c(4,6,4,1))

boxplot( c(evol[c1,c1]), c(evol[c1,-c1]),  c(evol[c2,c2]), c(evol[c2,-c2]),  c(evol[c3,c3]), c(evol[c3,-c3]),  c(evol[c4,c4]), c(evol[c4,-c4]) , ylim=c(0,1), at=c(0.2, 0.7, 1.4, 1.9, 2.6, 3.1, 3.8, 4.3), boxwex=0.4, xlim=c(0, 4.8) , col=c("#29ABE2", "#555555"), axes=F, xlab="", ylab="")

axis(3, at=c(0.45, 1.65, 2.85, 4.05), labels=c("Fungi", "Metazoa", "Plantae", "Protista"), tick=F, line=0)
axis(2, at=c(0, 0.5, 1), labels=c(0, 0.5, 1))
mtext("Inferred evolutionary distance", side=2, line=2.5)

legend("bottomleft", inset=c(0, -0.2), legend=c("within cluster", "outside cluster"), pch=15, col=c("#29ABE2", "#555555"), bty='n', xpd=TRUE)

dev.off()


