# Figure 1
# generate panels of Figure 1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plot tree

library(ape)
library(ggtree)

tree <- read.tree("../../data/Phylobayes/eucaryote_from_concensus_5_11.newick")
names <- read.table("../../data/list_organisms_final_tree_short.txt")
tax <- read.table("../../data/taxonomy/tree_taxonomy_corrected.txt", header=T, sep='\t')

tree$tip.label <- as.character(names[,1])


#groupInfo <- tax$phylum
groupInfo <- tax$kingdom

group <- list()

u <- unique(groupInfo)
for (i in 1:length(u)){	
	c <- which(groupInfo == u[i])
	group[[i]] <- c
}

tree <- groupOTU(tree,group ) 



postscript("../../figures/Figure1/Fig1D.ps", width=5, height=5, paper="special", horizontal=T, onefile=F)

ggtree(tree, layout='circular', aes(color=group), branch.length="none" ) + geom_tiplab2(size=1.8, family="Arial", aes(angle=angle)) + ggplot2::xlim(0, 18) #+ scale_colour_hue()

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot bechmark barplots

postscript("../../figures/Figure1/Fig1B.ps",  height=5, width=3, paper="special", horizontal=F, onefile=F)

par(mfrow=c(1,1))
par(mar=c(7,6,2,1))

pourcentage = c(15.84, 19.31 ,35.87)

b <- barplot(pourcentage, ylim = c(0,40) , axes=F, main="", xlab="", ylab="", col=c("#888888", "#555555", "#222222"), border=F )
 
axis(1, at=b, labels=c( "Blast","RNAmmer", "Combined"), cex.axis=1.2, tick=F, las=2)    
axis(2, at=c(0, 20, 40), cex.axis=1.3)
mtext("% identified 18S rRNA", side=2, line=3, cex=1.3)
 
dev.off()



postscript("../../figures/Figure1/Fig1C.ps", height=5, width=3, paper="special", horizontal=F, onefile=F)
 
sequences_rrna = c(524,188,158)

par(mfrow=c(1,1))
par(mar=c(7,6,2,1))

b <- barplot(sequences_rrna, ylim = c(0,550) , axes=F, main="", xlab="", ylab="", col=c("#222222", "#29ABE2", "red"), border=F)

axis(1, at=b, labels=c("all", "identified", "filtered"), cex.axis=1.2, tick=F, las=2)    
axis(2, at=c(0, 250, 500), cex.axis=1.3)

mtext("Number of genomes", side=2, line=3, cex=1.3)

dev.off()

         




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get taxonomy info

library(myTAI)


getinfo <- function(i){
	
	t <- taxonomy(names[i])

	ranks <- c(which(t$rank == "superkingdom"), which(t$rank == "kingdom"), which(t$rank == "phylum"), which(t$rank == "subphylum"), which(t$rank == "class"), which(t$rank == "order"), which(t$rank == "family") )

	classes <- t$name[ranks]

	output <- c(names[i], classes)

	output
}


out <- getinfo(1)

#df <-data.frame(out)
#names(df) <- c("name", "superkingdom", "kingdom", "phylum", "subphylum", "class", "order", "family")



for (i in 2:162){

out2 <- getinfo(i)

out <- rbind(out, out2)

}

df <- data.frame(out)
names(df) <- c("name", "superkingdom", "kingdom", "phylum", "subphylum", "class", "order", "family")
rownames(df) <- 1:nrow(df)
write.table(df,"../../data/taxonomy/tree_taxonomy.txt",sep="\t",row.names=FALSE, quote=F)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Analyses of chaperone networks on tree
# Figure S1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# This script plot matrix of gaps by three thresholds -  Yasmine Draceni October 2018
# Load data


matrice_gaps <- read.csv("../../data/pipeline_mapping_reads_genome/dataframe_gaps/matrice_gaps_17_sep.csv")



postscript("../../figures/Supplement/FigS1A.ps", width=5, height=4, paper="special", horizontal=T, onefile=F)

par(mfrow=c(1,1))
par(mar=c(6,6,2,1))

plot(matrice_gaps$somme[c(1:176)], type="p", pch=3, axes=F, xlab="", ylab="", main="", ylim=c(0, 6500))

#axis(1)
axis(2, at=c(0, 2000, 4000, 6000))
mtext("Species", side=1, line=1, cex=1.2)
mtext("Total number of gaps", side=2, line=3, cex=1.2)

abline(h =1000 , col = "red")
abline(h = 2000 , col = "red")
abline(h=3000 , col = "red")

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Convergence of MCMC PhyloBayes analysis


# Load data
gaps_2000.fungiS_CATGTR.4000_chains_B.3.10 <- read.delim("../../data/Phylobayes/gaps_2000-fungiS_CATGTR-4000_chains_B-3-10.trace")
#View(gaps_2000.fungiS_CATGTR.4000_chains_A.3.10)
Eukaryote_with_gaps_v2 <- read.delim("../../data/Phylobayes/Eukaryote_with_gaps_v2.trace", header=FALSE, comment.char="#")
Eukaryote_with_gaps <- read.delim("../../data/Phylobayes/Eukaryote_with_gaps.trace", header=FALSE, comment.char="#")
gaps_2000.fungiS_CATGTR.4000_chains_A.3.10 <- read.delim("../../data/Phylobayes/gaps_2000-fungiS_CATGTR-4000_chains_A-3-10.trace")

# save it in file
postscript("../../figures/Supplement/FigS1B.ps", height=4, width=4, paper="special", horizontal=F, onefile=T)

par(mfrow=c(1,1))
par(mar=c(6,6,2,1))

plot((gaps_2000.fungiS_CATGTR.4000_chains_A.3.10$loglik) , type = "l" , col="green4" , ylab = "" , xlab = "" ,axes = F, ylim=c(-110000 , -20000), lwd=3)
lines(gaps_2000.fungiS_CATGTR.4000_chains_B.3.10$loglik , col= "springgreen", lwd=3, lty=3)
lines(Eukaryote_with_gaps$V4 , type = "l" , col="royalblue4" , ylab = "loglik" , xlab = "")
lines(Eukaryote_with_gaps_v2$V4 , col="blue", lty=3)
axis(1, cex.axis = 1.2)
axis(2,at=c(-110000 ,-60000 , -20000),  cex.axis = 1.2)
mtext("log likelihood", side=2, line=2.5, cex=1.2)
mtext("Monte Carlo steps", side=1, line=2.5, cex=1.2)
legend(2000,-90000, c("Filtred (chain 1)", "Filtred (chain 2)", "Not filtred (chain 1)", "Not filtred (chain 2)"), lty=c(1,3), lwd=c(2.5,2.5),col=c("green4","springgreen" , "royalblue4" , "blue"), cex=.6, bty='n')

dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(proxy)

# chaperone networks
hsp <- read.csv("../../data/Heatmap_HSP/HSP_phylogeny_5_november.csv")
names <- read.table("../../data/list_organisms_final_tree_short.txt")
chaperones <- c("Hsp20","Hsp40", "Hsp60", "Hsp70", "Hsp90", "Hsp100")

hsp <- hsp[,-1]
colnames(hsp) <- c("species", chaperones, "proteome")


evol <- as.matrix(read.table("../../data/Heatmap_HSP/dataframe_evolution_distances.csv", sep=',', header=T, row.names=1))


chap.cor  <- matrix(NA, nrow=nrow(hsp), ncol=nrow(hsp))
chap.dist <- matrix(NA, nrow=nrow(hsp), ncol=nrow(hsp))

for (i in 1:nrow(hsp)){	
	for (j in 1:nrow(hsp)){
		
		a <- as.numeric(hsp[i,3:8])
		b <- as.numeric(hsp[j,3:8])
		
		chap.cor[i,j]  <- cor(a,b)
		chap.dist[i,j] <- as.numeric( simil(list(a,b),method="Euclidean") )
	}	
}



library(ape)

tree <- read.tree("../../data/Phylobayes/eucaryote_from_concensus_5_11.newick")
names <- read.table("../../data/list_organisms_final_tree_short.txt")
tree$tip.label <- as.character(names[,1])


c1 <- tax$phylum == "Ascomycota"
c2 <- tax$phylum == "Basidiomycota"
c3 <- tax$kingdom == "Protista"
c4 <- tax$kingdom == "Animalia"
c5 <- tax$kingdom == "Plantae"


postscript("../../figures/Supplement/FigS1C.ps", width=6, height=4, paper="special", horizontal=T, onefile=F)

par(mfrow=c(1,1))
par(mar=c(4,6,4,1))

boxplot( c(evol[c1,c1]), c(evol[c1,-c1]),  c(evol[c2,c2]), c(evol[c2,-c2]),  c(evol[c3,c3]), c(evol[c3,-c3]),  c(evol[c4,c4]), c(evol[c4,-c4]) ,  c(evol[c5,c5]), c(evol[c5,-c5]) , ylim=c(0,1), at=c(0.2, 0.7, 1.4, 1.9, 2.6, 3.1, 3.8, 4.3, 5, 5.5), boxwex=0.4, xlim=c(0, 5.7) , col=c("#29ABE2", "#555555"), axes=F, xlab="", ylab="")

axis(3, at=c(0.45, 1.65, 2.85, 4.05, 5.25), labels=c("FungiA", "FungiB", "Protists", "Animalia", "Plantae"), tick=F, line=0)
axis(2, at=c(0, 0.5, 1), labels=c(0, 0.5, 1))
mtext("Inferred evolutionary distance", side=2, line=2.5)

legend("bottomleft", inset=c(0, -0.2), legend=c("within cluster", "outside cluster"), pch=15, col=c("#29ABE2", "#555555"), bty='n', xpd=TRUE)

dev.off()


