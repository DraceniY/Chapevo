# Analyses of chaperone networks on tree

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
		
		a <- as.numeric(hsp[i,2:7])
		b <- as.numeric(hsp[j,2:7])
		
		chap.cor[i,j]  <- cor(a,b)
		chap.dist[i,j] <- as.numeric( simil(list(a,b),method="Euclidean") )
	}	
}

c1 <- tax$phylum == "Ascomycota"
c2 <- tax$phylum == "Basidiomycota"
c3 <- tax$kingdom == "Protista"
c4 <- tax$kingdom == "Animalia"
c5 <- tax$kingdom == "Plantae"



df.plantae <- data.frame( 	species=names(rowMeans(evol[c5,c5])), 
								evodist=as.vector(rowMeans(evol[c5,c5])) ,
								correlation=as.vector(rowMeans(chap.cor[c5,c5])),
								class=rep("3 - Plantae", sum(c5))   )

df.fungiA <- data.frame( 	species=names(rowMeans(evol[c1,c1])), 
								evodist=as.vector(rowMeans(evol[c1,c1])) ,
								correlation=as.vector(rowMeans(chap.cor[c1,c1])),
								class=rep("1 - FungiA", sum(c1))   )

df.animalia <- data.frame( 	species=names(rowMeans(evol[c4,c4])), 
								evodist=as.vector(rowMeans(evol[c4,c4])) ,
								correlation=as.vector(rowMeans(chap.cor[c4,c4])),
								class=rep("2 - Animalia", sum(c4))   )


df.evol <- rbind(df.fungiA, df.animalia, df.plantae)


postscript("../../figures/Figure3/Fig3A.ps", width=6, height=2.5, paper="special", horizontal=T, onefile=F)

p <- ggplot(df.evol, aes(x=evodist, y=correlation) ) + geom_point() + theme_classic( ) 
p + facet_grid(cols = vars(class), scales = "free_x" ) 

dev.off()



# ------------------

# DETAILS ON CHAPERONE NETWORK OF INTERESTING ANIMALS

hsp <- read.csv("../../data/Heatmap_HSP/HSP_phylogeny_5_november.csv")
names <- read.table("../../data/list_organisms_final_tree_short.txt")
chaperones <- c("Hsp20","Hsp40", "Hsp60", "Hsp70", "Hsp90", "Hsp100")

data <- as.matrix( hsp[,c(3:8)])
data <- data.frame(cbind(names,data))
colnames(data) <- c("species", chaperones)

hsp.animalia <- data[c4,]

oyster <- 	hsp.animalia$species == "C.gigas"
nmrat <- 	hsp.animalia$species == "H.glaber"
killi <- 	hsp.animalia$species == "N.furzeri"
sponge <- 	hsp.animalia$species == "S.purpuratus"
amphi <- 	hsp.animalia$species == "A.queenslandica"
loa <- 		hsp.animalia$species == "L.loa"
worm <- 	hsp.animalia$species == "C.elegans"



hsp.animalia.summary <- data.frame(class=colnames(hsp.animalia)[-1], mean=as.vector(colMeans(hsp.animalia[,-1])), sd=as.vector( apply(hsp.animalia[,-1], 2, sd) ) )
hsp.animalia.summary$ix <- 1:6
rownames(hsp.animalia.summary) <- 1:6

hsp.animalia.summary$oyster <- as.numeric(hsp.animalia[oyster,-1])
hsp.animalia.summary$nmrat <- as.numeric(hsp.animalia[nmrat,-1])
hsp.animalia.summary$killi <- as.numeric(hsp.animalia[killi,-1])
hsp.animalia.summary$sponge <- as.numeric(hsp.animalia[sponge,-1])
hsp.animalia.summary$amphi <- as.numeric(hsp.animalia[amphi,-1])
hsp.animalia.summary$loa <- as.numeric(hsp.animalia[loa,-1])
hsp.animalia.summary$worm <- as.numeric(hsp.animalia[worm,-1])


postscript("../../figures/Figure3/Fig3B.ps", width=6, height=4, paper="special", horizontal=T, onefile=F)

ggplot(hsp.animalia.summary, aes(x=ix, y=mean)) +
	geom_ribbon(aes(ymax=mean+sd, ymin=mean-sd), fill="#eeeeee", alpha=1) +
	geom_point(aes(x=ix, y=mean), shape=16) +
	theme_classic() +
	geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0, position=position_dodge(0.9)) +
	scale_y_continuous(expand = c(0, 0), limits=c(0,88)) +
	scale_x_continuous(breaks=c(1,2,3,4,5,6), labels=c("Hsp20", "Hsp40", "Hsp60", "Hsp70","Hsp90", "Hsp100")) +
	geom_point(aes(x=ix, y=oyster), col="#025074", size=4, shape=16) + geom_line(aes(x=ix, y=oyster), col='#025074') + 
	geom_point(aes(x=ix, y=nmrat), col="#b54a03", size=4, shape=15) + geom_line(aes(x=ix, y=nmrat), col='#b54a03') +
	geom_point(aes(x=ix, y=killi), col="#03b582", size=4, shape=18) + geom_line(aes(x=ix, y=killi), col='#03b582') +
    geom_point(aes(x=ix, y=sponge), col="#991e7b", size=4, shape=17) + geom_line(aes(x=ix, y=sponge), col='#991e7b') +
    
    theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14) )
    
   



dev.off()



hsp.animalia2 <- melt(hsp.animalia)


postscript("../../figures/Figure3/Fig3C.ps", width=6, height=4, paper="special", horizontal=T, onefile=F)

library(cowplot)

nmrat <- hsp.animalia[hsp.animalia$species=="H.glaber",]
killi <- hsp.animalia[hsp.animalia$species=="N.furzeri",]
oyster <- hsp.animalia[hsp.animalia$species=="C.gigas",]
sponge <- hsp.animalia[hsp.animalia$species=="S.purpuratus",]

plot.hsp40 <- ggplot(subset(hsp.animalia2, variable == "Hsp40"), aes(y=value) ) + 	
	geom_boxplot(aes(fill="#bbbbbb"), show.legend = FALSE) + 
	theme_classic( ) + 
	coord_flip() + 
	scale_fill_manual(values = c("grey")) +
	theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
	labs(x = "Hsp40", y="Counts") +

	geom_point(aes(x=0.2, y=oyster$Hsp40), col="#025074", shape=16, size=6) + 			#oyster
	geom_point(aes(x=0.2, y=nmrat$Hsp40), col="#b54a03", shape=15, size=6) +			#nmrat
	geom_point(aes(x=0.2, y=killi$Hsp40), col="#03b582", shape=18, size=6) +			#killi
	geom_point(aes(x=0.2, y=sponge$Hsp40), col="#991e7b", shape=17, size=6) 				#sponge


plot.hsp70 <- ggplot(subset(hsp.animalia2, variable == "Hsp70"), aes(y=value) ) + 	
	geom_boxplot(aes(fill="#bbbbbb"), show.legend = FALSE) + 
	theme_classic( ) + 
	coord_flip() + 
	scale_fill_manual(values = c("grey")) +
	theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
	labs(x = "Hsp70", y="Counts") +
	
	geom_point(aes(x=0.2, y=oyster$Hsp70), col="#025074", shape=16, size=6) + 			#oyster
	geom_point(aes(x=0.2, y=nmrat$Hsp70), col="#b54a03", shape=15, size=6) +			#nmrat
	geom_point(aes(x=0.2, y=killi$Hsp70), col="#03b582", shape=18, size=6) +			#killi
	geom_point(aes(x=0.2, y=sponge$Hsp70), col="#991e7b", shape=17, size=6) 				#sponge


plot.hsp90 <- ggplot(subset(hsp.animalia2, variable == "Hsp90"), aes(y=value) ) + 	
	geom_boxplot(aes(fill="#bbbbbb"), show.legend = FALSE) + 
	theme_classic( ) + 
	coord_flip() + 
	scale_fill_manual(values = c("grey")) +
	theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
	labs(x = "Hsp90", y="Counts") +
	
	geom_point(aes(x=0.2, y=oyster$Hsp90), col="#025074", shape=16, size=6) + 			#oyster
	geom_point(aes(x=0.2, y=nmrat$Hsp90), col="#b54a03", shape=15, size=6) +			#nmrat
	geom_point(aes(x=0.2, y=killi$Hsp90), col="#03b582", shape=18, size=6) +			#killi
	geom_point(aes(x=0.2, y=sponge$Hsp90), col="#991e7b", shape=17, size=6)


plot_grid(plot.hsp40, plot.hsp70, plot.hsp90, labels ="", ncol = 1, align = 'v')


dev.off()










postscript("../../figures/Supplement/FigS3.ps", width=6, height=4, paper="special", horizontal=T, onefile=F)

library(cowplot)

nmrat <- hsp.animalia[hsp.animalia$species=="H.glaber",]
killi <- hsp.animalia[hsp.animalia$species=="N.furzeri",]
oyster <- hsp.animalia[hsp.animalia$species=="C.gigas",]
sponge <- hsp.animalia[hsp.animalia$species=="S.purpuratus",]



plot.hsp20 <- ggplot(subset(hsp.animalia2, variable == "Hsp20"), aes(y=value) ) + 	
	geom_boxplot(aes(fill="#bbbbbb"), show.legend = FALSE) + 
	theme_classic( ) + 
	coord_flip() + 
	scale_fill_manual(values = c("grey")) +
	theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x = element_text(size=16)) +
	labs(x = "Hsp20", y="Counts") +

	geom_point(aes(x=0.2, y=oyster$Hsp20), col="#025074", shape=16, size=6) + 			#oyster
	geom_point(aes(x=0.2, y=nmrat$Hsp20), col="#b54a03", shape=15, size=6) +			#nmrat
	geom_point(aes(x=0.2, y=killi$Hsp20), col="#03b582", shape=18, size=6) +			#killi
	geom_point(aes(x=0.2, y=sponge$Hsp20), col="#991e7b", shape=17, size=6) 				#sponge


plot.hsp60 <- ggplot(subset(hsp.animalia2, variable == "Hsp60"), aes(y=value) ) + 	
	geom_boxplot(aes(fill="#bbbbbb"), show.legend = FALSE) + 
	theme_classic( ) + 
	coord_flip() + 
	scale_fill_manual(values = c("grey")) +
	theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x = element_text(size=16)) +
	labs(x = "Hsp60", y="Counts") +
	
	geom_point(aes(x=0.2, y=oyster$Hsp60), col="#025074", shape=16, size=6) + 			#oyster
	geom_point(aes(x=0.2, y=nmrat$Hsp60), col="#b54a03", shape=15, size=6) +			#nmrat
	geom_point(aes(x=0.2, y=killi$Hsp60), col="#03b582", shape=18, size=6) +			#killi
	geom_point(aes(x=0.2, y=sponge$Hsp60), col="#991e7b", shape=17, size=6) 				#sponge


plot.hsp100 <- ggplot(subset(hsp.animalia2, variable == "Hsp100"), aes(y=value) ) + 	
	geom_boxplot(aes(fill="#bbbbbb"), show.legend = FALSE) + 
	theme_classic( ) + 
	coord_flip() + 
	scale_fill_manual(values = c("grey")) +
	theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x = element_text(size=16)) +
	labs(x = "Hsp100", y="Counts") +
	
	geom_point(aes(x=0.2, y=oyster$Hsp100), col="#025074", shape=16, size=6) + 			#oyster
	geom_point(aes(x=0.2, y=nmrat$Hsp100), col="#b54a03", shape=15, size=6) +			#nmrat
	geom_point(aes(x=0.2, y=killi$Hsp100), col="#03b582", shape=18, size=6) +			#killi
	geom_point(aes(x=0.2, y=sponge$Hsp100), col="#991e7b", shape=17, size=6)


plot_grid(plot.hsp20, plot.hsp60, plot.hsp100, labels ="", ncol = 1, align = 'v')


dev.off()


