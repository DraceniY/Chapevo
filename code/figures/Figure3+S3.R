# Analyses of chaperone counts

library(proxy)

# chaperone counts
hsp <- read.csv("data/chaperones/hsp.txt", sep='\t', stringsAsFactors = F)
names <- read.table("data/tree/tree_list_nodes.txt")
chaperones <- c("Hsp20","Hsp40", "Hsp60", "Hsp70", "Hsp90", "Hsp100")
tax <- read.table("data/taxonomy/tree_taxonomy.txt", header=T, sep='\t')
colnames(hsp) <- c("species", chaperones)

evol <- as.matrix(read.table("data/tree/tree_distancematrix.txt", sep='\t', header=T, row.names=1))

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

c1 <- tax$kingdom == "Fungi"
c2 <- tax$kingdom == "Metazoa"
c3 <- tax$kingdom == "Viridiplantae"


df.plantae <- data.frame( 	species=names(rowMeans(evol[c3,c3])), 
								evodist=as.vector(rowMeans(evol[c3,c3])) ,
								correlation=as.vector(rowMeans(chap.cor[c3,c3])),
								class=rep("3 - Plantae", sum(c3, na.rm=T))   )

df.fungi <- data.frame( 	species=names(rowMeans(evol[c1,c1])), 
								evodist=as.vector(rowMeans(evol[c1,c1])) ,
								correlation=as.vector(rowMeans(chap.cor[c1,c1])),
								class=rep("1 - Fungi", sum(c1, na.rm=T))   )

df.animalia <- data.frame( 	species=names(rowMeans(evol[c2,c2])), 
								evodist=as.vector(rowMeans(evol[c2,c2])) ,
								correlation=as.vector(rowMeans(chap.cor[c2,c2])),
								class=rep("2 - Metazoa", sum(c2, na.rm=T))   )

df.evol <- rbind(df.fungi, df.animalia, df.plantae)


postscript("figures/Figure3/Fig3A.ps", width=6, height=2.5, paper="special", horizontal=T, onefile=F)

p <- ggplot(df.evol, aes(x=evodist, y=correlation) ) + geom_point() + theme_classic( ) 
p + facet_grid(cols = vars(class), scales = "free_x" ) 

dev.off()




# DETAILS ON CHAPERONE NETWORK OF INTERESTING ANIMALS
data <- as.matrix( hsp[,c(2:7)])
data <- data.frame(cbind(names,data))
colnames(data) <- c("species", chaperones)

c2 <- tax$kingdom == "Metazoa"
hsp.animalia <- data[c2,]

oyster <- hsp.animalia$species == "Crassostrea_gigas"
nmrat <- 	hsp.animalia$species == "Heterocephalus_glaber"
killi <- 	hsp.animalia$species == "Nothobranchius_furzeri"
urchin <- hsp.animalia$species == "Strongylocentrotus_purpuratus"

hsp.animalia.summary <- data.frame(class=colnames(hsp.animalia)[-1], mean=as.vector(colMeans(hsp.animalia[,-1])), sd=as.vector( apply(hsp.animalia[,-1], 2, sd) ) )
hsp.animalia.summary$ix <- 1:6
rownames(hsp.animalia.summary) <- 1:6

hsp.animalia.summary$oyster <- as.numeric(hsp.animalia[oyster,-1])
hsp.animalia.summary$nmrat <- as.numeric(hsp.animalia[nmrat,-1])
hsp.animalia.summary$killi <- as.numeric(hsp.animalia[killi,-1])
hsp.animalia.summary$urchin <- as.numeric(hsp.animalia[urchin,-1])


postscript("figures/Figure3/Fig3B.ps", width=6, height=4, paper="special", horizontal=T, onefile=F)

ggplot(hsp.animalia.summary, aes(x=ix, y=mean)) +
	geom_ribbon(aes(ymax=mean+sd, ymin=mean-sd), fill="#eeeeee", alpha=1) +
	geom_point(aes(x=ix, y=mean), shape=16) +
	theme_classic() +
	geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0, position=position_dodge(0.9)) +
	scale_y_continuous(expand = c(0, 0), limits=c(0,88)) +
	scale_x_continuous(breaks=c(1,2,3,4,5,6), labels=c("Hsp20", "Hsp40", "Hsp60", "Hsp70","Hsp90", "Hsp100")) +
	geom_point(aes(x=ix, y=oyster), col="#025074", size=4, shape=16) + geom_line(aes(x=ix, y=oyster), col='#025074') + 
	geom_point(aes(x=ix, y=nmrat), col="#b54a03", size=4, shape=15) + geom_line(aes(x=ix, y=nmrat), col='#b54a03') +
  geom_point(aes(x=ix, y=urchin), col="#991e7b", size=4, shape=17) + geom_line(aes(x=ix, y=urchin), col='#991e7b') +
  geom_point(aes(x=ix, y=killi), col="#03b582", size=4, shape=18) + geom_line(aes(x=ix, y=killi), col='#03b582') +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14) )
    
dev.off()



hsp.animalia2 <- melt(hsp.animalia)

postscript("figures/Figure3/Fig3C.ps", width=6, height=4, paper="special", horizontal=T, onefile=F)

library(cowplot)

nmrat <- hsp.animalia[hsp.animalia$species=="Heterocephalus_glaber",]
killi <- hsp.animalia[hsp.animalia$species=="Nothobranchius_furzeri",]
oyster <- hsp.animalia[hsp.animalia$species=="Crassostrea_gigas",]
urchin <- hsp.animalia[hsp.animalia$species=="Strongylocentrotus_purpuratus",]

plot.hsp40 <- ggplot(subset(hsp.animalia2, variable == "Hsp40"), aes(y=value) ) + 	
	geom_boxplot(aes(fill="#bbbbbb"), show.legend = FALSE) + 
	theme_classic( ) + 
	coord_flip() + 
	scale_fill_manual(values = c("grey")) +
	theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
	labs(x = "Hsp40", y="Counts") +

	geom_point(aes(x=0.2, y=oyster$Hsp40), col="#025074", shape=16, size=6) + 			#oyster
	geom_point(aes(x=0.2, y=nmrat$Hsp40), col="#b54a03", shape=15, size=6) +			#nmrat
  geom_point(aes(x=0.2, y=urchin$Hsp40), col="#991e7b", shape=17, size=6) +				#urchin
	geom_point(aes(x=0.2, y=killi$Hsp40), col="#03b582", shape=18, size=6) 			#killi
	

plot.hsp70 <- ggplot(subset(hsp.animalia2, variable == "Hsp70"), aes(y=value) ) + 	
	geom_boxplot(aes(fill="#bbbbbb"), show.legend = FALSE) + 
	theme_classic( ) + 
	coord_flip() + 
	scale_fill_manual(values = c("grey")) +
	theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
	labs(x = "Hsp70", y="Counts") +
	
	geom_point(aes(x=0.2, y=oyster$Hsp70), col="#025074", shape=16, size=6) + 			#oyster
	geom_point(aes(x=0.2, y=nmrat$Hsp70), col="#b54a03", shape=15, size=6) +			#nmrat
	geom_point(aes(x=0.2, y=urchin$Hsp70), col="#991e7b", shape=17, size=6) +				#urchin
  geom_point(aes(x=0.2, y=killi$Hsp70), col="#03b582", shape=18, size=6) 			#killi


plot.hsp90 <- ggplot(subset(hsp.animalia2, variable == "Hsp90"), aes(y=value) ) + 	
	geom_boxplot(aes(fill="#bbbbbb"), show.legend = FALSE) + 
	theme_classic( ) + 
	coord_flip() + 
	scale_fill_manual(values = c("grey")) +
	theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
	labs(x = "Hsp90", y="Counts") +
	
	geom_point(aes(x=0.2, y=oyster$Hsp90), col="#025074", shape=16, size=6) + 			#oyster
	geom_point(aes(x=0.2, y=nmrat$Hsp90), col="#b54a03", shape=15, size=6) +			#nmrat
	geom_point(aes(x=0.2, y=urchin$Hsp90), col="#991e7b", shape=17, size=6) +
  geom_point(aes(x=0.2, y=killi$Hsp90), col="#03b582", shape=18, size=6) 			#killi


plot_grid(plot.hsp40, plot.hsp70, plot.hsp90, labels ="", ncol = 1, align = 'v')


dev.off()





postscript("figures/Supplement/FigS3.ps", width=6, height=4, paper="special", horizontal=T, onefile=F)

library(cowplot)

nmrat <- hsp.animalia[hsp.animalia$species=="Heterocephalus_glaber",]
killi <- hsp.animalia[hsp.animalia$species=="Nothobranchius_furzeri",]
oyster <- hsp.animalia[hsp.animalia$species=="Crassostrea_gigas",]
urchin <- hsp.animalia[hsp.animalia$species=="Strongylocentrotus_purpuratus",]


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
	geom_point(aes(x=0.2, y=urchin$Hsp20), col="#991e7b", shape=17, size=6) 				#urchin


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
	geom_point(aes(x=0.2, y=urchin$Hsp60), col="#991e7b", shape=17, size=6) 				#urchin


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
	geom_point(aes(x=0.2, y=urchin$Hsp100), col="#991e7b", shape=17, size=6)


plot_grid(plot.hsp20, plot.hsp60, plot.hsp100, labels ="", ncol = 1, align = 'v')


dev.off()


