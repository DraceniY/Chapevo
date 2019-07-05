# Analyses of chaperone counts

library(proxy)
library(reshape)
library(ggplot2)

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggdendro)

hsp <- read.csv("data/chaperones/hsp.txt", sep='\t', stringsAsFactors = F)

hsp.name <- hsp[,1]
hsp.norm <- as.matrix(hsp[,-1])
hsp.norm <- hsp.norm / rowSums(hsp.norm)

tax <- read.table("data/taxonomy/tree_taxonomy.txt", header=T, sep='\t')
c1 <- tax$kingdom == "Fungi"
c2 <- tax$kingdom == "Metazoa"
c3 <- tax$kingdom == "Viridiplantae"


cosinesimilarity <- function(a, b){
	d = sum(a*b)/sqrt(sum(a^2)*sum(b^2))
	return(d)
}


cos.simil <- matrix(NA, nrow=nrow(hsp.norm), ncol=nrow(hsp.norm))

for (i in 1:nrow(hsp.norm)){	
	for (j in 1:nrow(hsp.norm)){
		
		s1 <- as.numeric(hsp.norm[i,])
		s2 <- as.numeric(hsp.norm[j,])
		
		cos.simil[i,j] <- cosinesimilarity(s1,s2)	}	
}


cosine.animalia <- cos.simil[c2,c2]
dist.cos.animalia <- as.dist(1-cosine.animalia)
hc.cos <- hclust(dist.cos.animalia)
hc.cos$labels = hsp.name[c2]

cos.data <- dendro_data(hc.cos, type = "rectangle")
ggplot(segment(cos.data)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  labs(x = "", y="Cosine distance")



postscript("figures/Figure3/Fig3D.ps", width=4, height=4, paper="special", horizontal=T, onefile=F)

ggdendrogram(hc.cos, rotate = TRUE, size = 2)

dev.off()



hsp.composition <- rbind(hsp.norm[c(45, 46, 50, 53, 1, 23, 37), ], colMeans((hsp.norm)))
hsp.composition <- as.data.frame(hsp.composition)
hsp.composition$Name <- c(hsp.name[c(45,46,50,53,1,23,37)], "avrg")

hsp.composition$Name <- c("1-S.purpuratus", "2-C.gigas", "3-S.ratti","4-C.remanei", "5-H.glaber", "6-P.alecto", "7-N.furzeri", "8-avrg")  



hsp.compo <- melt(hsp.composition)

postscript("figures/Figure3/Fig3D2.ps", width=4, height=5, paper="special", horizontal=T, onefile=F)


ggplot(hsp.compo, aes(x=Name, y=value, fill=variable)) + 
	geom_bar(stat="identity") + 
	scale_fill_brewer(palette="Blues") + 
	theme_classic()

dev.off()
	
	





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hsp <- read.csv("data/chaperones/hsp.txt", sep='\t', stringsAsFactors = F)
agg <- read.table("data/chaperones/aggregation_summary.txt", sep='\t', header=1)
hyd <- read.table("data/chaperones/hydrophobicity_summary.txt", sep='\t', header=1)
diso <- read.table("data/chaperones/disorder_summary_complete.txt", sep='\t', header=1)

tax <- read.table("data/taxonomy/tree_taxonomy.txt", header=T, sep='\t')
c1 <- tax$kingdom == "Fungi"
c2 <- tax$kingdom == "Metazoa"
c3 <- tax$kingdom == "Viridiplantae"


data.lm <- hsp
data.lm$agg <- agg$Median
data.lm$hyd <- hyd$Median
data.lm$diso <- diso$Median



postscript("figures/Figure3/Fig3E.ps", width=6, height=3, paper="special", horizontal=T, onefile=F)

diso1 <- ggplot(data.lm[c2,], aes(x=Hsp40, y=diso)) + geom_point() +
	theme_classic() +
	theme(axis.text.y=element_text(size=12), axis.text.x = element_text(size=12)) +
	labs(x = "Hsp40", y="Protein disorder (%)")

diso2 <- ggplot(data.lm[c2,], aes(x=Hsp70, y=diso)) + geom_point() + 
	theme(axis.text.y=element_text(size=12), axis.text.x = element_text(size=12)) +
	labs(x = "Hsp70", y="Protein disorder (%)") +
	theme_classic()

plot_grid(diso1, diso2, labels ="", ncol = 2, align = 'h')

dev.off()




library(nnls)
A.hsp <- as.matrix(data.lm[,c(2:7)])
nnls.animalia.agg <- nnls(A.hsp[c2,], data.lm$agg[c2])
nnls.animalia.diso <- nnls(A.hsp[c2,], data.lm$diso[c2])


nnls_Rsq <- function(response, feature){
	
	fit <- nnls(as.matrix(response), as.matrix(feature) )
	coef.fit <- coef(fit)
	
	Rsq = rep(0, length(coef.fit) + 1)
	
	pred <- as.vector(coef.fit %*% t(response))
	lm.dummy <- lm(feature ~ pred)
	out <- display(lm.dummy)
	
	Rsq[1] <- out$r.squared
		
	for (i in 1:length(coef.fit)){		
		pred <- as.vector(coef.fit[i] %*% t(response[,i]))
		lm.dummy <- lm(feature ~ pred)
		out <- display(lm.dummy)
		Rsq[i+1] <- out$r.squared	
	}
	names <- c("7-All", "6-Hsp20", "5-Hsp40", "4-Hsp60", "3-Hsp70", "2-Hsp90", "1-Hsp100")
	cf <- c(0, coef.fit)
	result <- data.frame(Name=names, Coeff=cf, R2=Rsq)
	return(result)
}
	
	
nnls.agg <- nnls_Rsq(A.hsp, data.lm$agg)
nnls.diso <- nnls_Rsq(A.hsp, data.lm$diso)
nnls.hyd <- nnls_Rsq(A.hsp, data.lm$hyd)

nnls.agg.animalia <- nnls_Rsq(A.hsp[c2,], data.lm$agg[c2])
nnls.diso.animalia <- nnls_Rsq(A.hsp[c2,], data.lm$diso[c2])
nnls.hyd.animalia <- nnls_Rsq(A.hsp[c2,], data.lm$hyd[c2])

nnls.agg.fungi <- nnls_Rsq(A.hsp[c1,], data.lm$agg[c1])
nnls.diso.fungi <- nnls_Rsq(A.hsp[c1,], data.lm$diso[c1])
nnls.hyd.fungi <- nnls_Rsq(A.hsp[c1,], data.lm$hyd[c1])


nnls.agg.plantae <- nnls_Rsq(A.hsp[c3,], data.lm$agg[c3])
nnls.diso.plantae <- nnls_Rsq(A.hsp[c3,], data.lm$diso[c3])
nnls.hyd.plantae <- nnls_Rsq(A.hsp[c3,], data.lm$hyd[c3])






#plot.fungi.hyd <- ggplot(nnls.hyd.fungi, aes(x=Name, y=R2)) + geom_bar(stat="identity") +
	labs(y = "R^2", x="Hydrophobicity") +
	coord_flip()	
plot.fungi.agg <- ggplot(nnls.agg.fungi, aes(x=Name, y=R2)) + geom_bar(stat="identity") +
	labs(y = "R^2", x="Aggregation propensity") +
	coord_flip()
plot.fungi.diso <- ggplot(nnls.diso.fungi, aes(x=Name, y=R2)) + geom_bar(stat="identity") +
	labs(y = "R^2", x="Protein disorder (%)") +
	coord_flip()
	
	
#plot.animalia.hyd <- ggplot(nnls.hyd.animalia, aes(x=Name, y=R2)) + geom_bar(stat="identity") +
	labs(y = "R^2", x="Hydrophobicity") +
	coord_flip()	
plot.animalia.agg <- ggplot(nnls.agg.animalia, aes(x=Name, y=R2)) + geom_bar(stat="identity") +
	labs(y = "R^2", x="Aggregation propensity") +
	coord_flip()
plot.animalia.diso <- ggplot(nnls.diso.animalia, aes(x=Name, y=R2)) + geom_bar(stat="identity") +
	labs(y = "R^2", x="Protein disorder (%)") +
	coord_flip()
	
	
#plot.plantae.hyd <- ggplot(nnls.hyd.plantae, aes(x=Name, y=R2)) + geom_bar(stat="identity") +
	labs(y = "R^2", x="Hydrophobicity") +
	coord_flip()
plot.plantae.agg <- ggplot(nnls.agg.plantae, aes(x=Name, y=R2)) + geom_bar(stat="identity") +
	labs(y = "R^2", x="Aggregation propensity") +
	coord_flip()
plot.plantae.diso <- ggplot(nnls.diso.plantae, aes(x=Name, y=R2)) + geom_bar(stat="identity") +
	labs(y = "R^2", x="Protein disorder (%)") +
	coord_flip()


postscript("figures/Supplement/FigS3_R2.ps", width=7, height=8, paper="special", horizontal=T, onefile=F)

plot_grid(plot.fungi.agg, plot.fungi.diso, plot.animalia.agg, plot.animalia.diso, plot.plantae.agg, plot.plantae.diso, labels=c("Fungi", "Fungi", "Animalia", "Animalia", "Plantae", "Plantae"), ncol = 2, align = 'h')

dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hsp <- read.csv("data/chaperones/hsp.txt", sep='\t', stringsAsFactors = F)
agg <- read.table("data/chaperones/aggregation_summary.txt", sep='\t', header=1)
hyd <- read.table("data/chaperones/hydrophobicity_summary.txt", sep='\t', header=1)
diso <- read.table("data/chaperones/disorder_summary_complete.txt", sep='\t', header=1)

tax <- read.table("data/taxonomy/tree_taxonomy.txt", header=T, sep='\t')
c1 <- tax$kingdom == "Fungi"
c2 <- tax$kingdom == "Metazoa"
c3 <- tax$kingdom == "Viridiplantae"


data.lm <- hsp
data.lm$agg <- agg$Median
data.lm$hyd <- hyd$Median
data.lm$diso <- diso$Median


data.lm$tax = rep("NA", nrow(agg))
data.lm$tax[c1] <- "1 - Fungi"
data.lm$tax[c2] <- "2 - Metazoa"
data.lm$tax[c3] <- "3 - Plantae"

data.lm <- data.lm[data.lm$tax!="NA",]

agg2  <- melt( data.frame(Name=data.lm$Name, Median=data.lm$agg, Tax=data.lm$tax) )
hyd2  <- melt( data.frame(Name=data.lm$Name, Median=data.lm$hyd, Tax=data.lm$tax) )
diso2 <- melt( data.frame(Name=data.lm$Name, Median=data.lm$diso, Tax=data.lm$tax) )



postscript("figures/Supplement/FigS3B.ps", width=5, height=10, paper="special", horizontal=T, onefile=F)

p1 <- ggplot(agg2, aes(x=value, fill=Tax) ) + 
	geom_density() + 
	theme_classic() + 
	theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) + 
	labs(x = "Median agg.score ", y="Density")	
dens.agg <- p1 + facet_grid(rows=vars(Tax) ) 


p2 <- ggplot(hyd2, aes(x=value, fill=Tax) ) + 
	geom_density() + 
	theme_classic() + 
	theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) + 
	labs(x = "Median hydrophobicity ", y="Density")	
dens.hyd <- p2 + facet_grid(rows=vars(Tax) ) 


p3 <- ggplot(diso2, aes(x=value, fill=Tax) ) + 
	geom_density() + 
	theme_classic() + 
	theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) + 
	labs(x = "Median % disorder ", y="Density")
dens.diso <- p3 + facet_grid(rows=vars(Tax) ) 


plot_grid(dens.hyd, dens.agg, dens.diso, labels="", ncol = 1, align = 'v')


dev.off()



