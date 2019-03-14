# Panels for Figure 4 and Figure S4

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

agg <- read.table("data/chaperones/aggregation_summary.txt", sep='\t', header=1)

postscript("figures/Figure4/Fig4A.ps", width=8, height=3, paper="special", horizontal=T, onefile=F)


ggplot(agg, aes(x=treepos, y=Median)) + 
	geom_hline(aes(yintercept = mean(Median)), size=0.5, linetype="twodash" ) + 
	geom_errorbar(aes(ymin=Q1, ymax=Q3), size = 0.1, width=0) +
	geom_point(color="red", shape=16, show.legend = FALSE) + 
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(
size=20)) + 
	labs(y = "AP [a.u.]", x="") 
		
dev.off()
	



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


tax <- read.table("data/taxonomy/tree_taxonomy.txt", header=T, sep='\t')

c1 <- tax$kingdom == "Fungi"
c2 <- tax$kingdom == "Metazoa"
c3 <- tax$kingdom == "Viridiplantae"

agg$tax = rep("NA", nrow(agg))
agg$tax[c1] <- "1 - Fungi"
agg$tax[c2] <- "2 - Metazoa"
agg$tax[c3] <- "3 - Plantae"

agg <- agg[agg$tax!="NA",]

agg2 <- melt( data.frame(Name=agg$Name, Median=agg$Median, Tax=agg$tax) )

postscript("figures/Figure4/Fig4B.ps", width=5, height=4, paper="special", horizontal=T, onefile=F)

p <- ggplot(agg2, aes(x=value, fill=Tax) ) + 
	geom_density() + 
	theme_classic() + 
	theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) + 
	labs(x = "Median agg.score ", y="Density")
	
p + facet_grid(rows=vars(Tax) ) 

dev.off()





#--------------------------------------------

postscript("figures/Figure4/Fig4C.ps", width=2, height=4, paper="special", horizontal=T, onefile=F)

ggplot(subset(agg2, Tax == "2 - Metazoa"), aes(x=1, y=value) ) + 	
	geom_violin(aes(fill="#bbbbbb"), draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) + 
	theme_classic( ) +
	scale_fill_manual(values = c("grey")) +
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=16)) +
	labs(x = "Animalia", y="Median aggregation score") +

	geom_point(aes(x=1, y=agg$Median[agg$Name=="Crassostrea_gigas"]), col="#025074", shape=16, size=6) + 			#oyster
	geom_point(aes(x=1, y=agg$Median[agg$Name=="Heterocephalus_glaber"]), col="#b54a03", shape=15, size=6) +			#nmrat
	geom_point(aes(x=1, y=agg$Median[agg$Name=="Nothobranchius_furzeri"]), col="#03b582", shape=18, size=6) +			#killi
	geom_point(aes(x=1, y=agg$Median[agg$Name=="Strongylocentrotus_purpuratus"]), col="#991e7b", shape=17, size=6) 			#urchin

dev.off()






#--------------------------------------------


agg <- as.data.frame(read.table("code/aggregation_summary.txt", sep='\t', header=1) )

hsp <- read.csv("data/chaperones/hsp.txt", stringsAsFactors = F, sep='\t')
names <- read.table("data/tree/tree_list_nodes.txt")
chaperones <- c("Hsp20","Hsp40", "Hsp60", "Hsp70", "Hsp90", "Hsp100")
colnames(hsp) <- c("species", chaperones)
hsp <- as.data.frame(hsp)

data <- data.frame(Name=as.character(agg$Name), Hsp20=hsp$Hsp20, Hsp40=hsp$Hsp40, 
                   Hsp60=hsp$Hsp60, Hsp70=hsp$Hsp70, Hsp90=hsp$Hsp90, Hsp100=hsp$Hsp100, 
                   Mean=agg$Mean, Sum=agg$Sum, Nprot=agg$Nprot)

tax <- read.table("data/taxonomy/tree_taxonomy.txt", header=T, sep='\t')
sel_animalia <- tax$kingdom == "Metazoa"

data.animalia <- data.frame(Name=data$Name[sel_animalia], Hsp=rowSums(data[sel_animalia,c(3,5,6)]), Sum=data$Sum[sel_animalia] )


postscript("figures/Figure4/Fig4D.ps", width=6, height=5, paper="special", horizontal=T, onefile=F)

ggplot(data.animalia, aes(x=Hsp, y=Sum)) + geom_point(size=2) + 
	labs(x="Core chaperone network", y="Proteome aggregation score [a.u.]") +
	scale_y_continuous(labels = scales::comma) + 
	#geom_smooth(method='lm') +
	#geom_density_2d(  )
	geom_point(aes(x=61, y=data.animalia$Sum[data.animalia$Name=="Nothobranchius_furzeri"]), col="#03b582", shape=18, size=6) + 		#killi
	geom_point(aes(x=111, y=data.animalia$Sum[data.animalia$Name=="Heterocephalus_glaber"]), col="#b54a03", shape=15, size=6) + 		#nmrat
	geom_point(aes(x=127, y=data.animalia$Sum[data.animalia$Name=="Crassostrea_gigas"]), col="#025074", shape=16, size=6) + 		#oyster
	geom_point(aes(x=78, y=data.animalia$Sum[data.animalia$Name=="Strongylocentrotus_purpuratus"]), col="#991e7b", shape=17, size=6) 			#urchin


dev.off()



#--------------------------------------------


data.animalia.2 <- data.frame(Name=data$Name[sel.animalia], Hsp=rowSums(data[sel.animalia,c(2:7)]), Sum=data$Sum[sel.animalia] )


postscript("figures/Supplement/FigS4A.ps", width=6, height=5, paper="special", horizontal=T, onefile=F)

nmrat <- data.animalia.2[data.animalia.2$Name=="Heterocephalus_glaber",]
killi <- data.animalia.2[data.animalia.2$Name=="Nothobranchius_furzeri",]
oyster <- data.animalia.2[data.animalia.2$Name=="Crassostrea_gigas",]
urchin <- data.animalia.2[data.animalia.2$Name=="Strongylocentrotus_purpuratus",]

ggplot(data.animalia.2, aes(x=Hsp, y=Sum)) + geom_point(size=2) + 
	labs(x="Complete chaperone network", y="Proteome aggregation score [a.u.]") +
	scale_y_continuous(labels = scales::comma) + 
	#geom_smooth(method='lm') +
	#geom_density_2d(  )
	geom_point(aes(x=killi$Hsp, y=killi$Sum), col="#03b582", shape=18, size=6) + 		#killi
	geom_point(aes(x=nmrat$Hsp, y=nmrat$Sum), col="#b54a03", shape=15, size=6) + 		#nmrat
	geom_point(aes(x=oyster$Hsp, y=oyster$Sum), col="#025074", shape=16, size=6) + 		#oyster
	geom_point(aes(x=urchin$Hsp, y=urchin$Sum), col="#991e7b", shape=17, size=6) 			#urchin


dev.off()





#--------------------------------------------

# all by all suppl


data.2 <- data.frame(Name=data$Name, Hsp=rowSums(data[,c(2:7)]), Sum=data$Sum)

tax <- read.table("data/taxonomy/tree_taxonomy.txt", header=T, sep='\t')

c1 <- tax$kingdom == "Fungi"
c2 <- tax$kingdom == "none"
c3 <- tax$kingdom == "Metazoa"
c4 <- tax$kingdom == "Viridiplantae"

data.2$tax = rep("NA", nrow(agg))
data.2$tax[c1] <- "1- Fungi"
data.2$tax[c2] <- "2- Protista"
data.2$tax[c3] <- "3- Animalia"
data.2$tax[c4] <- "4- Plantae"

data.2 <- data.2[data.2$tax!="NA",]


postscript("figures/Supplement/FigS4B.ps", width=7, height=5, paper="special", horizontal=T, onefile=F)

ggplot(data.2, aes(x=Hsp, y=Sum)) + geom_point(aes(color=tax), size=2) + 
	labs(x="Complete chaperone network", y="Proteome aggregation score [a.u.]") +
	scale_y_continuous(labels = scales::comma) 

dev.off()