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
	geom_point(aes(x=1, y=agg$Median[agg$Name=="Strongylocentrotus_purpuratus"]), col="#991e7b", shape=17, size=6) 	#urchin
	
dev.off()




#--------------------------------------------


agg <- as.data.frame(read.table("data/chaperones/aggregation_summary.txt", sep='\t', header=1) )

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

fit.rq.animalia <- rq(Sum ~ Hsp, data=data.animalia)


postscript("figures/Figure4/Fig4D_revised.ps", width=6, height=5, paper="special", horizontal=T, onefile=F)

ggplot(data.animalia, aes(x=Hsp, y=Sum)) + 
	geom_abline(slope= fit.rq.animalia$coeff[2], intercept = fit.rq.animalia$coeff[1], color="orange", size=7) +
	geom_point(size=2) + 
	labs(x="Core chaperone network", y="Proteome aggregation score [a.u.]") +
	scale_y_continuous(labels = scales::comma) + 
	theme_classic( ) +
	#geom_smooth(method='lm') +
	#geom_density_2d(  )
	geom_point(aes(x=61, y=data.animalia$Sum[data.animalia$Name=="Nothobranchius_furzeri"]), col="#03b582", shape=18, size=6) + 		#killi
	geom_point(aes(x=111, y=data.animalia$Sum[data.animalia$Name=="Heterocephalus_glaber"]), col="#b54a03", shape=15, size=6) + 		#nmrat
	geom_point(aes(x=127, y=data.animalia$Sum[data.animalia$Name=="Crassostrea_gigas"]), col="#025074", shape=16, size=6) + 		#oyster
	geom_point(aes(x=78, y=data.animalia$Sum[data.animalia$Name=="Strongylocentrotus_purpuratus"]), col="#991e7b", shape=17, size=6)		#urchin
	

dev.off()



#--------------------------------------------


age <- read.table("data/chaperones/summary_age.txt", header=T)

age.animalia <- age[sel_animalia,]
age.animalia$residuals <- round(fit.rq.animalia$residuals / 1000000, 5)

fit.rq.age <- rq(Age ~ residuals, data=age.animalia)

age.nmrat <- age.animalia[age.animalia$Name=="Heterocephalus_glaber",]
age.killi <- age.animalia[age.animalia$Name=="Nothobranchius_furzeri",]
age.oyster <- age.animalia[age.animalia$Name=="Crassostrea_gigas",]
age.urchin <- age.animalia[age.animalia$Name=="Strongylocentrotus_purpuratus",]

postscript("figures/Figure4/Fig4E.ps", width=5, height=4, paper="special", horizontal=T, onefile=F)

ggplot(age.animalia, aes(x=Age, y=residuals, size=Age)) +
	#geom_abline(slope=0, intercept=0, size=0.3, color="#555555") +
	geom_point() + 
	theme_classic() + 
	labs(x="Age", y="Residuals [a.u.]") +
	#geom_abline(slope= fit.rq.age$coeff[2], intercept = fit.rq.age$coeff[1], color="orange", size=0.5) +	
	theme(axis.text.x=element_text(size=16), axis.text.y = element_text(size=16)) +
	geom_point(aes(x=age.killi$Age, y=age.killi$residuals), col="#03b582", shape=18, size=3) + 	
	geom_point(aes(x=age.nmrat$Age, y=age.nmrat$residuals), col="#b54a03", shape=15, size=4) + 	
	#geom_point(aes(x=age.oyster$Age, y=age.oyster$residuals), col="#025074", shape=16, size=2) + 	
	geom_point(aes(x=age.urchin$Age, y=age.urchin$residuals), col="#991e7b", shape=17, size=5) 

dev.off()






postscript("figures/Supplement/FigS4_age.ps", width=8, height=4, paper="special", horizontal=T, onefile=F)

p1 <- ggplot(age.animalia, aes(x=Age/Weight, y=residuals, size=Weight)) +
	geom_abline(slope=0, intercept=0, size=0.3, color="#AAAAAA") +
	geom_point() + 
	theme_classic() + 
	labs(x="Age/Weight", y="Residuals [a.u.]") + 
	theme(legend.position="none")

p2 <- ggplot(age.animalia, aes(x=log(Age/Weight), y=residuals, size=Weight)) +
	geom_abline(slope=0, intercept=0, size=0.3, color="#AAAAAA") +
	geom_point() + 
	theme_classic() + 
	labs(x="log(Age/Weight)", y="Residuals [a.u.]") + 
    theme(legend.position="none")

plot_grid(p1, p2, labels="", ncol=2, align="h")

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


postscript("figures/Figure4/Fig4B_new.ps", width=7, height=5, paper="special", horizontal=T, onefile=F)

ggplot(data.2, aes(x=Hsp, y=Sum/1000000)) + geom_point(aes(color=tax), size=2) + 
	labs(x="Size of chaperone network", y="Proteome aggregation score [a.u.]") +
	theme_classic() +
	scale_y_continuous(labels = scales::comma) +
	geom_smooth(aes(col=tax), method="rlm", se=F)

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~

data.3 <- data.frame(Name=data$Name, Hsp=rowSums(data[,c(2:7)]), Nprot=data$Nprot)

tax <- read.table("data/taxonomy/tree_taxonomy.txt", header=T, sep='\t')

c1 <- tax$kingdom == "Fungi"
c2 <- tax$kingdom == "none"
c3 <- tax$kingdom == "Metazoa"
c4 <- tax$kingdom == "Viridiplantae"

data.3$tax = rep("NA", nrow(agg))
data.3$tax[c1] <- "1- Fungi"
data.3$tax[c2] <- "2- Protista"
data.3$tax[c3] <- "3- Animalia"
data.3$tax[c4] <- "4- Plantae"

data.3 <- data.3[data.3$tax!="NA",]


postscript("figures/Supplement/FigS4B_new.ps", width=7, height=5, paper="special", horizontal=T, onefile=F)

ggplot(data.3, aes(x=Hsp, y=Nprot)) + geom_point(aes(color=tax), size=2) + 
	labs(x="Size of chaperone network", y="Size of proteome") +
	theme_classic() +
	scale_y_continuous(labels = scales::comma) +
	geom_smooth(aes(col=tax), method="rlm", se=F)

dev.off()
	

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## coefficients (slope, R2) for table rq + rlm for both aggregation score + size

get_coeffs <- function(feature, response){
	
	fit.rq <- rq(feature ~ response)
	slope.rq <- fit.rq$coeff[2]
	R2.rq <- cor(feature, fit.rq$fitted.values)^2
	
	fit.rlm <- rlm(feature ~ response)
	slope.rlm <- fit.rlm$coeff[2]
	R2.rlm <- cor(feature, fit.rlm$fitted.values)^2
	
	output <- cbind( c(slope.rq, slope.rlm), c(R2.rq, R2.rlm))
	row.names(output) <- c("rq", "rlm")
	colnames(output) <- c("slope", "R2")
	
	output
	
}

rel_agg.fungi <- get_coeffs(data.2$Sum[c1], data.2$Hsp[c1])
rel_agg.protists <- get_coeffs(data.2$Sum[c2], data.2$Hsp[c2])
rel_agg.animalia <- get_coeffs(data.2$Sum[c3], data.2$Hsp[c3])
rel_agg.plantae <- get_coeffs(data.2$Sum[c4], data.2$Hsp[c4])

rel_size.fungi <- get_coeffs(data.3$Nprot[c1], data.3$Hsp[c1])
rel_size.protists <- get_coeffs(data.3$Nprot[c2], data.3$Hsp[c2])
rel_size.animalia <- get_coeffs(data.3$Nprot[c3], data.3$Hsp[c3])
rel_size.plantae <- get_coeffs(data.3$Nprot[c4], data.3$Hsp[c4])


fit.animalia.agg.rq <- rq(data.2$Sum[c3] ~ data.2$Hsp[c3])
fit.animalia.size.rq <- rq(data.3$Nprot[c3] ~ data.3$Hsp[c3])




#~~~ coeffs on correlations

heavies <- c(5, 14, 34, 35)

#proteome aggregation score
cor(age.animalia$Age, age.animalia$residuals, use="complete")
cor.test(age.animalia$Age, age.animalia$residuals, use="complete")

cor(age.animalia$Age[-heavies], age.animalia$residuals[-heavies], use="complete")
cor.test(age.animalia$Age[-heavies], age.animalia$residuals[-heavies], use="complete")


# proteome size

data.animalia.nprot <- data.frame(Name=data$Name, Hsp=rowSums(data[,c(2:7)]), Nprot=data$Nprot)
tax <- read.table("data/taxonomy/tree_taxonomy.txt", header=T, sep='\t')
c3 <- tax$kingdom == "Metazoa"
data.animalia.nprot$tax = rep("NA", nrow(agg))
data.animalia.nprot$tax[c3] <- "3- Animalia"
data.animalia.nprot <- data.animalia.nprot[c3,]

fit.rq.protsize <- rq(Nprot ~ Hsp, data=data.animalia.nprot)

sel_animalia <- tax$kingdom == "Metazoa"
age <- read.table("data/chaperones/summary_age.txt", header=T)
age.animalia.nprot <- age[sel_animalia,]
age.animalia.nprot$residuals <- round(fit.rq.protsize$residuals / 1000000, 5)

fit.rq.age.nprot <- rq(Age ~ residuals, data=age.animalia.nprot)

cor(age.animalia.nprot$Age, age.animalia.nprot$residuals, use="complete")
cor.test(age.animalia.nprot$Age, age.animalia.nprot$residuals, use="complete")

cor(age.animalia.nprot$Age[-heavies], age.animalia.nprot$residuals[-heavies], use="complete")
cor.test(age.animalia.nprot$Age[-heavies], age.animalia.nprot$residuals[-heavies], use="complete")

