# Load data and library
# generate panels for Figure 2

library(RColorBrewer)
library(gplots)


#-----------------

# TREE FOR HEATMAP

library(ape)
library(ggtree)

tree <- read.tree("../../data/Phylobayes/eucaryote_from_concensus_5_11.newick")
names <- read.table("../../data/list_organisms_final_tree_short.txt")
tax <- read.table("../../data/taxonomy/tree_taxonomy_corrected.txt", header=T, sep='\t')

tree$tip.label <- as.character(names[,1])

postscript("../../figures/Figure2/Fig2_tree.ps", width=5, height=10, paper="special", horizontal=T, onefile=F)

plot(tree, use.edge.length = FALSE, main = "")

dev.off()



#--------------------

library(reshape)

hsp <- read.csv("../../data/Heatmap_HSP/HSP_phylogeny_5_november.csv")
names <- read.table("../../data/list_organisms_final_tree_short.txt")
chaperones <- c("Hsp20","Hsp40", "Hsp60", "Hsp70", "Hsp90", "Hsp100")

data <- as.matrix( hsp[,c(3:8)])
data <- data.frame(cbind(names,data))
colnames(data) <- c("species", chaperones)
data$species <- as.character(data$species)
data$species <- factor(data$species, levels=unique(data$species))


data.lp <- melt(data)


# PLOT OF THE LINES FOR THE HEATMAP


postscript("../../figures/Figure2/Fig2_lines.ps", width=6, height=6, paper="special", horizontal=T, onefile=F)

p <- ggplot(data.lp, aes(x=species, y=value) ) + geom_line(aes(x=as.numeric(species), y=value)) + theme_classic() + geom_hline(aes(yintercept = mean(value)) )
p + facet_grid(rows = vars(variable),scales="free_y") 

dev.off()


#reverse order for heatmap


data.hm <- as.matrix( hsp[,c(3:8)])
data.hm <- data.frame(cbind(names,data.hm))
data.hm <- data.hm[162:1,]
colnames(data.hm) <- c("species", chaperones)
rownames(data.hm) <- c(1:162)

data.hm$species <- as.character(data.hm$species)
data.hm$species <- factor(data.hm$species, levels=unique(data.hm$species))
data.hm <- melt(data.hm)
data.hm$breaks <- cut(as.numeric(data.hm$value), c(0,1,2,3,5,10,15,20,25,30,40,50,60,75,90,200), right=F)


# PLOT OF THE HEATMAP

postscript("../../figures/Figure2/Fig2_heatmap.ps", width=6, height=6, paper="special", horizontal=T, onefile=F)

ggplot(data.hm, aes(x = variable, y = species)) + geom_tile(aes(fill=breaks),color="white", size=0.1) + labs(x=NULL, y=NULL) + theme(axis.text.y = element_blank(), axis.ticks=element_blank() ) + scale_fill_manual(values=c("#29ABE2", "#3E99CA", "#5488B3", "#6A769C", "#96536E", "#C13040", "#ED0D12", "#FF1A00", "#FF4600", "#FF7200", "#FF9D00", "#FFB300", "#FFC900", "#FFD700", 'yellow'))
      
dev.off()




#---------------------------------

# PLOT OF CHAPERONE DENSITIES BY CATEGORY

library(ggplot2)

hsp <- read.csv("../../data/Heatmap_HSP/HSP_phylogeny_5_november.csv")
names <- read.table("../../data/list_organisms_final_tree_short.txt")
chaperones <- c("Hsp20","Hsp40", "Hsp60", "Hsp70", "Hsp90", "Hsp100")
tax <- read.table("../../data/taxonomy/tree_taxonomy_corrected.txt", header=T, sep='\t')


data <- as.matrix( hsp[,c(3:8)])
data <- data.frame(cbind(names,data))
colnames(data) <- c("species", chaperones)
data$species <- as.character(data$species)
data$species <- factor(data$species, levels=unique(data$species))

data$tax <- rep("NA", nrow(data))
data$tax[tax$phylum == "Ascomycota"] <-"1 - FungiA"
data$tax[tax$phylum == "Basidiomycota"] <-"2 - FungiB"
data$tax[tax$kingdom == "Protista"] <- "3 - Protista"
data$tax[tax$kingdom == "Animalia"] <- "4 - Animalia"
data$tax[tax$kingdom == "Plantae"] <- "5 - Plantae"

data <- data[data$tax!="NA",]

data.grid <- melt(data)


postscript("../../figures/Supplement/FigS2.ps", width=9, height=6, paper="special", horizontal=T, onefile=F)

p <- ggplot(data.grid, aes(x=value, fill=variable) ) + geom_density() + theme_classic()
p + facet_grid(vars(tax),vars(variable), scales = 'free') 

dev.off()

