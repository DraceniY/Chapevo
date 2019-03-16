# get taxonomy info
library(myTAI)


list_species <- scan("data/tree/tree_list_nodes.txt", what="", sep="\n")



getinfo <- function(name){
  # function that retrieves taxonomy classification	
	t <- taxonomy(name)
	print(t)
	
	if (length(which(t$rank == "superkingdom") > 0)){superkingdom <- t$name[which(t$rank=="superkingdom")]}
	else {superkingdom <- "none"}
	
	if (length(which(t$rank == "kingdom") > 0)){kingdom <- t$name[which(t$rank=="kingdom")]}
	else {kingdom <- "none"}
	
	if (length(which(t$rank == "phylum") > 0)){phylum <- t$name[which(t$rank=="phylum")]}
	else {phylum <- "none"}
	
	if (length(which(t$rank == "subphylum") > 0)){subphylum <- t$name[which(t$rank=="subphylum")]}
	else {subphylum <- "none"}
	
	if (length(which(t$rank == "family") > 0)){family <- t$name[which(t$rank=="family")]}
	else {family <- "none"}
	
	
	
	classes <- c(superkingdom, kingdom, phylum, subphylum, family)
	output <- c(name, classes)
	output
}


#initialize DF
L = length(list_species)
tax_df <-data.frame(name=character(L), superkingdom=character(L), kingdom=character(L), phylum=character(L), subphylum=character(L), family=character(L), stringsAsFactors = FALSE)

#loop through species
for (i in 1:length(list_species)){
  out <- getinfo(list_species[i])
  tax_df[i,] = out
}

rownames(tax_df) <- 1:nrow(tax_df)
write.table(tax_df,"../../../data/taxonomy/tree_taxonomy.txt",sep="\t",row.names=FALSE, quote=F)
