import sys, os
import pandas as pd
import numpy as np
from ete3 import Tree


accession = pd.read_csv("../18SrRNA/refseq_download.txt", sep='\t', header=0, index_col=0) 


def load_taxdict():
    """
    load mapping from tax id to species name 
    """
    tax = {}
    with open("../../data/taxonomy/tree_taxid.txt", 'r') as file:
        for line in file:
            current_line = line.split() 
            current_taxid = current_line[0]
            current_name = current_line[1]
            tax[current_taxid] = current_name 

    return tax


def get_uniprot(tax, acc):
    """
    get matching uniprot id given taxid 
    """

    taxid1 = list(acc['taxid'])
    taxid2 = list(acc['taxid2'])
    taxid3 = list(acc['species_taxid'])

    if int(tax) in taxid1:
        unip = acc[acc['taxid']==int(tax)]['uniprot'].item()

    elif int(tax) in taxid2:
        unip = acc[acc['taxid2']==int(tax)]['uniprot'].item()

    elif int(tax) in taxid3:
        unip = acc[acc['species_taxid']==int(tax)]['uniprot'].item()
    
    else:
        unip = 'NA'

    return unip


def update_tip_names(tree, taxdict):
    """
    convert taxonomy id to species name in tip lables
    """

    list_nodes = []
    uniprot_mapping = pd.DataFrame(columns=['taxid', 'name', 'uniprot'])

    counter = 0
    for node in tree.traverse("postorder"):
        current_name = node.name

        if 'NMR' in current_name:
            new_name = "Heterocephalus_glaber"
            node.name = new_name
            list_nodes.append(node.name)
            taxid = "NA" 
            uniprot_mapping.loc[counter] = (taxid, new_name, "UP000006813")
            counter += 1

        elif 'Nfurzer' in current_name:
            new_name = "Nothobranchius_furzeri"
            node.name = new_name
            list_nodes.append(node.name)
            taxid = "NA"
            uniprot_mapping.loc[counter] = (taxid, new_name, new_name)
            counter += 1

        elif 'TAX' in current_name:
            taxid = current_name[3:].split('x')[0]
            new_name = taxdict.get(taxid, taxid) 
            node.name = new_name 
            list_nodes.append(node.name)
            unip = get_uniprot(taxid, accession)
            uniprot_mapping.loc[counter] = (taxid, new_name, unip)
            counter += 1


 
    tree.write(outfile="../../data/tree/tree.nw")

    nodes_df = pd.DataFrame(list_nodes)
    nodes_df.to_csv("../../data/tree/tree_list_nodes.txt", index=False, header=False)

    uniprot_mapping.to_csv("../../data/tree/tree_uniprot.txt", sep='\t', index=False, header=True)

    return tree, list_nodes



if __name__ == '__main__':

    tax = load_taxdict()
    tree = Tree("../../data/tree/arbre_15_10000.con.tre")
    new_tree, nodes = update_tip_names(tree, tax)

