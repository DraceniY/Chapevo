import sys, os
import pandas as pd
import numpy as np
from ete3 import Tree
import ete3



if __name__ == "__main__":

    # use ete3's get_distance function to compute pairwise additive distances between leaves in tree

    tree = Tree("../../data/tree/tree.nw")
    list_nodes = list(tree.get_leaves() )

    print(len(list_nodes))


    list_names = []

    dmat = np.zeros(( len(list_nodes), len(list_nodes) ))

    for i in range( len(list_nodes) ):
        list_names.append( list_nodes[i].name )
        for j in range(i, len(list_nodes) ):

            d = tree.get_distance(list_nodes[i], list_nodes[j], topology_only=False)
            dmat[i,j] = dmat[j,i] = round(d, 5)

    dist_df = pd.DataFrame(data=dmat, index=list_names, columns=list_names )
    dist_df.to_csv("../../data/tree/tree_distancematrix.txt", sep='\t', header=True)
 
