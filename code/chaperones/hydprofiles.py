import sys, os
import numpy as np
import pandas as pd
import pickle
import glob
from Bio import SeqIO


# Kyte-Doolittle hydrophobicity scale, normalized to have zero mean and unitary standard deviation (Pechmann et al. PNAS 2009)
aa_hydrophobicity = {"I": 4.5,  "V": 4.2,  "L": 3.8,  "F": 2.8,  "C": 2.5,  "M": 1.9,  "A": 1.8,  "G": -0.4, "T": -0.7, "S": -0.8, 
                     "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5}


def compute_hydprofile(proteome):
    """
    loop through fasta and get pkl of hydrophobicity profiles
    """

    seqs = SeqIO.to_dict(SeqIO.parse(proteome, 'fasta'))

    num_seqs = len(seqs)

    scores = np.zeros(( num_seqs ))

    for dx, d in enumerate( seqs.keys() ):
        current_seq = seqs[d].seq

        current_hyd = np.zeros(( len(current_seq) ))

        for ax, aa in enumerate(current_seq):
            current_hyd[ax] = aa_hydrophobicity.get(aa, np.nan)

        scores[dx] = np.around( np.nanmean( current_hyd ), 3)

    scores = scores[scores != np.nan]

    current_mean = np.mean(scores)
    current_median = np.median(scores)
    current_sd   = np.std(scores)
    current_q1   = np.percentile(scores, 25)
    current_q3   = np.percentile(scores, 75)
    current_min  = np.min(scores)
    current_max  = np.max(scores)
    current_sum  = np.sum(scores)

    return num_seqs, round(current_median,2), round(current_mean,2), round(current_sd,2), round(current_q1,2), round(current_q3,2), round(current_min,2), round(current_max,2), round(current_sum, 2)





if __name__ == '__main__':

    # set to directory that contains all unzipped uniprot proteome fasta files
    UNIPROT_DIR = '/home/sebastian/Denali/Yasmine/uniprot/'

    list_proteomes = glob.glob(UNIPROT_DIR + '*.fasta')
    list_organisms = [line.strip().decode("utf-8") for line in open('../../data/tree/tree_list_nodes.txt', 'rb')]


    df_tree = pd.read_csv("../../data/tree/tree_uniprot.txt", sep='\t', header=0)

    hyd_stats = pd.DataFrame(columns=['Name', 'treepos', 'Nprot', 'Median', 'Mean', 'Std', 'Q1', 'Q3', 'Min', 'Max', 'Sum'])

    counter = 0
    for i in range( len(df_tree) ):
        current_name = df_tree.iloc[i]['name']
        current_unip = df_tree.iloc[i]['uniprot']

        if "Notho" in current_unip :
            current_proteome = 'Nothobranchius_furzeri.fa'
        else:
            current_proteome = glob.glob(UNIPROT_DIR + current_unip+"*.fasta")[0]

        if current_name in list_organisms:
            treepos = int(list_organisms.index(current_name)) + 1
        else:
            treepos = "NA"

        L, Md, Mn, St, Q1, Q3, Min, Max, Sum = compute_hydprofile(current_proteome)
        hyd_stats.loc[counter] = (current_name, treepos, L, Md, Mn, St, Q1, Q3, Min, Max, Sum)

        print(current_name, treepos)

        counter += 1

    hyd_stats.sort_values('treepos', inplace=True)
    hyd_stats.to_csv("hydrophobicity_summary.txt", sep='\t', header=True, index=False)






    
