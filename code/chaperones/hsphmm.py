import os, sys
import numpy as np
import pandas as pd
import subprocess
import glob



def read_taxid():
    """
    read taxid into dict
    """

    taxid = {} 

    with open('taxid.txt', 'r') as f:
        for line in f:
            current_line = line.split()
            current_id = current_line[0]
            current_name = current_line[1] 
            taxid[current_id] = current_name

    return taxid


taxid = read_taxid()



hsp20='hmmer/hsp20.hmm'
hsp40='hmmer/hsp40.hmm'
hsp60='hmmer/hsp60.hmm'
hsp70='hmmer/hsp70.hmm'
hsp90='hmmer/hsp90.hmm'
hsp100='hmmer/hsp100.hmm'




 
def get_counts(hmm, proteome):
    """
    get counts of chaperone type in proteome
    """

    output = subprocess.check_output(['hmmsearch', '-E', '1e-5', '--cpu', '4', hmm, proteome], stderr=subprocess.STDOUT).decode('UTF-8')
    hits = [ line for line in output.split('\n') if '>>' in line]

    N = len(hits)

    return N



def proteome_counts(proteome):
    """
    get chaperone counts for proteome
    """

    N_hsp20 = get_counts(hsp20, proteome)
    N_hsp40 = get_counts(hsp40, proteome)
    N_hsp60 = get_counts(hsp60, proteome)
    N_hsp70 = get_counts(hsp70, proteome)
    N_hsp90 = get_counts(hsp90, proteome)
    N_hsp100 = get_counts(hsp100, proteome)

    return N_hsp20, N_hsp40, N_hsp60, N_hsp70, N_hsp90, N_hsp100







def run_all(proteomes):
    """
    analyze all proteomes
    """

    chaperone_df = pd.DataFrame(columns=['Tax', 'Name', 'Hsp20', 'Hsp40', 'Hsp60','Hsp70', 'Hsp90', 'Hsp100']) 
    counter = 0

    for i in proteomes:
        name = i.split('/')[-1].split('.')[0]
        ID = name.split('_')[1] 
        tax = taxid.get(ID, 'none')
        print(name, ID, tax)
        count_hsp20, count_hsp40, count_hsp60, count_hsp70, count_hsp90, count_hsp100 = proteome_counts(i)

        chaperone_df.loc[counter] = [ID, tax, count_hsp20, count_hsp40, count_hsp60, count_hsp70, count_hsp90, count_hsp100] 
        counter += 1

        if counter % 100 == 0:
            print(counter, "proteomes analyzed")

    chaperone_df.to_csv("hsp.txt", sep='\t', header=True, index=False)

    return chaperone_df


if __name__ == '__main__':

    list_proteomes = glob.glob('/home/sebastian/Denali/Yasmine/uniprot/*.fasta')

    df_tree = pd.read_csv("../../data/tree/tree_uniprot.txt", sep='\t', header=0)
    
    chaperone_df = pd.DataFrame(columns=['Name', 'Hsp20', 'Hsp40', 'Hsp60','Hsp70', 'Hsp90', 'Hsp100']) 
    counter = 0

    for i in range( len(df_tree) ):
        current_name = df_tree.iloc[i]['name']
        current_unip = df_tree.iloc[i]['uniprot']

        if "Notho" in current_unip :
            current_proteome = 'Nothobranchius_furzeri.fa'
        else:
            current_proteome = glob.glob("/home/sebastian/Denali/Yasmine/uniprot/"+current_unip+"*.fasta")[0]

        count_hsp20, count_hsp40, count_hsp60, count_hsp70, count_hsp90, count_hsp100 = proteome_counts(current_proteome)
        chaperone_df.loc[counter] = [current_name, count_hsp20, count_hsp40, count_hsp60, count_hsp70, count_hsp90, count_hsp100] 
        counter += 1
    
    chaperone_df.to_csv("hsp.txt", sep='\t', header=True, index=False)


#    run_all(list_proteomes)

