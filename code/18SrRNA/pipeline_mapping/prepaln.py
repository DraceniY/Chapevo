import sys, os
import subprocess
import numpy as np
import pandas as pd
from Bio.PDB import *
from Bio import SeqIO
from Bio import AlignIO
from Bio import Align
import itertools as it




def filterandparse_sequences(fastaOUT, theta):
    """
    filter for gaps and N characters
    """

    data = pd.read_csv("../../../data/18SrRNA/pipeline_mapping/gaps_matrix_final.txt", header=0, sep='\t', index_col=0)
    seqs = SeqIO.parse("refseq_18S_final.fasta", "fasta")

    records = []

    for s in seqs:
        current_id = s.id
    
        current_N = np.sum( np.array(s.seq) == "N" )
        current_gaps = np.array( data[current_id] ) 

        if (int(current_N) == 0 and np.mean(current_gaps) < theta) or ("NMR" in current_id): 		# NMR seq is from scaffold, but keep it anyways
            records.append(s)


    SeqIO.write(records, fastaOUT, "fasta") 




    



def ssu_aln(fastaIN, alnOUT):
    """
    align filtered sequences
    """

    TMPDIR = 'TMPSSU'

    cmd_ssualign = 'ssu-align --dna -f' + ' ' + fastaIN + ' ' + TMPDIR 
    subprocess.call(cmd_ssualign, shell=True)

    cmd_ssumask1 = 'ssu-mask' + ' ' + TMPDIR
    subprocess.call(cmd_ssumask1, shell=True)

    cmd_ssumask2 = 'ssu-mask --stk2afa' + ' ' + TMPDIR
    subprocess.call(cmd_ssumask2, shell=True)

    aln = AlignIO.parse(TMPDIR + "/" + TMPDIR+".eukarya.mask.stk", "stockholm")
    AlignIO.write(aln, alnOUT, "phylip")





if __name__ == "__main__":

    filterandparse_sequences("tmp.15.fasta", 0.15) 
    ssu_aln("tmp.15.fasta", "18Saln_final_15.phy")

    filterandparse_sequences("tmp.20.fasta", 0.20) 
    ssu_aln("tmp.20.fasta", "18Saln_final_20.phy")
    
    filterandparse_sequences("tmp.25.fasta", 0.25) 
    ssu_aln("tmp.25.fasta", "18Saln_final_25.phy")

    filterandparse_sequences("tmp.30.fasta", 0.30) 
    ssu_aln("tmp.30.fasta", "18Saln_final_30.phy")

    filterandparse_sequences("tmp.50.fasta", 0.50) 
    ssu_aln("tmp.50.fasta", "18Saln_final_50.phy")

    filterandparse_sequences("tmp.100.fasta", 1) 
    ssu_aln("tmp.100.fasta", "18Saln_final_100.phy")

