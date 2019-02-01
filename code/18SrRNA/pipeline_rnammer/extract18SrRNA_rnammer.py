 # This script run rnammer software by analysing all genome files.
import os, sys
import pandas as pd
import numpy as np
import glob 
import subprocess

# requires to have a functioning install of RNAmmer in the same directory or path



def run_rnammer(genome):
    """
    run RNAmmer on annotated genomic RNA sequences
    """

    prefix = genome.split('/')[-1].split('.')[0]

    if not os.path.exists('result'):
        os.makedirs('result')


    cmd_rnammer = "perl rnammer -S euk -m ssu -gff result/"+prefix+"_gff -h result/"+prefix+"_report -f result/"+prefix+"_18S_output.fasta " + genome 
    subprocess.call(cmd_rnammer, shell=True)

    if os.path.exists("result/"+prefix+"_18S_output.fasta"):
        print("18S rRNA sequence written to file", prefix+"_18S_output.fasta")


    if os.path.exists("temp.*"):
        os.remove("temp.*")
    if os.path.exists("1000000*"):
        os.remove("1000000*")
    if os.path.exists("*.hmmsearchresult"):
        os.remove("*.hmmsearchresult")





if __name__ == '__main__':


    list_genomes = glob.glob('genomes/*.fna')
    #list_genomes = [line.rstrip('\n') for line in open('list_genomes.txt')]

    for ix, i in enumerate(list_genomes):
        print("Analyzing", i)
        run_rnammer(list_genomes[ix])


        if ix % 10 == 0:
            print(ix, "genomes analyzed")

