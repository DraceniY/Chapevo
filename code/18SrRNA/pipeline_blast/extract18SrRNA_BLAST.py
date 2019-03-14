import os, sys
import numpy as np
import pandas as pd
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



# this runs on the annotated genomic RNA files from RefSeq

# set to directory where the genomic annotated RNA fasta files from RefSeq are stored
PATH_TO_GENOMES = './rna/' 

# set path to local BLAST install
BLASTN = '~/sw/ncbi-blast-2.2.31+/bin/blastn'



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def run_blast(genome):
    """
    blast input genome against eukaryotic 18S sequences
    """

    prefix = genome.split('/')[-1].split('.')[0]

    if not os.path.exists('results'):
        os.makedirs('results')

    cmd_blastn = BLASTN + ' ' + '-query' + ' ' + genome + ' ' + '-db ../../../data/18SrRNA/pipeline_blast/blast/18S_Euka -outfmt 6 -out' + ' ' 'results/'+prefix+'.blast -evalue 1e-6' 
    subprocess.call(cmd_blastn, shell=True)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_output(genome):
    """
    read in BLAST output, and extract 18S rRNA sequence (if found)
    """

    prefix = genome.split('/')[-1].split('.')[0]
    results_dir = 'results/'

    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    blastout = results_dir + prefix + '.blast'


    out = None

    if os.path.exists(blastout) and os.stat(blastout).st_size > 0 :

        data = pd.read_table(blastout , header=None)
        hits = list(data[0].unique())

        genome_fasta = SeqIO.to_dict(SeqIO.parse(open(genome),'fasta'))

        records = []
        for i in hits:
            current_rec =  genome_fasta[i] 
            current_des = current_rec.description

            if '18S' in current_des and 'protein' not in current_des and 'mRNA' not in current_des:
                records.append(current_rec)
                    
        SeqIO.write(records, results_dir + prefix + "_18SrRNA.fasta", "fasta") 


    return out



def collect_sequences(fasta_dir):
    """
    merges all extracted putative 18S rRNA sequences
    """

    list_seqs = glob.glob(fasta_dir + '/*.fasta')

    records = []
    for i in list_seqs:
        if os.path.exists(i) and os.stat(i).st_size > 0:
            current_seqs = SeqIO.parse(i, "fasta")
        
            added = False 
            for s in current_seqs:
                if not added:
                    records.append(s)
                    added = True

    SeqIO.write(records, "blast_18SrRNA.fasta", "fasta")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    result_df = pd.DataFrame(columns=['genome', '18S'])
    list_genomes = glob.glob(PATH_TO_GENOMES+"/*.fna")


    for ix, i in enumerate(list_genomes):
        print("Analyzing", i)
        prefix = i.split('/')[-1].split('.')[0]

        run_blast(list_genomes[ix])
        out_status = parse_output(list_genomes[ix]) 

        result_df.loc[ix] = (prefix, out_status)


    collect_sequences('results')
    result_df.to_csv('blast_result.txt', sep='\t', header=True, index=False)
