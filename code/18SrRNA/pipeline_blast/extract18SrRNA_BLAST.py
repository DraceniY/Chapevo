import os, sys
import numpy as np
import pandas as pd
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



# this runs on the annotated genomic RNA files from RefSeq


PATH_TO_GENOMES = '../genomes'
BLASTN = '~/sw/ncbi-blast-2.2.31+/bin/blastn'



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def run_blast(genome):
    """
    blast input genome against eukaryotic 18S sequences
    """

    prefix = genome.split('/')[-1].split('.')[0]

    if not os.path.exists('results'):
        os.makedirs('results')

    cmd_blastn = BLASTN + ' ' + '-query' + ' ' + genome + ' ' + '-db ../../../data/pipeline_blast/blast/18S_Euka -outfmt 6 -out' + ' ' 'results/'+prefix+'.blast -evalue 1e-6' 
    subprocess.call(cmd_blastn, shell=True)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_output(genome):
    """
    read in BLAST output, and extract 18S rRNA sequence (if found)
    """

    prefix = genome.split('/')[-1].split('.')[0]

    blastout = prefix+'.blast'

    if not os.path.exists('results'):
        os.makedirs('results')


    out = None

    if os.path.exists(blastout):
        if(os.stat(blastout).st_size > 0): # Check if the file is empty

            data = pd.read_table(blastout , header=None)
            hits = list(data[0].unique())

            genome_fasta = SeqIO.to_dict(SeqIO.parse(open(genome),'fasta'))

            records = []
            for i in hits:
                current_rec =  genome_fasta[i] 
                current_des = current_rec.description

                if '18S' in current_des and 'protein' not in current_des and 'mRNA' not in current_des:
                    records.append(current_rec)
            
            SeqIO.write(records, prefix+"_18SrRNA.fasta", "fasta") 

            out = "yes"
        else:
            out = "no"

    else:
        out = "no"

    return out






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



    result_df.to_csv('blast_analysis.txt', sep='\t', header=True, index=False)
