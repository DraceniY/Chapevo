import os , sys , glob, argparse , itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import subprocess



accession_table = pd.read_csv("../refseq_download.txt", sep='\t', header=0)
for i in range(len(accession_table)):
    accession_table.loc[i,'refseq'] = accession_table.loc[i,'refseq'].split('.')[0]



taxid_dict = {}
with open("../../../data/taxonomy/taxid.txt", 'r') as file:
    for line in file:
        current_line = line.split()
        current_id = int(current_line[0])
        current_name = current_line[1]
        taxid_dict[current_id] = current_name




def load_sequences(seqdir):
    """
    load all 18S sequences
    change name
    parse to one fasta file
    """

    list_seqs = glob.glob(seqdir+'/*_18S_output.fasta')

    records = []

    list_added = []

    for i in list_seqs:
        prefix = i.split('/')[-1].split('_')[0:2]
        prefix2 = "_".join(prefix)

        taxid = accession_table[accession_table['refseq']==prefix2][['taxid', 'taxid2', 'species_taxid']].values
        taxid2 = list(set(list(taxid[0]))  )

        good_taxid = []
        names = []
        for j in taxid2:
            if j in taxid_dict.keys():
                current_name = taxid_dict.get(j)
                names.append(current_name)
                good_taxid.append(j)
        names = list(set(names))
        current_name = names[0]
        current_short = "TAX" + str(good_taxid[0]) + "xxxxxxxxxx"
        current_short = current_short[0:10]
#        current_short = str(current_name.split('_')[0][0]) + str(current_name.split('_')[1][0:3]) + str(good_taxid[0])
        current_fname = 'result/' + prefix2 + '_18S_output.fasta'

        if os.path.exists(current_fname) and os.path.getsize(current_fname) > 0:
            current_fasta = SeqIO.parse(open(current_fname, mode='r'), 'fasta')
            for current_seq in current_fasta:
                current_seq.id = current_short 
                current_seq.description = "Taxid:"+str(good_taxid[0]) + " " + "Organism:"+current_name + " " + "Type:18S_rRNA" 

                if current_name not in list_added:
                    records.append(current_seq)
                    list_added.append(current_name)


    SeqIO.write(records , 'refseq_18S.fasta' , "fasta")


if __name__ == '__main__':

    load_sequences('result_new') 
