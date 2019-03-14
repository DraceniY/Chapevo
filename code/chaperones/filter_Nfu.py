import sys
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess


file_in ='Nfu_20150522.genes_20150922.proteins.fa'
file_out='Nothobranchius_furzeri.fa'

url = 'http://nfingb.leibniz-fli.de/data/raw/notho4/Nfu_20150522.genes_20150922.proteins.fa.gz'

subprocess.call('wget' + ' ' + url, shell=True)
subprocess.call('gunzip Nfu_20150522.genes_20150922.proteins.fa.gz', shell=True)

list_geneID = []

with open(file_out, 'w') as f_out:
    for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):
        seq_record.description=' '.join(seq_record.description.split()[1:])
        current_description = seq_record.description
        current_geneid = (current_description.split()[0].split('=')[1])

        if current_geneid not in list_geneID:
            r=SeqIO.write(seq_record, f_out, 'fasta')
            if r!=1: print('Error while writing sequence:  ' + seq_record.id)
            list_geneID.append(current_geneid)
