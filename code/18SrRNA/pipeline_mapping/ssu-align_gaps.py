import os, sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pandas as pd
import itertools as it
import subprocess
import glob

# expects a local copy of ssu-align in the path


def gaps_ssualign(fastaIN):
    """
    runs ssu-align on input fasta, returns number of gaps in sequence alignment
    """

    TMPDIR = 'TMP_SSU'

    cmd_ssualign = 'ssu-align -f' + ' ' + fastaIN + ' ' + TMPDIR + ' ' + '> log'
    subprocess.call(cmd_ssualign, shell=True)

    cmd_ssumask1 = 'ssu-mask' + ' ' + TMPDIR + ' ' + '> log'
    subprocess.call(cmd_ssumask1, shell=True)

    cmd_ssumask2 = 'ssu-mask --stk2afa' + ' ' + TMPDIR + ' ' + '> log'
    subprocess.call(cmd_ssumask2, shell=True)

    aln_file = TMPDIR + '/' + TMPDIR + '.eukarya.afa'    

    gaps1 = 0
    gaps2 = 0
    if os.path.exists(aln_file):
        alignment = list(SeqIO.parse(aln_file, 'fasta'))

        seq1 = alignment[0]
        seq1 = np.array(seq1.seq)
        gaps1 += sum( seq1 == '.')
        gaps1 += sum( seq1 == '-')
        gaps1 /= float(len(seq1))
        gaps1 = round(gaps1, 2)
        
        seq2 = alignment[1]
        seq2 = np.array(seq2.seq)
        gaps2 += sum( seq2 == '.')
        gaps2 += sum( seq2 == '-')
        gaps2 /= float(len(seq2))
        gaps2 = round(gaps2, 2)

    else:
        gaps1 = np.nan
        gaps2 = np.nan

    return gaps1, gaps2




if __name__ == "__main__":

    sequences = SeqIO.to_dict(SeqIO.parse("refseq_18S_final.fasta", "fasta"))
    seqids = list(sequences.keys())
    sequence_pairs = list(it.combinations(np.arange( len(seqids) ), 2) )

    gapmat = np.zeros(( len(seqids), len(seqids) ), dtype=float)

    print( len(sequence_pairs) )
    counter = 0
    for idx1, idx2 in sequence_pairs:
        seqid1 = seqids[idx1]
        seqid2 = seqids[idx2]

        # 1. write sequence pair to unaligned FASTA
        records = []
        records.append(sequences[seqid1])
        records.append(sequences[seqid2])
        SeqIO.write(records, "tmp.fasta", "fasta")
        
        # 2. align sequences with SSU-ALIGN and get gaps
        gaps1, gaps2 = gaps_ssualign('tmp.fasta')
        print(gaps1, gaps2, idx1, idx2)
        gapmat[idx1, idx2] = float(gaps1)
        gapmat[idx2, idx1] = float(gaps2)

        if os.path.exists("tmp.fasta"):
            os.remove("tmp.fasta") 


        counter += 1

        if counter % 1000 == 0:
            print(counter, "of", len(sequence_pairs), "sequence pairs analyzed")

    print(gapmat)
    gap_df = pd.DataFrame(data=gapmat, index=seqids, columns=seqids )
    gap_df.to_csv("gaps_matrix_final.txt", sep='\t', header=True)




