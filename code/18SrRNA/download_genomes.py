import os, sys
import numpy as np
import pandas as pd
import glob
import subprocess


# set UNIPROT_DIR to local directory that contains all uniprot proteome fasta files (only those without isoforms)
UNIPROT_DIR = './uniprot/' 

# set GENOME_DIR to local directory where all RefSeq genomes will be stores (several GB)
GENOME_DIR = './genome/'


list_taxid = []
with open("../../data/taxonomy/taxid.txt", 'r') as file:
    for line in file:
        current_line = line.split()
        current_taxid = current_line[0]
        list_taxid.append(int(current_taxid))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def retrieve_refseq_summary():
    """
    download the RefSeq assembly summary and load into DF
    """

    if not os.path.exists("assembly_summary_refseq.txt"):
        subprocess.call('wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt', shell=True)

    refseq = pd.read_csv("assembly_summary_refseq.txt", sep='\t', header=1, index_col=False)

    refseq = refseq[['# assembly_accession', 'organism_name', 'taxid', 'species_taxid', 'ftp_path']]


    return refseq



rs = retrieve_refseq_summary()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def refseq_with_uniprot(refseqDF):
    """
    find refseq genomes to uniprot proteoms
    """

    refseqDF = refseqDF.drop_duplicates(subset='taxid', keep="first")


    download_df = pd.DataFrame(columns=['taxid', 'taxid2', 'species_taxid', 'refseq', 'uniprot', 'ftp'])
    counter = 0

    list_uniprot = glob.glob(UNIPROT_DIR + '*.fasta')

    for i in list_uniprot:

        current_name = i.split('/')[-1]
        current_uniprot = current_name.split('_')[0]
        current_taxid = current_name.split('_')[1].split('.')[0]


        current_tx1 = refseqDF[refseqDF['taxid']==int(current_taxid)]
        current_tx2 = refseqDF[refseqDF['species_taxid']==int(current_taxid) ]

        if len(current_tx1.index) > 0:

            refseq_accession = current_tx1['# assembly_accession'].item()
            refseq_taxid = current_tx1['taxid'].item()
            refseq_speciestaxid = current_tx1['species_taxid'].item()
            refseq_ftp = current_tx1['ftp_path'].item()
            if refseq_taxid in list_taxid or refseq_speciestaxid in list_taxid:
                download_df.loc[counter] = [current_taxid, refseq_taxid, refseq_speciestaxid, refseq_accession, current_uniprot, refseq_ftp  ]
                counter += 1

        elif len(current_tx2.index) > 0:

            refseq_accession = current_tx2['# assembly_accession'].item()
            refseq_taxid = current_tx2['taxid'].item()
            refseq_speciestaxid = current_tx2['species_taxid'].item()
            refseq_ftp = current_tx2['ftp_path'].item()
            if refseq_taxid in list_taxid or refseq_speciestaxid in list_taxid:
                download_df.loc[counter] = [current_taxid, refseq_taxid, refseq_speciestaxid, refseq_accession, current_uniprot, refseq_ftp  ]
                counter += 1


    download_df = download_df.drop_duplicates('refseq')	# use unambigous matches
    download_df.to_csv("refseq_download.txt", sep='\t', header=True )

    return download_df

df = refseq_with_uniprot(rs)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def download_genomes(DF):
    """
    downloads genome files from NCBI site
    """

    RNA_dir = GENOME_DIR + '/rna/'
    GENOMIC_dir = GENOME_DIR + '/genomic/'

    if not os.path.exists(RNA_dir):
        os.makedirs(RNA_dir)
    if not os.path.exists(GENOMIC_dir):
        os.makedirs(GENOMIC_dir)

    for genome_link in DF['ftp']:
        prefix = genome_link.split('/')[-1]
        rna_fname = prefix + '_rna.fna.gz'
        genomic_fname = prefix + '_genomic.fna.gz' 

        cmd_rna     = 'wget' + ' ' + '-O' + ' ' + RNA_dir+rna_fname + ' ' + genome_link + '/' + rna_fname
        cmd_genomic = 'wget' + ' ' + '-O' + ' ' + GENOMIC_dir + genomic_fname + ' ' + genome_link + '/' + genomic_fname

        if not os.path.exists(RNA_dir + rna_fname):
            subprocess.call(cmd_rna, shell=True)

        if not os.path.exists(GENOMIC_dir + genomic_fname):
            subprocess.call(cmd_genomic, shell=True)

download_genomes(df)



