import os , sys , glob, argparse , itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import subprocess






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def create_index(genome_fasta, nameOUT):
    """
    creates HISAT index for target genome
    aligns synthetic rRNA reads
    """



    SYN_rRNA_FASTA = '../../../data/pipeline_mapping_reads_genome/fastq_files/reads_total.fastq' 


    # extract species prefix
    prefix = genome_fasta.split("/")[-1].split(".")[0]


    # build Hisat2 index for target genome
    cmd_index_build = './hisat2-build' + ' ' +  genome_fasta + ' ' + prefix
    subprocess.call(cmd_index_build, shell=True)

    # align synthetic rRNA mix to target genome
    if not os.path.exists('TMP'):
        os.makedirs('TMP')

    cmd_align = './hisat2' + ' ' + '-f' + ' ' + '-a' + ' ' + '-x' + prefix + ' ' + '-U' + ' ' + SYN_rRNA_FASTA + ' ' + '>' + ' ' + nameOUT
    subprocess.call(cmd_align, shell=True) 

    # remove index
    if os.path.exists(prefix + '*.ht2'): 
        os.remove(prefix + '*.ht2') 





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_sam_file(sam_file):
    """
    load SAM file into data frame
    """

    sam_df = pd.DataFrame(columns=['reads', 'position', 'chromosome'])
    counter=0

    with open(sam_file, "r") as f:

        for line in f:
            if line[0] != "@" and line[0] != "/":
                current_line = line.split()
                chromo = current_line[2]
                pos = int(current_line[3])
                read= current_line[0]

                if pos > 0 and chromo != "*":
                    sam_df.loc[counter]= (read , pos , chromo)
                    counter += 1

                    if counter % 1000 == 0:
                        print(counter, "reads read into DF")
    return sam_df







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_sorted_sam_file(sam_file):
    """
    load SAM file into data frame
    """

    best_df = None
    best_cov = 0
    last_chromo = None
    
    sam_df = pd.DataFrame(columns=['reads', 'position', 'chromosome'])
    counter=0

    with open(sam_file, "r") as f:

        for line in f:
            if line[0] != "@" and line[0] != "/":
                current_line = line.split()
                chromo = current_line[2]
                pos = int(current_line[3])
                read= current_line[0]

                if chromo == last_chromo: 
                
                    if pos > 0 and chromo != "*":            
                        sam_df.loc[counter]= (read , pos , chromo)
                        last_chromo = chromo
                        counter += 1
                   
                elif chromo != last_chromo:
                    if len(sam_df['position']) > best_cov:
                        best_df = sam_df
                        best_cov = len(best_df['position'])

                    sam_df = pd.DataFrame(columns=['reads', 'position', 'chromosome'])
                    counter = 0
                    sam_df.loc[counter] = (read, pos, chromo)
                    last_chromo = chromo
                    counter += 1

    if len(sam_df['position']) > best_cov:
        best_df = sam_df
 
    return best_df





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def analyze_sam(sam_df):
    """
    extracts regions of candidate 18S rRNA with best 'coverage'
    """

    list_chromosome = set(sam_df["chromosome"].values)

    best_chromosome = "None"
    best_coverage = 0
    best_start = 10000000000
    best_end = 0
    best_len = -1

    nchromo = 0
    for chromo in list_chromosomes:

        current_chr = sam_df[sam_df["chromosome"] == chromo]
        nchromo += 1
        print("Chromosome nr. ", nchromo, "is being analyzed")

        # get max , min and length of chromosome
        current_end   = np.max( np.array(current_chr["position"] ) )
        current_start = np.min( np.array(current_chr["position"] ) )
        current_len   = int(current_end) - int(current_start)
        nb_reads = len(set(current_chr["reads"].tolist()))

        if nb_reads > best_coverage and current_len > 500:        
            best_coverage = np.copy(nb_reads)
            best_chromosome = chromo
            best_start = np.copy(current_start)
            best_end  = np.copy(current_end)
            best_len  = int(best_end) - (best_start)



    if best_chromosome != "None" and best_len > 0 and best_len < 3000 :

        # add 500 bases to start and end position
        start = int(best_start) - np.minimum(500, int(best_start) ) 
        end   = int(best_end) + 500


    elif best_chromosome != "None" and best_len > 3000:
        # find a smaller candidate region, here e.g. of max 3000nt
        best_df = sam_df[sam_df["chromosome"] == best_chromosome]
        current_pos = np.array(best_df['position'])
        new_start = best_start
        new_end = best_end
        new_coverage = 0
 
        for i in range(int(best_start), int(best_end) - 3000):
            sel_lo = current_pos > i
            sel_hi = current_pos < i + 3000
            sel = sel_lo * sel_hi
            no_reads = len(current_pos[sel])

            if no_reads > new_coverage:
                new_start = i
                new_end = i + 3000
                new_coverage = np.copy(no_reads)

        start = new_start
        end   = new_end


    else:
        best_chromosome = 'none'
        start = 'none'
        end   = 'none' 

    print(best_chromosome, start, end)
    return best_chromosome, start, end




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def analyze_sorted_sam(sam_df):
    """
    extracts regions of candidate 18S rRNA with best 'coverage'
    """

    best_chromosome = list(set(list(sam_df["chromosome"])))
    best_chromosome =  best_chromosome[0]


    best_coverage = 0
    best_start = 10000000000
    best_end = 0
    best_len = -1


    # get max , min and length of chromosome
    current_end   = np.max( np.array(sam_df["position"] ) )
    current_start = np.min( np.array(sam_df["position"] ) )
    current_len   = int(current_end) - int(current_start)
    nb_reads = len(set(sam_df["reads"].tolist()))

    if nb_reads > best_coverage and current_len > 500:
        best_coverage = np.copy(nb_reads)
        best_start = np.copy(current_start)
        best_end  = np.copy(current_end)
        best_len  = int(best_end) - (best_start)


    if best_len > 0 and best_len < 3000 :

        # add 500 bases to start and end position
        start = int(best_start) - np.minimum(500, int(best_start) )
        end   = int(best_end) + 500

    elif best_len > 3000:
        # find a smaller candidate region, here e.g. of max 3000nt
        current_pos = np.array(sam_df['position'])
        new_start = best_start
        new_end = best_end
        new_coverage = 0

        for i in range(int(best_start), int(best_end) - 3000):
            sel_lo = current_pos > i
            sel_hi = current_pos < i + 3000
            sel = sel_lo * sel_hi
            no_reads = len(current_pos[sel])

            if no_reads > new_coverage:
                new_start = i
                new_end = i + 3000
                new_coverage = np.copy(no_reads)

        start = new_start
        end   = new_end


    else:
        best_chromosome = 'none'
        start = 'none'
        end   = 'none'


    return best_chromosome, start, end




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def run_pipeline(genome_file):
    """
    get 18S rRNA from genome fasta
    """

    prefix = genome_file.split("/")[-1].split(".")[0]
    sam_file = 'TMP/'+prefix+'.sam'

    # create index and align rRNA reads
    if not os.path.exists(sam_file):
        create_index(genome_file, sam_file)

    if os.path.exists(sam_file) and os.path.getsize(sam_file) > 0:

        cmd_sort = 'samtools sort' + ' ' + sam_file + ' ' + '-o tmp.sam'
        subprocess.call(cmd_sort, shell=True)

        # parse SAM file
        #current_sam = parse_sam_file(sam_file)
        current_sam = parse_sorted_sam_file('tmp.sam')

        
        #print("Identify chromosomal region with highest coverage")
        #chromo, start, end = analyze_sam(current_sam)
        chromo, start, end = analyze_sorted_sam(current_sam)

        if chromo != 'none':

            # load genome fasta and extract candidate seq
            current_fasta = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
            print(chromo)
            candidate_seq = current_fasta[chromo][start:end]
            SeqIO.write(candidate_seq , prefix + "_candidate.fasta" , "fasta")

            if not os.path.exists('result'):
                os.makedirs('result')

            # run RNAmmer
            print("Running RNAmmer")
            if os.path.exists(prefix+"_candidate.fasta"):

                cmd_rnammer = "perl rnammer -S euk -m ssu -gff result/"+prefix+"_gff -h result/"+prefix+"_report -f result/"+prefix+"_18S_output.fasta "+prefix+"_candidate.fasta"
                subprocess.call(cmd_rnammer, shell=True)

                print("18S rRNA sequence written to file", prefix+"_18S_output.fasta") 

                if os.path.exists("temp.*"):
                    os.remove("temp.*")
                if os.path.exists("1000000*"):
                    os.remove("1000000*")
                if os.path.exists("*.hmmsearchresult"):
                    os.remove("*.hmmsearchresult")

    else:
        print("no rRNA found")


    if os.path.exists('tmp.sam'):
        os.remove('tmp.sam')


if __name__ == '__main__':


    list_genomes = glob.glob('genomes/*.fna')

    for ix, i in enumerate(list_genomes):
        print("Analyzing", i)
        run_pipeline(i)

 
        if ix % 10 == 0: 
            print(ix, "genomes analyzed")
