import sys, os
import numpy as np
import pandas as pd
import pickle
import glob





list_organisms = [line.strip().decode("utf-8") for line in open('../data/tree/tree_list_nodes.txt', 'rb')]





def pickle_agg(aggdir, fileOUT):
    """
    writes individual agg files into dict
    """

    files = glob.glob(aggdir + "*.txt")
    agg_dict = {}

    for i in files:
        current_rec = pd.read_csv(i, header=0, sep='\t', index_col=None)
        current_agg = np.array(current_rec['Aggregation'])
        current_name = i.split('/')[-1].split('.')[0]
        agg_dict[current_name] = current_agg

    with open(fileOUT, 'wb') as handle:
        pickle.dump(agg_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


# after running Tango aggregation prediction software on all uniprot proteomes, read in output here
files = glob.glob("aggregation/*.agg.pkl")






def summary_stats(fname):
    """
    get summary stats from pickled dict
    """
    with open(fname, 'rb') as f:
        agg = pickle.load(f)

    name = fname.split('/')[-1].split('.')[0]

    if name in list_organisms:
        treepos = int(list_organisms.index(name)) + 1
    else:
        treepos = "NA"

    L = len(agg.keys() )

    scores = np.zeros(( L ))
    for ix, i in enumerate( agg.keys() ):
        scores[ix] = np.sum( agg[i] )

    current_mean = np.mean(scores)
    current_median = np.median(scores)
    current_sd   = np.std(scores)
    current_q1   = np.percentile(scores, 25)
    current_q3   = np.percentile(scores, 75)
    current_min  = np.min(scores)
    current_max  = np.max(scores)
    current_sum  = np.sum(scores)

    return name, treepos, L, round(current_median,2), round(current_mean,2), round(current_sd,2), round(current_q1,2), round(current_q3,2), round(current_min,2), round(current_max,2), round(current_sum, 2)



#out = summary_stats(files[1])



def agg_stats(flist):
    """
    go through pickled aggregation dicts and extract summary stats
    """

    agg_stats = pd.DataFrame(columns=['Name', 'treepos', 'Nprot', 'Median', 'Mean', 'Std', 'Q1', 'Q3', 'Min', 'Max', 'Sum'])

    counter = 0
    for i in flist:
        # set directory to local DIR with output of aggregation predictions
        filename = 'aggregation/' + i + '.agg.pkl'
        out = summary_stats(filename)
        print(out)
        agg_stats.loc[counter] = out
        counter += 1

    agg_stats.sort_values('treepos', inplace=True)
    agg_stats.to_csv("aggregation_summary.txt", sep='\t', header=True, index=False)


agg_stats(list_organisms)




    
