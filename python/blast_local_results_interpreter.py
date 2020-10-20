import pandas
import os
import sys
import re
import numpy
import argparse

def get_options():
    purpose = '''This file will take an input from local blast output and look to return
     the species of the top 5% of bitscoe hits. Usage:
     blast_local_results_interpreter.py blast_results_csv out_result'''

    parser = argparse.ArgumentParser(description=purpose, prog='blast_local_results_interpreter.py')

    parser.add_argument('--results_csv',required=True, help='Blast out csv from local search', type=str)
    parser.add_argument('--out_dir', required=True, help='out_directory to store results', type=str)

    args = parser.parse_args()

    return args

###############################################################################
## First we'll check whether the current results file is empty or not, then ###
## load it up if it is ########################################################
###############################################################################
if __name__ == '__main__':

    input_args = get_options()

    if os.path.exists(input_args.results_csv) == False:
        bassy_name = os.path.basename(input_args.results_csv)

        bassy_name = re.split("_", bassy_name, maxsplit=2)[:2]

        bassy_name = bassy_name[0] + "_" + bassy_name[1]
        print(input_args.results_csv)

        sys.exit("No Blast results for this isolate: %s" % bassy_name)

    elif os.stat(input_args.results_csv).st_size == 0:
        sys.exit("This is an empty file")

    elif os.stat(input_args.results_csv).st_size != 0:

       results_csv = pandas.read_csv(input_args.results_csv,
                                  names=['subjectid','qstart','qend','sstart',
                                          'send','pident','align','evalue',
                                          'bitscore'])

    ###############################################################################
    ## Now we'll subset the results to include only the top 5% of hits ############
    ###############################################################################

    results_csv = results_csv.sort_values('bitscore', ascending=False)
    top_bitscore = results_csv.iloc[0, 8]
    top_5_percent = top_bitscore - (top_bitscore / 100 * 5)

    results_csv = results_csv[results_csv['bitscore'] >= top_5_percent]

    ###############################################################################
    ## Now we'll extract the species of these hits ################################
    ###############################################################################

    species_hits = results_csv['subjectid']

    species = []
    bitscore = []

    for k in range(len(species_hits.index)):
        current_spec = species_hits.iloc[k]
        spec = os.path.basename(current_spec)
        bitscore_current = results_csv.iloc[k, 8]
        bitscore.append(bitscore_current)

        species.append(spec)

    bassy_name = os.path.basename(input_args.results_csv)
    bassy_name = re.split("_",bassy_name, maxsplit=2)[:2]
    bassy_name = bassy_name[0] + "_" + bassy_name[1]


    output_df = pandas.DataFrame()
    isolates = numpy.repeat(bassy_name, len(species))
    output_df['isolate_id'] = pandas.Series(data=isolates)
    output_df['species'] = pandas.Series(data=species, index=output_df.index)
    output_df['bitscore'] = pandas.Series(data= bitscore, index = output_df.index)

    out_name = input_args.out_dir + bassy_name + "_species_list.csv"


    output_df.to_csv(out_name,
                     index=False)