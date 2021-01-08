import time
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

    parser.add_argument('--results_list',required=True, help='Blast out csv from local search', type=str)
    parser.add_argument('--out_dir', required=True, help='out_directory to store results', type=str)

    args = parser.parse_args()

    return args

def get_species(results_csv, out_dir):
    if os.path.exists(results_csv) == False:
        bassy_name = os.path.basename(results_csv)

        bassy_name = re.split("_", bassy_name, maxsplit=2)[:2]

        bassy_name = bassy_name[0] + "_" + bassy_name[1]
        print(results_csv)

        sys.exit("No Blast results for this isolate: %s" % bassy_name)

    elif os.stat(results_csv).st_size == 0:
        bassy_name = os.path.basename(results_csv)

        bassy_name = re.split("_", bassy_name, maxsplit=2)[:2]

        bassy_name = bassy_name[0] + "_" + bassy_name[1]
        print("No Blast results for this isolate: %s" % bassy_name)
        sys.exit(0)

    elif os.stat(results_csv).st_size != 0:

       blast_results_csv = pandas.read_csv(results_csv,
                                  names=['subjectid','qstart','qend','sstart',
                                          'send','pident','align','evalue',
                                          'bitscore'])

    end_file_check = time.perf_counter()

    ###############################################################################
    ## Now we'll subset the results to include only the top 50% of hits ############
    ###############################################################################
    start_subset = time.perf_counter()
    blast_results_csv = blast_results_csv.sort_values('bitscore', ascending=False)
    top_bitscore = blast_results_csv.iloc[0, 8]
    top_5_percent = top_bitscore - (top_bitscore * 1)

    blast_results_csv = blast_results_csv[blast_results_csv['bitscore'] >= top_5_percent]
    end_file_check = time.perf_counter()

    ###############################################################################
    ## Now we'll extract the species of these hits ################################
    ###############################################################################


    start_subset = time.perf_counter()
    species_hits = blast_results_csv['subjectid']

    species = []
    bitscore = []

    for k in range(len(species_hits.index)):
        current_spec = species_hits.iloc[k]
        spec = os.path.basename(current_spec)
        bitscore_current = blast_results_csv.iloc[k, 8]
        bitscore.append(bitscore_current)

        species.append(spec)

    end_subset = time.perf_counter()


    bassy_name = os.path.basename(results_csv)
    if bassy_name.__contains__("!"):
        bassy_name = re.split("_control",bassy_name, maxsplit=2)[:1]
        bassy_name = bassy_name[0]
    else:
        bassy_name = re.split("_",bassy_name, maxsplit=2)[:2]
        bassy_name = bassy_name[0] + "_" + bassy_name[1]


    output_df = pandas.DataFrame()
    isolates = numpy.repeat(bassy_name, len(species))
    output_df['isolate_id'] = pandas.Series(data=isolates)
    output_df['species'] = pandas.Series(data=species, index=output_df.index)
    output_df['bitscore'] = pandas.Series(data= bitscore, index = output_df.index)

    out_name = out_dir + bassy_name + "_species_list.csv"


    output_df.to_csv(out_name,
                     index=False)

###############################################################################
## First we'll check whether the current results file is empty or not, then ###
## load it up if it is ########################################################
###############################################################################
if __name__ == '__main__':
    species_start = time.perf_counter()
    input_args = get_options()

    res_list = input_args.results_list

    out_dir = input_args.out_dir

    start_file_check = time.perf_counter()

    with open(res_list) as myfile:
        blast_res = myfile.read().splitlines()

    for file in blast_res:
        get_species(file, out_dir)
    species_end = time.perf_counter()
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Species list creator done ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Species finder for BLAST hits took: %s (seconds)" % (species_end - species_start))
