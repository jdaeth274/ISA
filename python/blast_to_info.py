import pandas
import os
import sys
import re
import numpy
import argparse

def get_options():

    purpose = '''This is a script to take the blast species results csvs 
    and then from there create a more informative outcome.
    Usage: blast_results_to_info.py <species_list_csvs> <python_hit_locs_csv> out_results'''

    parser = argparse.ArgumentParser(prog='blast_to_info.py', description=purpose)

    parser.add_argument('--list_file', required=True, help='Species list ot from blast local results interprete', type=str)
    parser.add_argument('--hit_locs', required=True, help='Hits df from hit allocator csv', type=str)
    parser.add_argument('--out_name', required=True, help='out_name for csv', type=str)

    input_args = parser.parse_args()

    return input_args

###############################################################################
## Now we'll look at getting the species results from the different results ###
###############################################################################


if __name__ == '__main__':

    input_args = get_options()
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Compiling BLAST results now ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


    list_file = open(input_args.list_file, "r")
    species_list = list_file.read().splitlines()

    python_hit_locs = pandas.read_csv(input_args.hit_locs)


    isolates = []
    top_species = []
    bitscores = []
    insertion_point = []

    for k in range(len(species_list)):
        current_file = species_list[k]
        current_csv = pandas.read_csv(current_file)
        top_hit = current_csv.iloc[0]
        isolates.append(top_hit.iloc[0])
        top_species.append(top_hit.iloc[1])
        bitscores.append(top_hit.iloc[2])

        izzy = top_hit.iloc[0]

        if "!" in izzy:
            insertion_point.append("reference")
        else:
            matching_python_record = python_hit_locs[python_hit_locs['id'] == top_hit.iloc[0]]
            clus_name = matching_python_record['insert_name'].values[0]


            #insertion_point.append(str(matching_python_record.iloc[0:1,matching_python_record.columns.get_loc("cluster_names")]))
            insertion_point.append(str(clus_name))


    out_df = pandas.DataFrame()
    out_df['isolate'] = pandas.Series(data=isolates)
    out_df['top_species'] = pandas.Series(data=top_species, index=out_df.index)
    out_df['bitscore'] = pandas.Series(data=bitscores, index=out_df.index)
    out_df['insertion_point'] = pandas.Series(data=insertion_point, index=out_df.index)

    out_df.to_csv(input_args.out_name,index=False)

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Finished ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

