import re
import pandas
import numpy
import time
import sys
import subprocess
import argparse


def get_options():
    purpose = '''This is a python script to intake a the embl branch recon file and output
    a 5 column table of: ['start_node','end_node','start_base','end_base','base_number'] .
    Usage: changing_embl_branch_recon_to_table.py <input_embl> <output_csv> '''
    parser = argparse.ArgumentParser(description=purpose, prog='embl_base_csv.py')

    parser.add_argument('--embl_file', required=True, help='Embl branch base from gubbins', type=str)
    parser.add_argument('--out_name', required=True, help='Out csv filename (required)', type=str)

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    input_args = get_options()
    ## get the size of the new table

    grep_command = "grep -c variation " + input_args.embl_file + " > grep_variation.txt"
    subprocess.run(grep_command, shell=True)

    with open("grep_variation.txt") as grepper:
        branch_and_bases = grepper.readlines()

    branches_and_bases = int(branch_and_bases[0])


    embl_file = open(input_args.embl_file, "r")



    new_table = pandas.DataFrame(data = numpy.zeros(shape=(branches_and_bases, 5)),
                                 columns=['start_node','end_node', 'start_base','end_base',
                                          'base_number'])

    counter = 0

    test_series = pandas.Series(data="G", index=['start_base'])

    print("")
    print("Beginning run through")
    print("")

    for line in embl_file:
        start_time = time.time()
        if re.search("variation", line):
            base_number = re.split("variation", line)
            base_number_only = base_number[1].rstrip()
            base_number_only = re.sub(" ", "", base_number_only)
            base_number_only = int(base_number_only)
            new_table.iloc[counter, 4] = base_number_only

        elif re.search("node", line):
            nodes_only = re.split("=\"",line)
            indiv_nodes = re.split("->",nodes_only[1])
            indiv_nodes[1] = indiv_nodes[1].rstrip()
            indiv_nodes[1] = indiv_nodes[1][:-1]
            indiv_nodes = pandas.Series(indiv_nodes, index=['start_node','end_node'])
            new_table.iloc[counter, [0, 1]] = indiv_nodes

        elif re.search("parent_base", line):
            parents_only = re.split("=\"", line)
            indiv_base = parents_only[1]
            indiv_base = indiv_base.rstrip()
            indiv_base = indiv_base[:-1]
            #indiv_base = pandas.Series(data = str(indiv_base), index=['start_base'])
            new_table.iloc[counter, 2] = indiv_base

        elif re.search("replace", line):
            new_only = re.split("=\"", line)
            replace_base = new_only[1]
            replace_base = replace_base.rstrip()
            replace_base = replace_base[:-1]
            #replace_base = pandas.Series(data = replace_base, index=['end_base'])
            new_table.iloc[counter, 3] = replace_base
            counter += 1
            end_time = time.time()
            duration = end_time - start_time
            time_left = duration * (branches_and_bases - counter)
            if counter > 0:
                if counter % 100 == 0:
                    print((counter / branches_and_bases) * 100, "%  ","This long left:", time_left, end='\r', flush=True)
                    
                    new_table.to_csv(input_args.out_name,
                                     index=False)


    new_table.to_csv(input_args.out_name,index=False)






