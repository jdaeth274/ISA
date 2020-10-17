import re
import pandas
import numpy
import time
import sys
import argparse


def get_options():

    purpose = '''This is a python script to intake a the embl recombination prediction file and output
    a 5 column table of: ['start_node','end_node','start_rec','end_rec','snp_number'] .
    Usage: changing_gubbines_recomb_embl_to_table.py <input_embl> <output_csv> '''
    parser = argparse.ArgumentParser(description=purpose, prog='embl_recombinations.py')

    parser.add_argument('--embl_file', required=True, help='Embl branch base from gubbins', type=str)
    parser.add_argument('--out_name', required=True, help='Out csv filename (required)', type=str)

    args = parser.parse_args()

    return args



if __name__ == '__main__':

    input_args = get_options()
    embl_file = open(input_args.embl_file, "r")

    embl_copy = embl_file
    branches_and_bases = 0

    for line in embl_file:
        if re.search("misc",line):
            branches_and_bases += 1

    embl_file.close()

    embl_file = open(input_args.embl_file, "r")

    print("There are %s recombination blocks" % branches_and_bases)

    new_table = pandas.DataFrame(data = numpy.zeros(shape=(branches_and_bases, 5)),
                                 columns=['start_node','end_node', 'start_rec','end_rec',
                                          'snp_number'])

    counter = 0
    print("")
    print("Running through embl file")
    print("")
    for line in embl_file:
        start_time = time.time()
        if re.search("misc_feature", line):
            base_number = re.split("misc_feature", line)
            base_number_only = base_number[1].rstrip()
            base_number_only = re.sub(" ", "", base_number_only)
            base_numbers = re.split("\.\.",base_number_only)
            base_numbers = pandas.Series(base_numbers, index=['start_rec','end_rec'])
            new_table.iloc[counter, [2,3]] = base_numbers

        elif re.search("node", line):
            nodes_only = re.split("=\"",line)
            indiv_nodes = re.split("->",nodes_only[1])
            indiv_nodes[1] = indiv_nodes[1].rstrip()
            indiv_nodes[1] = indiv_nodes[1][:-1]
            indiv_nodes = pandas.Series(indiv_nodes, index=['start_node','end_node'])
            new_table.iloc[counter, [0, 1]] = indiv_nodes

        elif re.search("SNP", line):
            snps_only = re.split("=\"", line)
            indiv_snp = snps_only[1]
            indiv_snp = indiv_snp.rstrip()
            indiv_snp = indiv_snp[:-1]
            indiv_snp = int(indiv_snp)
            #indiv_base = pandas.Series(data = str(indiv_base), index=['start_base'])
            new_table.iloc[counter, 4] = indiv_snp

            counter += 1
            end_time = time.time()
            duration = end_time - start_time
            time_left = duration * (branches_and_bases - counter)
            if counter > 0:
                if counter % 100 == 0:
                    print("This long left: ", time_left, " ", (counter / branches_and_bases) * 100, "%  ", end='\r', flush=True)
                    new_table.to_csv(input_args.out_name,
                                     index=False)


    new_table.to_csv(input_args.out_name,
                     index=False)

