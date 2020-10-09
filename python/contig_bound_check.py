import re
import pandas
import numpy
import sys
import os
import argparse

def get_options():

    purpose = ''' This is a script to get the bounds of contigs from an initial grep search for
    the > among a multiple fa document.
    Usage: python contig_bound_check.py <grep_output> <fasta_file>'''
    parser = argparse.ArgumentParser(description=purpose,
                                     prog='contig_bound_check.py')

    parser.add_argument('--grep_file', required=True, help='Out file from grep search of fasta file for ">" ', type=str)
    parser.add_argument('--fasta_file', required=True, help='Isolate fasta file',type=str)

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    files_for_input = get_options()

    grep_file = open(files_for_input.grep_file, "r")

    grep_file_lines = grep_file.readlines()

    contig_starts = []

    for k in range(len(grep_file_lines)):
        current_line = grep_file_lines[k]
        current_line_start = re.split(":",current_line)[0]
        current_line_start = re.sub("\'", "",current_line_start)
        current_line_start = int(current_line_start)
        contig_starts.append(current_line_start)

    ###############################################################################
    ## Now we'll get the fasta file in ############################################
    ###############################################################################

    fasta_file = open(files_for_input.fasta_file, "r")

    temp = fasta_file.read().splitlines()

    contig_bounds = pandas.DataFrame(data=numpy.zeros(shape=(len(contig_starts), 2)))

    for k in range(len(contig_starts)):
        if k == 0:
            contig_bounds.iloc[0,0] = 1
        else:
            contig_bounds.iloc[k, 0] = contig_bounds.iloc[k-1, 1] + 1

        current_start = contig_starts[k]
        if (k + 1) != len(contig_starts):
            current_end = contig_starts[k + 1] - 1
        else:
            current_end = len(temp)

        if (k + 1) != len(contig_starts):
            dist = sum(len(i) for i in temp[current_start:current_end])
            contig_bounds.iloc[k, 1] = contig_bounds.iloc[k, 0] + dist - 1
        else:
            dist = sum(len(i) for i in temp[current_start:current_end])
            contig_bounds.iloc[k, 1] = contig_bounds.iloc[k, 0] + dist - 1



    # contig_name = re.split(">\.", temp[0])[-1]
    # if contig_name == "":
    #     contig_name = re.split(">", temp[0])[-1]
    #
    # contig_name = re.sub("\.1$","",contig_name)
    #
    # if contig_name[0] == ">":
    #     print(contig_name)
    #     contig_name = re.split("\.", contig_name)[-1]
    #
    # if contig_name[0] == ">":
    #     print(contig_name)
    #     contig_name = re.split("^>",contig_name)[1]

    contig_name = os.path.basename(files_for_input.fasta_file)
    contig_name = re.sub("\..*[a-z,A-Z]$", "", contig_name)
    contig_name = re.sub("#","_", contig_name)
    print(contig_name)


    file_path = "./" + contig_name + "#contig_bounds.csv"


    contig_bounds.to_csv(path_or_buf=file_path,
                         index=False)


