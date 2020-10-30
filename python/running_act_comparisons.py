import sys
import os
import pandas
import subprocess
import argparse
import re

def get_options():
    purpose = ''' Script to take in a csv of all the locations of fastas for all isolates in first 
    column and then the location of the reference gff for that particular strain. Usage:
    python running_act_comparisons.py location_csv perl_script_loc act_compo_dir'''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='pen_checker_cdc.py')

    parser.add_argument('--csv', required=True, help='csv of isolate fasta loc and reference fasta ', type=str)
    parser.add_argument('--perl_dir', required=True, help='directory where compare_genomes_orig script is', type=str)
    parser.add_argument('--act_dir', required=True, help='Out directory to store act comparisons', type=str)

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    input_args = get_options()
    input_names = ['isolate','reference', 'cluster_name']
    input_csv = pandas.read_csv(input_args.csv)
    input_csv.columns = input_names
    perl_dir = input_args.perl_dir
    act_dir = input_args.act_dir

    ## Get the unique references in the dataframe and work through those

    unique_refs = input_csv.reference.unique()
    print(unique_refs)

    for k in range(len(unique_refs)):
        ## lets narrow down the df to just for these references.
        current_data_set = input_csv[input_csv['reference'] == unique_refs[k]]
        reference_gff = str(current_data_set.iloc[0,1])
        ref_base = os.path.basename(reference_gff)

        print("On reference %s" % ref_base)
        print("")

        if bool(re.search("\.fasta$", ref_base)) or bool(re.search("\.fa$", ref_base)):
            ## first make .dna then run through command
            dna_file = re.sub("\..*$","",ref_base)
            dna_file = re.sub("#","_",dna_file)
            dna_file = dna_file + ".dna"
            dna_command = "perl " + perl_dir + "converting_velvet_contigs_to_dna.pl " + reference_gff
            mv_command = "mv " + dna_file + " referoo.fasta"
            subprocess.call(dna_command, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
            subprocess.call(mv_command, shell=True)
        else:
            seqret_command = "seqret -sequence " + reference_gff + " -outseq referoo.fasta"
            subprocess.call(seqret_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        for l in range(len(current_data_set.index)):
            print("On isolate number %s" % l, end="\r", flush=True)
            current_isolate = str(current_data_set.iloc[l, 0])
            perl_command = "perl " + perl_dir + "compare_genomes_orig.pl referoo.fasta " + current_isolate
            move_command = "mv *.crunch.gz " + act_dir
            subprocess.call(perl_command, shell=True, stdout=subprocess.DEVNULL)
            subprocess.call(move_command, shell=True, stdout=subprocess.DEVNULL)



