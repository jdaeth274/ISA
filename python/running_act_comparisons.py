import sys
import os
import pandas


purpose = ''' Script to take in a csv of all the locations of fastas for all isolates in first 
column and then the location of the reference gff for that particular strain. Usage:
python running_act_comparisons.py location_csv perl_script_loc act_compo_dir'''

input_args = sys.argv

if len(input_args) != 4:
    print("Not correct number of arguments, need 3, you have %s" % (len(input_args) -1 ))
    print(purpose)
    sys.exit()

input_names = ['isolate','reference']
input_csv = pandas.read_csv(input_args[1])
input_csv.columns = input_names
perl_dir = input_args[2]
act_dir = input_args[3]

## Get the unique references in the dataframe and work through those

unique_refs = input_csv.reference.unique()

for k in range(len(unique_refs)):
    ## lets narrow down the df to just for these references.
    current_data_set = input_csv[input_csv['reference'] == unique_refs[k]]
    reference_gff = str(current_data_set.iloc[0,1])
    print("On reference %s" % os.path.basename(reference_gff))
    print("")
    seqret_command = "seqret -sequence " + reference_gff + " -outseq referoo.fasta"
    os.system(seqret_command)
    for l in range(len(current_data_set.index)):
        print("On isolate number %s" % l, end="\r", flush=True)
        current_isolate = str(current_data_set.iloc[l, 0])
        perl_command = "perl " + perl_dir + "compare_genomes_orig.pl referoo.fasta " + current_isolate
        move_command = "mv *.crunch.gz " + act_dir
        os.system(perl_command)
        os.system(move_command)



