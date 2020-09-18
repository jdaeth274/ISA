import re
import pandas
import numpy
import sys
import dendropy
from dendropy.model.parsimony import fitch_down_pass
from dendropy.model.parsimony import fitch_up_pass
import csv
from statistics import stdev
from statistics import mean
import os
import string
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq
import math
import argparse




def get_options():
    purpose = '''This is a python script to intake a csv of hit locations for an MGE for a large collection of 
     genomes and then create a library of hits for which to search against 
    Usage: python library_creator.py <hit_csv> <reference_and_isolate_loc.csv> <align_cutoff> <path to act comparisons> <contig_checker_path>
    <out_csv>'''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='pen_checker_cdc.py')

    parser.add_argument('--hit_csv', required=True, help='Out merged BLAST csv (required)', type=str)
    parser.add_argument('--refrerence_csv', required=True, help='Isolates and references csv (required)', type=str)
    parser.add_argument('--align_cutoff', required=True, help='Alignment length cut_off for hit  (required)', type=int)
    parser.add_argument('--act_loc', required=True, help='Directory where act comparisons stored (required)', type=str)
    parser.add_argument('--contig_loc', required=True, help='Directory where contig_numbers stored  (required)', type=str)
    parser.add_argument('--output', required=True, help='Out Library csv name (required)', type=str)

    args = parser.parse_args()

    return args


def contig_end_checker(hitters, contig_tab, act_path, isolate_id):

    ## This function takes in the hitters and checks if near enough to the end of a
    ## contig to be considered a split up hit.
    last_occy = isolate_id.rfind('_')
    isolate_id_list = list(isolate_id)
    isolate_id_list[last_occy] = "#"
    isolate_id_2 = "".join(isolate_id_list)

    compo_file = act_path + isolate_id_2 + ".crunch.gz"
    compo_names = ['query', 'subject', 'pid', 'align', 'gap', 'mismatch', 'qstart',
                   'qend', 'sstart', 'send', 'eval', 'bitscore']
    compo_table = pandas.read_csv(compo_file, sep='\t', names=compo_names)

    min_hit_val = hitters.values.min()
    max_hit_val = hitters.values.max()

    if hitters.iloc[1] > hitters.iloc[0]:
        mge_ori = "forward"
    elif hitters.iloc[0] > hitters.iloc[1]:
        mge_ori = "reverse"



    contig = contig_tab[(contig_tab['0'] <= (min_hit_val + 15)) & (contig_tab['1'] >= (max_hit_val - 15))]


    start_pos = "None"
    end_pos = "None"

    ## Now we need to check if there are any act hits across this bound


    if len(contig.index) != 1:
        pos = "None"
    else:
        pos = "None"




        if (min_hit_val - 65) <= contig.iloc[0,0] <= (min_hit_val + 15) :
            if min_hit_val == hitters.iloc[0]:
                start_pos = "start"
            else:
                end_pos = "end"

            if  (max_hit_val - 15) <= contig.iloc[0,1] <= (max_hit_val + 50):
                if start_pos == "start":
                    end_pos = "end"
                else:
                    start_pos = "start"
        elif (max_hit_val - 15) <= contig.iloc[0,1] <= (max_hit_val + 50):
            if max_hit_val == hitters.iloc[0]:
                start_pos  = "start"
            else:
                end_pos = "end"

        if (start_pos == "None") or (end_pos == "None"):
            if mge_ori == "forward":
                before_region_to_check = [contig.iloc[0,0], hitters.iloc[0]]
                after_region_to_check = [hitters.iloc[1], contig.iloc[0,1]]
            elif mge_ori == "reverse":
                before_region_to_check = [hitters.iloc[ 0], contig.iloc[0, 1]]
                after_region_to_check = [contig.iloc[0, 0] , hitters.iloc[ 1]]

            before_spans = compo_table[(compo_table['qstart'] <= before_region_to_check[0]) & (compo_table['qend'] >= before_region_to_check[1])]
            after_spans = compo_table[(compo_table['qstart'] <= after_region_to_check[0]) & (compo_table['qend'] >= after_region_to_check[1])]

            before_within = compo_table[(compo_table['qstart'] >= before_region_to_check[0]) & (compo_table['qend'] <= before_region_to_check[1])]
            after_within = compo_table[(compo_table['qstart'] >= after_region_to_check[0]) & (compo_table['qend'] <= after_region_to_check[1])]

            before_overlap_1 = compo_table[(compo_table['qstart'] >= before_region_to_check[0]) & (compo_table['qstart'] <= before_region_to_check[1]) & (compo_table['qend'] >= before_region_to_check[1])]
            before_overlap_2 = compo_table[(compo_table['qstart'] <= before_region_to_check[0]) & (compo_table['qend'] <= before_region_to_check[1]) & (compo_table['qend'] >= before_region_to_check[0])]

            after_overlap_1 = compo_table[(compo_table['qstart'] >= after_region_to_check[0]) & (compo_table['qend'] >= after_region_to_check[1]) & (compo_table['qstart'] <= after_region_to_check[1])]
            after_overlap_2 =  compo_table[(compo_table['qstart'] <= after_region_to_check[0]) & (compo_table['qend'] <= after_region_to_check[1]) & (compo_table['qend'] >= after_region_to_check[0])]





            before_df = pandas.concat([before_spans, before_within, before_overlap_1, before_overlap_2])
            after_df = pandas.concat([after_spans, after_within, after_overlap_1, after_overlap_2])


            if start_pos == "None":
                if before_df.empty:
                    start_pos = "start"
                else:
                    high_pid_only = before_df[before_df['pid'] >= 98.0]
                    if high_pid_only.empty:
                        start_pos = "start"
                    else:
                        acceptable_contig_gap = round((before_region_to_check[1] - before_region_to_check[0]) * 0.1)
                        aligns = []

                        before_region_to_check_set = set(range(int(before_region_to_check[0]), int(before_region_to_check[1])))
                        for k in range(len(high_pid_only.index)):
                            current_range = range(high_pid_only.iloc[k, 6], high_pid_only.iloc[k, 7])
                            aligno = len(before_region_to_check_set.intersection(current_range))
                            aligns.append(aligno)

                        divisibleBySeven = [num for num in aligns if num > acceptable_contig_gap]
                        if len(divisibleBySeven) == 0:
                            start_pos = "start"
                        else:
                            acceptable_contig_gap = round((before_region_to_check[1] - before_region_to_check[0]) * 0.1)
                            high_pid_high_align = high_pid_only[high_pid_only['align'] >= acceptable_contig_gap]

                            if high_pid_high_align.empty:
                                start_pos = "start"

            if end_pos == "None":
                if after_df.empty:
                    end_pos = "end"
                else:
                    high_pid_only = after_df[after_df['pid'] >= 98.0]
                    high_pid_only.reset_index(drop=True)
                    if high_pid_only.empty:
                        end_pos = "end"
                    else:
                        acceptable_contig_gap = round((after_region_to_check[1] - after_region_to_check[0]) * 0.1)
                        aligns = []
                        after_region_to_check_set = set(range(int(after_region_to_check[0]), int(after_region_to_check[1])))
                        for k in range(len(high_pid_only.index)):
                            current_range = range(high_pid_only.iloc[k, 6], high_pid_only.iloc[k, 7])
                            aligno = len(after_region_to_check_set.intersection(current_range))
                            aligns.append(aligno)

                        divisibleBySeven = [num for num in aligns if num > acceptable_contig_gap]
                        if len(divisibleBySeven) == 0:
                            end_pos = "end"
                        else:
                            high_pid_high_align = high_pid_only[high_pid_only['align'] >= acceptable_contig_gap]
                            if high_pid_high_align.empty:
                                end_pos = "end"

    pos = start_pos + "_" + end_pos



    return  pos


def row_merger(narrowed_rows):
    ## This function takes in the narrowed rows from merged_contig_checker and then
    ## works out how best to merge them
    narrowed_rows['end_loc'] = narrowed_rows['contig_pos'].str.find("end")
    narrowed_rows['start_loc'] = narrowed_rows['contig_pos'].str.find("start")

    end_rows = narrowed_rows[narrowed_rows['end_loc'] != -1]


    if len(end_rows.index) == 0:
        # print(first_row.iloc[0])
        returned_row = narrowed_rows.drop(['merged_index', 'contig_pos', 'end_loc', 'start_loc'], axis=1)
        merged_indexers = []
        merged_locs = []
        return returned_row, merged_indexers, merged_locs
    else:
        start_rows = narrowed_rows[narrowed_rows['start_loc'] != -1]
        start_rows = start_rows.reset_index(drop=True)
        first_row = end_rows.iloc[0, :].copy()
        ## while loop to incorporate rest of hits
        if len(start_rows.index) == 0:
            #print(first_row.iloc[0])
            returned_row = narrowed_rows.drop(['merged_index', 'contig_pos', 'end_loc', 'start_loc'], axis = 1 )
            merged_indexers = []
            merged_locs = []
            return returned_row, merged_indexers, merged_locs
        else:



            first_start = start_rows.iloc[0,:]
            first_start_index = start_rows.index[0]
            if first_start['subject'] == first_row['subject']:
                first_start = start_rows.iloc[1, :]
                first_start_index = start_rows.index[1]

            counter = 0
            merged_row = first_row
            merged_locs = [[first_row['sstart'], first_row['ssend']]]
            merged_indexers = [merged_row.loc['merged_index']]
            added_start_indexes = []

            while counter < len(start_rows.index):

                if counter > 1000:
                    print("\n","Could be stuck in a while loop here")
                    print(narrowed_rows.iloc[0,0])
                    print(narrowed_rows)
                    sys.exit(1)

                current_range = merged_row.iloc[[3,4]]
                current_align = merged_row.iloc[2]
                poss_align = first_start.iloc[2]
                poss_range = first_start.iloc[[3,4]]
                overlap_test = current_range[0] <= poss_range[0] and poss_range[1] <= current_range[1] and \
                               current_align > poss_align
                if overlap_test == True:
                    counter += 1
                    added_start_indexes.append(first_start_index)
                    first_start_index += 1
                    if first_start_index <= (len(start_rows.index) - 1):
                        first_start = start_rows.iloc[first_start_index]


                else:
                    if first_start['contig_pos'] == "start_None":
                        ## check for other start_nones in the narrowed_df
                        start_nones = start_rows[start_rows['contig_pos'] == "start_None"]
                        if len(start_nones.index) == 1:
                                new_align = merged_row['align'] + start_nones.iloc[0,2]
                                new_qend = start_nones.iloc[0, 4]
                                new_send = start_nones.iloc[0, 6]
                                merged_row.loc['align'] = new_align
                                merged_row.loc['qend'] = new_qend
                                merged_row.loc['ssend'] = new_send
                                merged_locs.append([first_start.iloc[5],first_start.iloc[6]])
                                merged_indexers.append(first_start.iloc[10])
                                added_start_indexes.append(first_start_index)
                                first_start_index += 1
                                if first_start_index <= (len(start_rows.index) - 1):
                                    first_start = start_rows.iloc[first_start_index]

                                counter = len(start_rows.index)

                        else:
                                current_pid = merged_row['pident']
                                poss_sstart = start_nones['sstart'].values.tolist()
                                #closest_pid_val = min(poss_pids, key=lambda x: abs(x - current_pid))
                                #closest_pid_index = poss_pids.index(closest_pid_val)
                                lowest_sstart_index = poss_sstart.index(min(poss_sstart))

                                new_align = merged_row['align'] + start_nones.iloc[lowest_sstart_index, 2]
                                new_qend = start_nones.iloc[lowest_sstart_index, 4]
                                new_send = start_nones.iloc[lowest_sstart_index, 6]
                                merged_row.loc['align'] = new_align
                                merged_row.loc['qend'] = new_qend
                                merged_row.loc['ssend'] = new_send
                                merged_locs.append([start_nones.iloc[lowest_sstart_index,5], start_nones.iloc[lowest_sstart_index, 6]])
                                merged_indexers.append(start_nones.iloc[lowest_sstart_index, 10])
                                added_start_indexes.append(start_nones.index)
                                first_start_index += 1
                                if first_start_index <= (len(start_rows.index) - 1):
                                    first_start = start_rows.iloc[first_start_index]

                                counter += len(start_nones.index)
                    else:
                            new_align = merged_row['align'] + first_start.iloc[2]
                            new_qend = first_start.iloc[4]
                            new_send = first_start.iloc[6]
                            merged_row.loc['align'] = new_align
                            merged_row.loc['qend'] = new_qend
                            merged_row.loc['ssend'] = new_send
                            merged_locs.append([first_start.iloc[5], first_start.iloc[6]])
                            merged_indexers.append(first_start.iloc[10])
                            added_start_indexes.append(first_start_index)
                            first_start_index += 1
                            if first_start_index <= (len(start_rows.index) - 1):
                                first_start = start_rows.iloc[first_start_index]
                            counter += 1


            returned_row = merged_row.drop(['merged_index','contig_pos','end_loc','start_loc'])

            return returned_row, merged_indexers, merged_locs


def merged_contig_checker(merged_csv, contig_file_abs_path, act_path):
    multi_rows = []



    for k in range(len(merged_csv.index)):
        current_id = merged_csv.iloc[k, 0]#.values.to_string()
        underscore_count = current_id.count("_")
        if underscore_count > 1:
            multi_rows.append(k)


    multi_row_db = merged_csv.iloc[multi_rows].copy()
    #print(multi_row_db['file_loc'])
    file_locs = set(multi_row_db.iloc[:,8].copy().to_list())
    file_locs = list(file_locs)

    print(len(file_locs))

    merged_rows_to_drop = []
    merged_rows_to_add = pandas.DataFrame(columns=merged_csv.columns)
    merged_locs = []


    for k in range(len(file_locs)):
        isolate_id = file_locs[k]
        isolate_rows = multi_row_db[multi_row_db['file_loc'] == isolate_id].copy()
        isolate_rows.columns = isolate_rows.columns.str.strip()
        isolate_rows['merged_index'] = isolate_rows.index
        isolate_rows = isolate_rows.reset_index(drop=True)
        isolate_rows = isolate_rows.sort_values(by = 'qstart', ascending=True)
        contig_suffix = "#contig_bounds.csv"
        contig_isolate = re.sub("#", "_", isolate_id)
        if ".contigs_velvet.fa" in contig_isolate:
            contig_isolate = re.sub(".contigs_velvet.fa","",contig_isolate)
        contig_file_path = contig_file_abs_path + contig_isolate + contig_suffix

        contig_tab = pandas.read_csv(contig_file_path)

        positions = []

        for l in range(len(isolate_rows.index)):
            current_hitters = isolate_rows.iloc[l, [5,6]]
            position = contig_end_checker(current_hitters, contig_tab, act_path, contig_isolate)
            positions.append(position)

        orig_positions = positions

        isolate_rows['contig_pos'] = pandas.Series(positions, index=isolate_rows.index)



        isolate_rows_narrow = isolate_rows[isolate_rows['contig_pos'] != "None_None"].copy()




        if len(isolate_rows_narrow.index) < 10 and len(isolate_rows_narrow.index) > 1:
            isolate_rows_narrow = isolate_rows_narrow.sort_values(by = 'qstart', ascending=True)
            isolate_rows_narrow = isolate_rows_narrow.reset_index(drop=True)
            merged_row, merged_row_to_drop, merged_loc = row_merger(isolate_rows_narrow)
            if len(merged_row_to_drop) > 0:
                merged_rows_to_add = merged_rows_to_add.append(merged_row)
                merged_rows_to_drop.extend(merged_row_to_drop)
                merged_locs.append(merged_loc)


        iter_val = "{0:0=3d}".format((k + 1))
        print("Completed %s of %s rows. Just finished: %s" % (iter_val, len(file_locs), isolate_id),
              end="\r",flush=True)

    ## Now we'll remove the merged rows from the df
    merged_csv = merged_csv.drop(merged_rows_to_drop)
    merged_csv['merged_index'] = numpy.NAN
    ## Now lets add in the rows
    merged_rows_to_add['merged_index'] = range(len(merged_rows_to_add.index))
    merged_csv = merged_csv.append(merged_rows_to_add)

    return merged_csv, merged_locs



## Ok so first lets load up the merged BLAST CSV and narrow it down to just those
## over the threshold length.

if __name__ == '__main__':
    tab = str.maketrans("ACTG", "TGAC")

    files_for_input = get_options()

    hit_csv = pandas.read_csv(files_for_input.hit_csv)
    contig_file_abs_path = files_for_input.contig_loc
    absolute_act_path = files_for_input.act_loc


    merged_csv , merged_locs = merged_contig_checker(hit_csv, contig_file_abs_path, absolute_act_path)
    is_2k = merged_csv['align'] >= int(files_for_input.align_cutoff)

    proper_hits = merged_csv[is_2k]
    proper_hits = proper_hits.reset_index(drop=True)

    ## Now lets load up the csv with the isolate names and their reference location




















