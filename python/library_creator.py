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
    parser.add_argument('--reference_csv', required=True, help='Isolates gff and references gff csv (required)', type=str)
    parser.add_argument('--align_cutoff', required=True, help='Alignment length cut_off for hit  (required)', type=int)
    parser.add_argument('--act_loc', required=True, help='Directory where act comparisons stored (required)', type=str)
    parser.add_argument('--contig_loc', required=True, help='Directory where contig_numbers stored  (required)', type=str)
    parser.add_argument('--output', required=True, help='Out Library csv name (required)', type=str)

    args = parser.parse_args()

    return args

def contig_checker(contig_csv, hit_to_check):
    ## This is a function to check what contig a BLAST hit appears on
    ## Input: contig_csv: The CSV file containing contig start in a single column and contig end in another
    ##        hit_to_check: Start and end of a BLAST hit
    ## Output: Contig number of BLAST hit.
    hit_contig = 0
    for j in range(len(contig_csv.index)):
        c_start_and_end = contig_csv.iloc[j]
        if j == 0:
            overhang_before = 0
            overhang_after = 10
        else:
            overhang_before = 15
            overhang_after = 15
        start_hit = hit_to_check[0] >= (c_start_and_end[0] - overhang_before) and \
                    hit_to_check[0] <= (c_start_and_end[1] + overhang_after)
        end_hit = hit_to_check[1] >= (c_start_and_end[0] - overhang_before) and \
                  hit_to_check[1] <= (c_start_and_end[1] + overhang_after)

        if start_hit == True and end_hit == True:
            hit_contig = j + 1

    return hit_contig

def bounds_of_contig(contig_tab, contig_mge):
    ## Function to get bounds of a contig
    contig_bounds = contig_tab.iloc[contig_mge - 1]

    return contig_bounds



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

def ref_contains_hit(compo_table, hitters, mge_bounds, isolate_id):
    ## This function checks if there are any hits within the mge position
    ## where compo_table is the act comparison, hitters the position of
    ## the mge in the query and mge bounds the contig of the mge.

    if hitters[0] < hitters[1]:
        mge_ori = "forward"
    elif hitters[0] > hitters[1]:
        mge_ori = "reverse"

    if mge_ori == "forward":
        whole_overlap = compo_table[(compo_table['qstart'] < hitters[0]) & (compo_table['qend'] > hitters[1])]
        whole_overlap = whole_overlap.sort_values(by=['qstart'], ascending=False)
        if whole_overlap.empty:
            hits_bef = compo_table[(compo_table['qend'] > hitters[0] + 100) & (compo_table['qstart'] < hitters[0])]
            hits_bef = hits_bef.sort_values(by=['qstart'], ascending=False)
            if hits_bef.empty:
                hit_bef = "Not"
                hit_aft = "Not"
            else:
                hit_bef = hits_bef.iloc[0]

                hits_aft = compo_table[(compo_table['qstart'] < hitters[1] - 100) & (compo_table['qend'] > hitters[1])]
                hits_aft = hits_aft.sort_values(by=['qend'], ascending=True)
                if hits_aft.empty:
                    hit_aft = "Not"
                else:
                    hit_aft = hits_aft.iloc[0]
        else:
            whole_match = whole_overlap.iloc[0]
            hit_bef, hit_aft = whole_match_splitter(whole_match, hitters)

    elif mge_ori == "reverse":
        whole_overlap = compo_table[(compo_table['qstart'] < hitters[1]) & (compo_table['qend'] > hitters[0])]
        whole_overlap = whole_overlap.sort_values(by=['qstart'], ascending=False)
        if whole_overlap.empty:
            hits_bef = compo_table[(compo_table['qstart'] < hitters[0] - 100) & (compo_table['qend'] > hitters[0])]
            hits_bef = hits_bef.sort_values(by=['qend'], ascending=True)
            if hits_bef.empty:
                hit_bef = "Not"
                hit_aft = "Not"
            else:
                hit_bef = hits_bef.iloc[0]

                hits_aft = compo_table[(compo_table['qend'] > hitters[1] + 100) & (compo_table['qstart'] < hitters[1])]
                hits_aft = hits_aft.sort_values(by=['qstart'], ascending=False)
                if hits_aft.empty:
                    hit_aft = "Not"
                else:
                    hit_aft = hits_aft.iloc[0]
        else:
            whole_match = whole_overlap.iloc[0]
            hit_bef, hit_aft = whole_match_splitter(whole_match, hitters)

    if type(hit_bef) == str or type(hit_aft) == str:
        overlap = "No"
    else:
        overlap = "Yes"


    return hit_bef, hit_aft, overlap

def whole_match_splitter(match, mge_locs):
    ## This is a function to split a whole match if found to a before hit that
    ## lines up with the start of the mge and then the after hit with the end

    if mge_locs[0] < mge_locs[1]:
        mge_ori = "forward"
    elif mge_locs[0] > mge_locs[1]:
        mge_ori = "reverse"

    if match['sstart'] < match['send']:
        match_ori = "forward"
    elif match['sstart'] > match['send']:
        match_ori = "reverse"

    if match_ori == "forward":
        if mge_ori == "forward":
            length_bef = mge_locs[0] - match['qstart']
            hit_bef = pandas.Series({'query': match['query'],
                                     'subject': match['subject'],
                                     'pid': match['pid'],
                                     'align': length_bef,
                                     'gap': match['gap'],
                                     'mismatch': match['mismatch'],
                                     'qstart': match['qstart'],
                                     'qend': mge_locs[0],
                                     'sstart': match['sstart'],
                                     'send': (match['sstart'] + length_bef),
                                     'eval': match['eval'],
                                     'bitscore': match['bitscore']})

            length_aft = match['qend'] - mge_locs[1]
            hit_aft = pandas.Series({'query': match['query'],
                                     'subject': match['subject'],
                                     'pid': match['pid'],
                                     'align': length_aft,
                                     'gap': match['gap'],
                                     'mismatch': match['mismatch'],
                                     'qstart': mge_locs[1],
                                     'qend': match['qend'],
                                     'sstart': (match['send'] - length_aft),
                                     'send': match['send'],
                                     'eval': match['eval'],
                                     'bitscore': match['bitscore']})
        elif mge_ori == "reverse":
            length_bef = match['qend'] - mge_locs[0]
            hit_bef = pandas.Series({'query': match['query'],
                                     'subject': match['subject'],
                                     'pid': match['pid'],
                                     'align': length_bef,
                                     'gap': match['gap'],
                                     'mismatch': match['mismatch'],
                                     'qstart': mge_locs[0],
                                     'qend': match['qend'],
                                     'sstart': (match['send'] - length_bef),
                                     'send': match['send'],
                                     'eval': match['eval'],
                                     'bitscore': match['bitscore']})

            length_aft = mge_locs[1] - match['qstart']
            hit_aft = pandas.Series({'query': match['query'],
                                     'subject': match['subject'],
                                     'pid': match['pid'],
                                     'align': length_aft,
                                     'gap': match['gap'],
                                     'mismatch': match['mismatch'],
                                     'qstart': match['qstart'],
                                     'qend': mge_locs[1],
                                     'sstart': match['sstart'],
                                     'send': (match['sstart'] + length_aft),
                                     'eval': match['eval'],
                                     'bitscore': match['bitscore']})
    elif match_ori == "reverse":
        if mge_ori == "forward":
            length_bef = mge_locs[0] - match['qstart']
            hit_bef = pandas.Series({'query': match['query'],
                                     'subject': match['subject'],
                                     'pid': match['pid'],
                                     'align': length_bef,
                                     'gap': match['gap'],
                                     'mismatch': match['mismatch'],
                                     'qstart': match['qstart'],
                                     'qend': mge_locs[0],
                                     'sstart': match['sstart'],
                                     'send': (match['sstart'] - length_bef),
                                     'eval': match['eval'],
                                     'bitscore': match['bitscore']})

            length_aft = match['qend'] - mge_locs[1]
            hit_aft = pandas.Series({'query': match['query'],
                                     'subject': match['subject'],
                                     'pid': match['pid'],
                                     'align': length_aft,
                                     'gap': match['gap'],
                                     'mismatch': match['mismatch'],
                                     'qstart': mge_locs[1],
                                     'qend': match['qend'],
                                     'sstart': (match['send'] + length_aft),
                                     'send': match['send'],
                                     'eval': match['eval'],
                                     'bitscore': match['bitscore']})
        elif mge_ori == "reverse":
            length_bef = match['qend'] - mge_locs[0]
            hit_bef = pandas.Series({'query': match['query'],
                                     'subject': match['subject'],
                                     'pid': match['pid'],
                                     'align': length_bef,
                                     'gap': match['gap'],
                                     'mismatch': match['mismatch'],
                                     'qstart': mge_locs[0],
                                     'qend': match['qend'],
                                     'sstart': (match['send'] + length_bef),
                                     'send': match['send'],
                                     'eval': match['eval'],
                                     'bitscore': match['bitscore']})

            length_aft = mge_locs[1] - match['qstart']
            hit_aft = pandas.Series({'query': match['query'],
                                     'subject': match['subject'],
                                     'pid': match['pid'],
                                     'align': length_aft,
                                     'gap': match['gap'],
                                     'mismatch': match['mismatch'],
                                     'qstart': match['qstart'],
                                     'qend': mge_locs[1],
                                     'sstart': match['sstart'],
                                     'send': (match['sstart'] - length_aft),
                                     'eval': match['eval'],
                                     'bitscore': match['bitscore']})

    return hit_bef, hit_aft

def before_and_after_hits(hit_info, compo_table, contig_bounds):
    ## This function looks for matches around the MGE loc based on their
    ## position in the query sequence. First looks for hits with align
    ## greater than 2000, then if this isn't on contig looks for hits
    ## simply on the contig, no matter the size.
    overhang = 50
    if hit_info[0] < hit_info[1]:
        hits_before = compo_table.loc[compo_table[compo_table.columns[7]] < (hit_info[0] + overhang)]
        hits_before = hits_before.sort_values(by=['qend'], ascending=False)
        hits_before_1k = hits_before.loc[hits_before['align'] > 2000]
        ## Just check if this align > 1500 hit is on the same contig, if not we
        ## use the original hits_before ordered df



        if hits_before_1k.empty:
            hits_before_1k = hits_before

        if hits_before_1k.iloc[0, 6] < (contig_bounds[0] - 10):
            hits_before_1k = hits_before
        else:
            hits_before_1k = hits_before_1k

        ## Check now if this align > 1500 df is empty, if so again we use the
        ## original hits before df
        if hits_before_1k.empty:
            hits_before_1k = hits_before
        if hits_before_1k.empty:
            hit_before = pandas.DataFrame(data=numpy.zeros(shape=(1, 12)),
                                          columns=hits_before.columns.values)
        else:
            hit_before = hits_before_1k.iloc[0]

        ## Now we get the major hit before the insertion

        hits_after = compo_table.loc[compo_table[compo_table.columns[6]] > (hit_info[1] - overhang)]
        hits_after = hits_after.sort_values(by=['qstart'], ascending=True)
        hits_after_1k = hits_after.loc[hits_after['align'] > 2000]

        if hits_after_1k.empty:
            hits_after_1k = hits_after

        if hits_after_1k.iloc[0, 7] > (contig_bounds[1] + 10):
            hits_after_1k = hits_after
        else:
            hits_after_1k = hits_after_1k

        if hits_after_1k.empty:
            hits_after_1k = hits_after
        if hits_after_1k.empty:
            hit_after = pandas.DataFrame(data=numpy.zeros(shape=(1, 12)),
                                         columns=hits_after.columns.values)
        else:
            hit_after = hits_after_1k.iloc[0]


    elif hit_info[0] > hit_info[1]:

        hits_before = compo_table.loc[compo_table[compo_table.columns[6]] > (hit_info[0] - overhang)]
        hits_before = hits_before.sort_values(by=['qstart'], ascending=True)
        hits_before_1k = hits_before.loc[hits_before['align'] > 2000]

        if hits_before_1k.empty:
            hits_before_1k = hits_before


        if hits_before_1k.iloc[0, 7] > (contig_bounds[1] + 10):
            hits_before_1k = hits_before
        else:
            hits_before_1k = hits_before_1k

        if hits_before_1k.empty:
            hits_before_1k = hits_before
        if hits_before_1k.empty:
            hit_before = pandas.DataFrame(data=numpy.zeros(shape=(1, 12)),
                                          columns=hits_before.columns.values)
        else:
            hit_before = hits_before_1k.iloc[0]

        hits_after = compo_table.loc[compo_table[compo_table.columns[7]] < (hit_info[1] + overhang)]
        hits_after = hits_after.sort_values(by=['qend'], ascending=False)
        hits_after_1k = hits_after.loc[hits_after['align'] > 2000]

        if hits_after_1k.iloc[0, 6] < (contig_bounds[0] - 10):
            hits_after_1k = hits_after
        else:
            hits_after_1k = hits_after_1k

        if hits_after_1k.empty:
            hits_after_1k = hits_after
        if hits_after_1k.empty:
            hit_after = pandas.DataFrame(data=numpy.zeros(shape=(1, 12)),
                                         columns=hits_after.columns.values)
        else:
            hit_after = hits_after_1k.iloc[0]


    is_aft_within_bef = False
    is_bef_within_aft = False
    hits_to_use = "both"

    if hit_before.iloc[3] > hit_after.iloc[3]:
        is_aft_within_bef = within_a_hit(hit_before, hit_after)
    elif hit_before.iloc[3] < hit_after.iloc[3]:
        is_bef_within_aft = within_a_hit(hit_after, hit_before)

    if is_aft_within_bef == True:
        hits_to_use = "before"
    if is_bef_within_aft == True:
        hits_to_use = "after"


    return hit_before, hit_after, hits_to_use

def within_a_hit(big_hit, little_hit):
    ## This only works with two hits of the same orientation
    if big_hit['sstart'] < big_hit['send']:
        orientation = "forward"
    elif big_hit['sstart'] > big_hit['send']:
        orientation = "reverse"

    if orientation == "forward":
        output = little_hit['sstart'] > big_hit['sstart'] and little_hit['send'] < big_hit['send']
    elif orientation == "reverse":
        output = little_hit['send'] > big_hit['send'] and little_hit['sstart'] < big_hit['sstart']

    return output

def gff_finder(gff_csv, isolate_id):
    ## Function to get the location of an isolates gff file
    isolate_check = isolate_id + "\."
    isolate_rows = gff_csv['isolate'].str.contains(isolate_check)
    isolate_row_indy = isolate_rows.where(isolate_rows == True)
    isolate_row_indy = isolate_row_indy.index[isolate_row_indy == True].tolist()
    if len(isolate_row_indy) != 1:
        print(isolate_id)
        print(isolate_row_indy)
    isolate_loc = gff_csv.iloc[isolate_row_indy,0]

    return isolate_loc

def gff_to_dna(gff_csv, contig_csv, isolate_id):
    ## Function to turn the contigs based values of gff gene locs into a continuous DNA
    ## based values (1 .. seq length). Needs there to be at least one gene on a contig
    ## Input gff_csv: gff for a contiged assembly
    ##       contig_csv: contig_bound csv

    finding_izzy = re.sub("#","_",isolate_id)


    if "contig" in gff_csv.iloc[1,0]:
        ## go through the contig_csv and add on differences
        for k in  range(1 , len(contig_csv.index)):
            current_contig = k + 1
            num_zeros = 6 - len(str(current_contig))
            empty_str = ""
            contig_finder = "contig" + (empty_str.join((["0"] * num_zeros))) + str(current_contig)
            num_to_add = contig_csv.iloc[k,0]
            contig_rows = gff_csv['seqid'].str.contains(contig_finder)
            contig_rows_indy = contig_rows.where(contig_rows == True)
            contig_rows_indy = contig_rows_indy.index[contig_rows_indy == True].tolist()

            gff_csv.iloc[contig_rows_indy,3] = gff_csv.iloc[contig_rows_indy, 3] + num_to_add - 1
            gff_csv.iloc[contig_rows_indy, 4] =  gff_csv.iloc[contig_rows_indy, 4] + num_to_add - 1
    elif finding_izzy in gff_csv.iloc[1,0]:
        for k in  range(1 , len(contig_csv.index)):
            current_contig = k + 1
            contig_finder = finding_izzy + "\." + str(current_contig)
            num_to_add = contig_csv.iloc[k,0]
            contig_rows = gff_csv['seqid'].str.contains(contig_finder)
            contig_rows_indy = contig_rows.where(contig_rows == True)
            contig_rows_indy = contig_rows_indy.index[contig_rows_indy == True].tolist()


            gff_csv.iloc[contig_rows_indy,3] = gff_csv.iloc[contig_rows_indy, 3] + num_to_add - 1
            gff_csv.iloc[contig_rows_indy, 4] =  gff_csv.iloc[contig_rows_indy, 4] + num_to_add - 1

    return gff_csv


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

    isolate_ref_gff = pandas.read_csv(files_for_input.reference_csv)

    ## Set up the library df

    library_id = []
    library_mge_start = []
    library_mge_end = []
    library_insert_start = []
    library_insert_end = []
    library_length = []
    library_num_genes = []
    library_mge_length = []
    library_flank_genes = []
    librarY_flank_gene_length = []





    ## Now loop through the blast results ##

    for k in range(len(proper_hits.index)):
        ## First we'll get the hit locs in place for each of the hits
        current_row = proper_hits.iloc[[k]]
        hitters = (list(current_row.iloc[0, [5, 6]]))
        if hitters[0] < hitters[1]:
            mge_ori = "forward"
        elif hitters[0] > hitters[1]:
            mge_ori = "reverse"

        ## Now we'll go into the tab files and get the compo files
        isolate_id_z = proper_hits.iloc[k, 0]
        if isolate_id_z.count('_') == 2:
            last_occy = isolate_id_z.rfind('_')
            isolate_id = isolate_id_z[0:last_occy]
        else:
            isolate_id = isolate_id_z

        print(isolate_id)
        current_gff_loc = gff_finder(isolate_ref_gff, isolate_id)
        compo_file = absolute_act_path + isolate_id + ".crunch.gz"
        compo_names = ['query', 'subject', 'pid', 'align', 'gap', 'mismatch', 'qstart',
                       'qend', 'sstart', 'send', 'eval', 'bitscore']
        compo_table = pandas.read_csv(compo_file, sep='\t', names=compo_names)

        ## Now we've got and opened the tab blast file, we should look for hits
        ## before the start of the MGE in the host genome.


        ## So now we've got the hits before and after the insertion site, we need to check if they're on the same
        ## contig as each other.

        contig_suffix = "#contig_bounds.csv"
        contig_isolate = re.sub("#", "_", isolate_id)
        contig_file_path = contig_file_abs_path + contig_isolate + contig_suffix

        contig_tab = pandas.read_csv(contig_file_path)

        contig_mge = contig_checker(contig_tab, hitters)

        mge_bounds = bounds_of_contig(contig_tab, contig_mge)

        ###########################################################################
        ## Now we'll do a check on the reference to see if it has the same ########
        ## identical hits as the current isolate in question ######################
        ###########################################################################

        hit_before, hit_after, overlap = ref_contains_hit(compo_table, hitters, mge_bounds, isolate_id)

        if overlap == "No":

            hit_before, hit_after, which_hit = before_and_after_hits(hitters, compo_table, mge_bounds)
            if isolate_id == "15608_5#22":
                print(hit_before, hit_after)
        else:
            which_hit = "both"

        if hit_before[0] == 0:
            contig_before = None
        else:
            hit_before_loc = hit_before.iloc[[6, 7]]
            hit_before_length = abs(hit_before_loc[1] - hit_before_loc[0])
            contig_before = contig_checker(contig_tab, hit_before_loc)

        if hit_after[0] == 0:
            contig_after = None
        else:
            hit_after_loc = hit_after.iloc[[6, 7]]
            hit_after_length = abs(hit_after_loc[1] - hit_after_loc[0])
            contig_after = contig_checker(contig_tab, hit_after_loc)

        all_one_tig = contig_before == contig_mge and contig_mge == contig_after and which_hit == "both"
        all_one_tig_5k = hit_before_length >= 5000 and hit_after_length >= 5000 and all_one_tig



        if all_one_tig_5k:

            ## If the before and after hits to the isolate all line up on one contig we can assume this is
            ## a good enough hit to be used in library creation. Also needs to have at least 5k bp either side
            ## in this hit

            current_gff = pandas.read_csv(current_gff_loc.iloc[0], sep='\t',
                                       names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                             'attributes'],
                                       header=None)
            #current_gff = gff_to_dna(current_gff, contig_tab, isolate_id)


            if mge_ori == "forward":
                ## Lets get the element length
                ## The insertion length
                ## Number of genes in the inserton
                ## number genes & average length 500bp either side.
                current_mge_length = hitters[1] - hitters[0]
                current_insert_locs = [hit_before_loc[1], hit_after_loc[0]]
                current_insert_length = int(hit_after_loc[0]) - int(hit_before_loc[1])
                genes = current_gff[(current_gff['start'] >= current_insert_locs[0]) & (current_gff['end'] <= current_insert_locs[1])]

                gene_num = len(genes.index)
                print(gene_num)
                ## 500 bp regions.





























