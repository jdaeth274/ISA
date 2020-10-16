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
from Bio import SeqIO
from Bio.Seq import Seq
import math
import argparse
import time
import pyfastx
import subprocess

def get_options():
    purpose = '''This is a python script to intake a csv of hit locations for an MGE for a large collection of 
     genomes and then create a library of hits for which to search against 
    Usage: python library_creator.py <hit_csv> <reference_and_isolate_loc.csv> <align_cutoff> <path to act comparisons> <contig_checker_path>
    <out_csv>'''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='pen_checker_cdc.py')

    parser.add_argument('--lib_csv', required=True, help='Out merged BLAST csv (required)', type=str)
    parser.add_argument('--reference_csv', required=True, help='Isolates gff and references gff csv (required)', type=str)
    parser.add_argument('--blast_csv', required=True, help='Merged blast csv (required)', type=str)
    parser.add_argument('--fasta_csv', required=True, help= "isolates fasta and reference fasta csv", type = str)
    parser.add_argument('--act_loc', required=True, help='Directory where act comparisons stored (required)', type=str)
    parser.add_argument('--contig_loc', required=True, help='Directory where contig_numbers stored  (required)', type=str)
    parser.add_argument('--align', required=True, help='Int value to cutoff align for blast hits (required)', type=int)
    parser.add_argument('--output', required=True, help='Out Library csv name (required)', type=str)

    args = parser.parse_args()

    return args

def n50_calc(fasta_list, reference_id):
    ## Function to calculate the largest n50 to take as new comparison for ACT
    ## given that the reference has the mge within
    ## Input: fasta_list: The list of the locations of the fastas for a cluster
    ##        reference_id: The old reference

    n50s = []
    ids = []

    for k in fasta_list:
        current_isolate = k
        current_basename = os.path.basename(current_isolate)
        current_id = re.sub("\..*[a-z,A-Z]*.$","",current_basename)
        if current_id != reference_id:
            current_fasta = pyfastx.Fasta(current_isolate)
            current_n50 = current_fasta.nl(50)[0]
            n50s.append(current_n50)
            ids.append(current_isolate)


    max_n50 = n50s.index(max(n50s))

    new_act_iso = ids[max_n50]

    return new_act_iso


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





            before_df = pandas.concat([before_spans, before_within, before_overlap_1, before_overlap_2], sort=False)
            after_df = pandas.concat([after_spans, after_within, after_overlap_1, after_overlap_2], sort=False)


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

    ## First we work out which hits have both starts and ends
    narrowed_rows['end_loc'] = narrowed_rows['contig_pos'].str.find("end")
    narrowed_rows['start_loc'] = narrowed_rows['contig_pos'].str.find("start")

    ## End rows
    end_rows = narrowed_rows[narrowed_rows['end_loc'] != -1]


    if len(end_rows.index) == 0:
        ## If no end rows just take the returned row as the end
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
            ## If no start rows just take the returned row as the end
            #print(first_row.iloc[0])
            returned_row = narrowed_rows.drop(['merged_index', 'contig_pos', 'end_loc', 'start_loc'], axis = 1 )
            merged_indexers = []
            merged_locs = []
            return returned_row, merged_indexers, merged_locs
        else:

            ## Now we'll loop through the hits that have a start at the begininning of a contig and then
            ## altered the merged row with this new data. The overlap test will then check if the
            ## new start row to be potentially merged is within the current merged hit.

            first_start = start_rows.iloc[0,:]
            first_start_index = start_rows.index[0]
            if first_start['subject'] == first_row['subject']:
                first_start = start_rows.iloc[1, :]
                first_start_index = start_rows.index[1]

            counter = 0
            merged_row = first_row
            merged_locs = [[first_row['sstart'], first_row['ssend']]]
            merged_row['align'] = abs(merged_locs[0][1] - merged_locs[0][0])
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
                                new_align = merged_row['align'] + abs(start_nones.iloc[0,6] - start_nones.iloc[0,5])
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

                                new_align = merged_row['align'] + abs(start_nones.iloc[lowest_sstart_index,6] - start_nones.iloc[lowest_sstart_index,5])
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
                            new_align = merged_row['align'] + abs(first_start.iloc[6] - first_start.iloc[5])
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
            returned_row['merged'] = "Yes"

            return returned_row, merged_indexers, merged_locs


def merged_contig_checker(merged_csv, contig_file_abs_path, act_path):
    multi_rows = []

    for k in range(len(merged_csv.index)):
        current_id = merged_csv.iloc[k, 0]  # .values.to_string()
        underscore_count = current_id.count("_")
        if underscore_count > 1:
            current_loc = merged_csv['file_loc'].iloc[k]
            loccies = merged_csv[merged_csv['file_loc'] == current_loc]
            if len(loccies.index) > 1:
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
                merged_rows_to_add = merged_rows_to_add.append(merged_row, sort=False, ignore_index=True)
                merged_rows_to_drop.extend(merged_row_to_drop)
                merged_locs.append(merged_loc)


        iter_val = "{0:0=5d}".format((k + 1))
        print("Completed %s of %s rows. Just finished: %s" % (iter_val, len(file_locs), isolate_id),
              end="\r",flush=True)

    ## Now we'll remove the merged rows from the df
    merged_csv = merged_csv.drop(merged_rows_to_drop)
    merged_csv['merged_index'] = numpy.NAN
    merged_csv['merged'] = pandas.Series(list(numpy.repeat("No", len(merged_csv.index))), index=merged_csv.index)
    ## Now lets add in the rows
    merged_rows_to_add['merged_index'] = range(len(merged_rows_to_add.index))
    merged_csv = merged_csv.append(merged_rows_to_add, sort=False, ignore_index=True)

    return merged_csv, merged_locs

def ref_contains_hit(compo_table, hitters, mge_bounds, isolate_id):
    ## This function checks if there are any hits within the mge position
    ## where compo_table is the act comparison, hitters the position of
    ## the mge in the query and mge bounds the contig of the mge.

    if hitters[0] < hitters[1]:
        mge_ori = "forward"
    elif hitters[0] > hitters[1]:
        mge_ori = "reverse"

    cut_off_val = abs(hitters[1] - hitters[0]) * 0.5

    if mge_ori == "forward":
        whole_overlap = compo_table[(compo_table['qstart'] < hitters[0]) & (compo_table['qend'] > hitters[1])]
        whole_overlap = whole_overlap.sort_values(by=['qstart'], ascending=False)

        if whole_overlap.empty:
            ## Check if one hit is within the mge bounday above the cut off val established above

            middle_hits = compo_table[(compo_table['qstart'] > hitters[0]) &\
                                      (compo_table['qend'] < hitters[1])]
            if middle_hits.empty:
                middle_length = 0
            else:
                middle_length = middle_hits['align'].sum()

            ## Check if there are hits that overlap with the MGE before the MGE insertion.

            hits_bef = compo_table[(compo_table['qend'] > hitters[0])  & (compo_table['qstart'] < hitters[0])]
            hits_bef = hits_bef.sort_values(by=['qstart'], ascending=False)

            if hits_bef.empty:
                hit_bef_length = 0
                hit_bef = before_and_after_hits(hitters, compo_table, mge_bounds, "before")
                if hit_bef.iloc[0] == 0:
                    hit_bef = "No"
            else:
                hit_bef = hits_bef.iloc[0]
                hit_bef_length = hit_bef['qend'] - hitters[0]
                hit_bef = whole_match_splitter(hit_bef, hitters, "before")

            ## Check if there are hits that overlap with the MGE and the after flank region.

            hits_aft = compo_table[(compo_table['qstart'] < hitters[1] ) & (compo_table['qend'] > hitters[1])]
            hits_aft = hits_aft.sort_values(by=['qend'], ascending=True)
            if hits_aft.empty:
                hit_aft_length = 0
                hit_aft = before_and_after_hits(hitters, compo_table, mge_bounds, "after")
                if hit_aft.iloc[0] == 0:
                    hit_aft = "No"
            else:
                hit_aft = hits_aft.iloc[0]
                hit_aft_length = hitters[1] - hit_aft['qstart']
                hit_aft = whole_match_splitter(hit_aft, hitters, "after")


            tot_overlap_length = hit_bef_length + middle_length + hit_aft_length
            bef_aft_length = hit_bef_length + hit_aft_length

            ## If the overall overlap is below the threshold and the hit_afts and befs don't
            ## have a big enough overlap we can use the methods used for other hits.
            ## If however there is a big overlap (could come from a big middle hit and smaller afters & before
            ## or big afters and before) then lets use these new hits.
            ## if however there is a below overlap, but one (or both) of the after and before hits have an overlap
            ## greater than 50 (which the normal methods would miss) then use these new hits.
            ## Finally if there is a big overlap but neither of the before or after hits are overlapping (so there is a
            ## big middle overlap) just use the estimated hits.

            if tot_overlap_length < cut_off_val and hit_bef_length < 50 and hit_aft_length < 50:
                hit_bef = "No"
                hit_aft = "No"
            elif tot_overlap_length > cut_off_val and hit_bef_length > 0 and hit_aft_length > 0:
                hit_bef = hit_bef
                hit_aft = hit_aft
            elif tot_overlap_length < cut_off_val and (hit_bef_length > 50 or hit_aft_length > 50):
                hit_bef = hit_bef
                hit_aft = hit_aft
            elif tot_overlap_length > cut_off_val and hit_bef_length == 0 and hit_aft_length == 0:
                hit_bef = hit_bef
                hit_aft = hit_aft



        else:
            whole_match = whole_overlap.iloc[0]
            hit_bef, hit_aft = whole_match_splitter(whole_match, hitters, "both")

    elif mge_ori == "reverse":
        whole_overlap = compo_table[(compo_table['qstart'] < hitters[1]) & (compo_table['qend'] > hitters[0])]
        whole_overlap = whole_overlap.sort_values(by=['qstart'], ascending=False)
        if whole_overlap.empty:
            if whole_overlap.empty:
                ## Check if one hit is within the mge bounday above the cut off val established above

                middle_hits = compo_table[(compo_table['qstart'] >= hitters[1]) & \
                                          (compo_table['qend'] <= hitters[0])]
                if middle_hits.empty:
                    middle_length = 0
                else:
                    middle_length = middle_hits['align'].sum()

                hits_bef = compo_table[(compo_table['qstart'] < hitters[0]) & (compo_table['qend'] > hitters[0])]
                hits_bef = hits_bef.sort_values(by=['qstart'], ascending=True)

                if hits_bef.empty:
                    hit_bef_length = 0
                    hit_bef = before_and_after_hits(hitters, compo_table, mge_bounds, "before")
                    if hit_bef.iloc[0] == 0:
                        hit_bef = "No"
                else:
                    hit_bef = hits_bef.iloc[0]
                    hit_bef_length = hitters[0] - hit_bef['qstart']
                    hit_bef = whole_match_splitter(hit_bef, hitters, "before")

                hits_aft = compo_table[(compo_table['qend'] > hitters[1]) & (compo_table['qstart'] < hitters[1])]
                hits_aft = hits_aft.sort_values(by=['qend'], ascending=False)
                if hits_aft.empty:
                    hit_aft_length = 0
                    hit_aft = before_and_after_hits(hitters, compo_table, mge_bounds, "after")
                    if hit_aft.iloc[0] == 0:
                        hit_aft = "No"
                else:
                    hit_aft = hits_aft.iloc[0]
                    hit_aft_length = hit_aft['qend'] - hitters[1]
                    hit_aft = whole_match_splitter(hit_aft, hitters, "after")

                tot_overlap_length = hit_bef_length + middle_length + hit_aft_length
                bef_aft_length = hit_bef_length + hit_aft_length

                if tot_overlap_length < cut_off_val and hit_bef_length < 50 and hit_aft_length < 50:
                    hit_bef = "No"
                    hit_aft = "No"
                elif tot_overlap_length > cut_off_val and hit_bef_length > 0 and hit_aft_length > 0:
                    hit_bef = hit_bef
                    hit_aft = hit_aft
                elif tot_overlap_length < cut_off_val and (hit_bef_length > 50 or hit_aft_length > 50):
                    hit_bef = hit_bef
                    hit_aft = hit_aft

        else:
            whole_match = whole_overlap.iloc[0]
            hit_bef, hit_aft = whole_match_splitter(whole_match, hitters, "both")

    if type(hit_bef) == str or type(hit_aft) == str:
        overlap = "No"
    else:
        overlap = "Yes"


    return hit_bef, hit_aft, overlap

def whole_match_splitter(match, mge_locs, hit_to_split):
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
            if hit_to_split == "both" or hit_to_split == "before":

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
            if hit_to_split == "both" or hit_to_split == "after":
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
            if hit_to_split == "both" or hit_to_split == "before":
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
            if hit_to_split == "both" or hit_to_split == "after":
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
            if hit_to_split == "both" or hit_to_split == "before":
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
            if hit_to_split == "both" or hit_to_split == "after":
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
            if hit_to_split == "both" or hit_to_split == "before":
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
            if hit_to_split == "both" or hit_to_split == "after":
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

    if hit_to_split == "both":
        return hit_bef, hit_aft
    elif hit_to_split == "after":
        return hit_aft
    elif hit_to_split == "before":
        return  hit_bef

def before_and_after_hits(hit_info, compo_table, contig_bounds, hits_to_search):
    ## This function looks for matches around the MGE loc based on their
    ## position in the query sequence. First looks for hits with align
    ## greater than 2000, then if this isn't on contig looks for hits
    ## simply on the contig, no matter the size.
    overhang = 50

    if hit_info[0] < hit_info[1]:
        if hits_to_search == "both" or hits_to_search == "before":
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
        if hits_to_search == "both" or hits_to_search == "after":
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
        if hits_to_search == "both" or hits_to_search == "before":
            hits_before = compo_table.loc[compo_table[compo_table.columns[6]] > (hit_info[0] - overhang)]
            hits_before = hits_before.sort_values(by=['qstart'], ascending=True)
            hits_before_1k = hits_before.loc[hits_before['align'] > 2000]

            if hits_before_1k.empty:
                hits_before_1k = hits_befores


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
        if hits_to_search == "both" or hits_to_search == "after":
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

    hits_to_use = "both"
    if hits_to_search == "both":
        is_aft_within_bef = False
        is_bef_within_aft = False


        if hit_before.iloc[3] > hit_after.iloc[3]:
            is_aft_within_bef = within_a_hit(hit_before, hit_after)
        elif hit_before.iloc[3] < hit_after.iloc[3]:
            is_bef_within_aft = within_a_hit(hit_after, hit_before)

        if is_aft_within_bef == True:
            hits_to_use = "before"
        if is_bef_within_aft == True:
            hits_to_use = "after"

    if hits_to_search == "both":
        return hit_before, hit_after, hits_to_use
    elif hits_to_search == "before":
        return  hit_before
    elif hits_to_search == "after":
        return  hit_after

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

def gff_finder(gff_csv, isolate_id, clus_name):
    ## Function to get the location of an isolates gff file
    isolate_check = isolate_id + "\."
    isolate_rows = gff_csv['isolate'].str.contains(isolate_check)
    isolate_row_indy = isolate_rows.where(isolate_rows == True)
    isolate_row_indy = isolate_row_indy.index[isolate_row_indy == True].tolist()
    if len(isolate_row_indy) != 1:
        print(isolate_id)
        print(isolate_row_indy)
    isolate_loc = gff_csv.iloc[isolate_row_indy,0]
    isolate_ref = gff_csv.iloc[isolate_row_indy,1]
    if clus_name:
        cluster_name = gff_csv.iloc[isolate_row_indy,2]
    else:
        cluster_name = "BOO!"

    return isolate_loc, isolate_ref, cluster_name

def gff_to_dna(gff_csv, contig_csv, isolate_id, input_k):
    ## Function to turn the contigs based values of gff gene locs into a continuous DNA
    ## based values (1 .. seq length). Needs there to be at least one gene on a contig
    ## Input gff_csv: gff for a contiged assembly
    ##       contig_csv: contig_bound csv

    finding_izzy = re.sub("#","_",isolate_id)
    contig_lengths = contig_csv.iloc[:,1] - contig_csv.iloc[:,0]
    index_to_keep = contig_lengths[contig_lengths > 10000].index
    ## remove DNA sequences ##
    starting_seq = gff_csv['seqid'].str.startswith("##FASTA", na = False)
    starting_seq_indy = starting_seq.where(starting_seq == True)
    starting_seq_indy = starting_seq_indy.index[starting_seq_indy == True].tolist()
    starting_seq_indy = min(starting_seq_indy)

    narrowed_gff = gff_csv.drop(range(starting_seq_indy,(len(gff_csv.index) - 1) ))
    narrowed_gff = narrowed_gff.reset_index(drop=True)


    if "contig" in narrowed_gff.iloc[1,0]:
        ## go through the contig_csv and add on differences
        for k in index_to_keep :
            if k == 0:
                current_contig = k + 1
                num_zeros = 6 - len(str(current_contig))
                empty_str = ""
                contig_finder = "contig" + (empty_str.join((["0"] * num_zeros))) + str(current_contig) + "$"
                contig_rows = narrowed_gff['seqid'].str.contains(contig_finder)
                contig_rows_indy = contig_rows.where(contig_rows == True)
                contig_rows_indy = contig_rows_indy.index[contig_rows_indy == True].tolist()
                last_indy = max(contig_rows_indy)

            else:
                current_contig = k + 1
                num_zeros = 6 - len(str(current_contig))
                empty_str = ""
                contig_finder = "contig" + (empty_str.join((["0"] * num_zeros))) + str(current_contig) + "$"
                num_to_add = contig_csv.iloc[k,0]
                contig_rows = narrowed_gff['seqid'].str.contains(contig_finder)
                contig_rows_indy = contig_rows.where(contig_rows == True)
                contig_rows_indy = contig_rows_indy.index[contig_rows_indy == True].tolist()
                if len(contig_rows_indy) == 0:
                    continue

                narrowed_gff.iloc[contig_rows_indy,3] = narrowed_gff.iloc[contig_rows_indy, 3] + num_to_add - 1
                narrowed_gff.iloc[contig_rows_indy, 4] =  narrowed_gff.iloc[contig_rows_indy, 4] + num_to_add - 1

                last_indy = max(contig_rows_indy)

        gff_index_to_drop = range((last_indy + 1), (len(narrowed_gff.index) - 1))
        out_gff_csv = narrowed_gff.drop(gff_index_to_drop)



    elif finding_izzy in narrowed_gff.iloc[1,0]:
        for k in  index_to_keep:
            if k == 0:
                current_contig = k + 1
                contig_finder = finding_izzy + "\." + str(current_contig) + "$"
                contig_rows = narrowed_gff['seqid'].str.contains(contig_finder)
                contig_rows_indy = contig_rows.where(contig_rows == True)
                contig_rows_indy = contig_rows_indy.index[contig_rows_indy == True].tolist()
                last_indy = max(contig_rows_indy)

            else:
                current_contig = k + 1
                contig_finder = finding_izzy + "\." + str(current_contig) + "$"
                num_to_add = contig_csv.iloc[k,0]


                contig_rows = narrowed_gff['seqid'].str.contains(contig_finder)
                contig_rows_indy = contig_rows.where(contig_rows == True)
                contig_rows_indy = contig_rows_indy.index[contig_rows_indy == True].tolist()
                if len(contig_rows_indy) == 0:
                    continue


                narrowed_gff.iloc[contig_rows_indy,3] = narrowed_gff.iloc[contig_rows_indy, 3] + num_to_add - 1
                narrowed_gff.iloc[contig_rows_indy, 4] =  narrowed_gff.iloc[contig_rows_indy, 4] + num_to_add - 1

                last_indy = max(contig_rows_indy)

        gff_index_to_drop = range((last_indy + 1), (len(narrowed_gff.index) ))
        out_gff_csv = narrowed_gff.drop(gff_index_to_drop)

    ## test out finder ##


    out_gff_csv = out_gff_csv.reset_index(drop=True)

    return out_gff_csv

def gene_name_tryer(prospective_csv, library_csv, out_hit, missing_isolate, mergio):
    ## Function to take in hits they have missed the first few cutoffs in the hit_integrator function
    ## And try to merge them instead using the flank gene names and then the mge charecteristics
    ## Input: prospective_csv: The single line isolate from the hit_integrator function
    ##        library_csv: The library csv
    ##        out_hit: The total hit df
    ##        missing_isolate: The df for the isolates with not library hit

    missing_copy = missing_isolate.copy()
    out_hit_copy = out_hit.copy()

    bef_and_aft = library_csv[(library_csv['before_gene_name'] == prospective_csv['before_gene_name'][0]) & \
                              (library_csv['after_gene_name'] == prospective_csv['after_gene_name'][0])]

    if not bef_and_aft.empty:

        current_hit_length = bef_and_aft[(bef_and_aft['mge_length'] >= (prospective_csv['mge_length'][0] - 100)) & \
                                         (bef_and_aft['mge_length'] <= (prospective_csv['mge_length'][0] + 100))]

        if mergio:
            tot_hit = current_hit_length
        else:
            current_hit_gene = bef_and_aft[(bef_and_aft['mge_genes'] >= (prospective_csv['mge_genes'][0] - 1)) & \
                                       (bef_and_aft['mge_genes'] <= (prospective_csv['mge_genes'][0] + 1))]

            tot_hit = pandas.concat([current_hit_length, current_hit_gene], sort=False, ignore_index=True)
            tot_hit = tot_hit.drop_duplicates()


        if len(tot_hit.index) == 1:
            prospective_csv['insert_name'] = pandas.Series(tot_hit['insert_name'], index=prospective_csv.index)
            out_hit_copy = out_hit_copy.append(prospective_csv, sort = False)
        elif len(tot_hit.index) > 1:
            tot_hit = tot_hit.reset_index(drop=True)
            length_target = prospective_csv['mge_length'][0]
            mge_dist = abs(tot_hit['mge_length'] - length_target)
            min_val = mge_dist.idxmin()
            single_hit = tot_hit.iloc[min_val].to_frame().transpose()
            prospective_csv['insert_name'] = pandas.Series(single_hit['insert_name'], index=prospective_csv.index)
            out_hit_copy = out_hit_copy.append(prospective_csv, sort = False)
        else:
            missing_current = pandas.DataFrame()
            missing_current['id'] = pandas.Series(prospective_csv['id'])
            missing_current['mge_start'] = pandas.Series(prospective_csv['mge_start'], index=missing_current.index)
            missing_current['mge_end'] = pandas.Series(prospective_csv['mge_end'], index=missing_current.index)
            missing_current['ref_name'] = pandas.Series(prospective_csv['ref_name'], index=missing_current.index)
            missing_current['cluster_name'] = pandas.Series(prospective_csv['cluster_name'], index=missing_current.index)
            missing_current['mge_length'] = pandas.Series(prospective_csv['mge_length'], index=missing_current.index)
            missing_current['reason'] = pandas.Series(["No gene name matches"], index=missing_current.index)
            missing_copy = missing_copy.append(missing_current, sort = False)
    else:
        missing_current = pandas.DataFrame()
        missing_current['id'] = pandas.Series(prospective_csv['id'])
        missing_current['mge_start'] = pandas.Series(prospective_csv['mge_start'], index=missing_current.index)
        missing_current['mge_end'] = pandas.Series(prospective_csv['mge_end'], index=missing_current.index)
        missing_current['ref_name'] = pandas.Series(prospective_csv['ref_name'], index=missing_current.index)
        missing_current['cluster_name'] = pandas.Series(prospective_csv['cluster_name'], index=missing_current.index)
        missing_current['mge_length'] = pandas.Series(prospective_csv['mge_length'], index=missing_current.index)
        missing_current['reason'] = pandas.Series(["No gene name matches"], index=missing_current.index)
        missing_copy = missing_copy.append(missing_current, sort = False)

    return out_hit_copy, missing_copy

def hit_detector(library_csv, prospective_csv, isolate_id, hit_csv, missing_isolates, mergio):
    ## Function to allocate hit location
    ## Input: library_csv: The current library csv to search for matches against
    ##        prospective_csv: The prospective set of results to be sorted
    ##        isolate_id: The current isolate id (mainly for debugging)
    ##        hit_csv: The total hit csv the isolate to be integrated into
    ##        missing_isolates: The df of isolates not assigned to a hit value
    ##        mergio: Whether or not the isolate was merged together.
    ## This will first work on the basis of hits sharing almost identical blast matches
    ## Then we'll look into flank region composition and finally the total insert composition
    ## to see if this is a novel hit or not.

    lib_new = library_csv.copy()
    out_hit = hit_csv.copy()
    missing_df = missing_isolates.copy()
    ids_to_drop = []

    ## lets narrow down first by number of genes in the element if this is not exact I think its
    ## fair to assume this would be novel.
    #'before_flank_gene', 'after_flank_gene', 'before_flank_avg',
    #'after_flank_avg'




    ## So there is a hit with the same number of mge genes, let now check the element length +- 2 bp
    ## Needs to be a hit with no length +- 2bp and no genes +- 1

    if mergio:
        out_hit, missing_df = gene_name_tryer(prospective_csv, library_csv, out_hit, missing_df)
    else:


        mge_hits = lib_new[(lib_new['mge_genes'] >= (prospective_csv['mge_genes'][0] - 1)) &\
                           (lib_new['mge_genes'] <= (prospective_csv['mge_genes'][0] + 1))]
        mge_length_hits = mge_hits[(mge_hits['mge_length'] >= (prospective_csv['mge_length'][0] - 2))\
                                     & (mge_hits['mge_length'] <= (prospective_csv['mge_length'][0] + 2))]


        if not mge_length_hits.empty:

                ## Lets check if the before or after hits match quite closely with the number of genes +- 1
                ## Before
                before_gene_num = prospective_csv['before_flank_gene'][0]
                after_gene_num = prospective_csv['after_flank_gene'][0]

                before_gene_hits = mge_length_hits[(mge_length_hits['before_flank_gene'] >= (before_gene_num - 1)) & \
                                                   (mge_length_hits['before_flank_gene'] <= (before_gene_num + 1))]

                after_gene_hits = mge_length_hits[(mge_length_hits['after_flank_gene'] >= (after_gene_num - 1)) & \
                                                  (mge_length_hits['after_flank_gene'] <= (after_gene_num + 1))]

                before_empty = before_gene_hits.empty
                after_empty = after_gene_hits.empty

                if not before_empty or not after_empty:
                    ## So either the before or the after hit or both match to ones already in the database
                    ## check matches +- 25 bp
                    before_flank_empty = True
                    after_flank_empty = True



                    if not before_empty and after_empty:
                        before_mean_flank = prospective_csv['before_flank_avg'][0]

                        before_flank_means = before_gene_hits[(before_gene_hits['before_flank_avg'] >= (before_mean_flank - 25)) &\
                                                              (before_gene_hits['before_flank_avg'] <= (before_mean_flank + 25))]

                        before_flank_empty = before_flank_means.empty

                    elif before_empty and not after_empty:
                        after_mean_flank = prospective_csv['after_flank_avg'][0]

                        after_flank_means = after_gene_hits[(after_gene_hits['before_flank_avg'] >= (after_mean_flank - 25)) & \
                                                             (after_gene_hits['before_flank_avg'] <= (after_mean_flank + 25))]

                        after_flank_empty = after_flank_means.empty
                    elif not before_empty and not after_empty:
                        before_mean_flank = prospective_csv['before_flank_avg'][0]

                        before_flank_means = before_gene_hits[
                            (before_gene_hits['before_flank_avg'] >= (before_mean_flank - 25)) & \
                            (before_gene_hits['before_flank_avg'] <= (before_mean_flank + 25))]

                        before_flank_empty = before_flank_means.empty

                        after_mean_flank = prospective_csv['after_flank_avg'][0]

                        after_flank_means = after_gene_hits[
                            (after_gene_hits['after_flank_avg'] >= (after_mean_flank - 25)) & \
                            (after_gene_hits['after_flank_avg'] <= (after_mean_flank + 25))]

                        after_flank_empty = after_flank_means.empty

                    if not before_flank_empty or not after_flank_empty:

                        ## So the before or after (or both) flanks seem to match in composition
                        ## Now lets check if the insert length is similar, if not likely a novel insertion in the
                        ## same insert site as before.


                        insert_genes = prospective_csv['insert_genes'][0]
                        insert_length = prospective_csv['insert_length'][0]



                        if "before_flank_means" in locals() and "after_flank_means" in locals():

                            remaining_hits = pandas.concat([before_flank_means, after_flank_means], ignore_index=True, sort=False)
                            remaining_hits = remaining_hits.drop_duplicates()
                        elif "before_flank_means" in locals() and "after_flank_means" not in locals():
                            remaining_hits = before_flank_means
                        else:
                            remaining_hits = after_flank_means

                        gene_hits = remaining_hits[(remaining_hits['insert_genes'] >= (insert_genes - 2)) & \
                                                   (remaining_hits['insert_genes'] <= (insert_genes + 2))]

                        length_hits = remaining_hits[(remaining_hits['insert_length'] >= (insert_length - 500)) & \
                                                     (remaining_hits['insert_length'] <= (insert_length + 500))]



                        if gene_hits.empty and length_hits.empty:
                            out_hit, missing_df = gene_name_tryer(prospective_csv, library_csv, out_hit, missing_df, mergio)
                        else:
                            remain_48 = pandas.concat([gene_hits, length_hits], ignore_index=True, sort = False)
                            remain_48 = remain_48.drop_duplicates()


                            if len(remain_48.index) == 1:
                                prospective_csv['insert_name'] = pandas.Series(remain_48['insert_name'], index=prospective_csv.index)
                                out_hit = out_hit.append(prospective_csv, sort = False)
                            elif len(remain_48.index) > 1:
                                ## If there are two hits with similar tendencies base on insert length
                                remain_48 = remain_48.reset_index(drop=True)
                                target_insert = prospective_csv['insert_length'][0]
                                lengers = abs(remain_48['insert_length'] - target_insert)
                                closest_index = lengers.idxmin()
                                prospective_csv['insert_name'] = pandas.Series(remain_48['insert_name'].iloc[closest_index], index = remain_48.index)
                                out_hit = out_hit.append(prospective_csv, sort = False)
                            else:
                                print("Odd behaviour here")

                    else:
                        out_hit, missing_df = gene_name_tryer(prospective_csv, library_csv, out_hit, missing_df, mergio)
                else:
                    out_hit, missing_df = gene_name_tryer(prospective_csv, library_csv, out_hit, missing_df, mergio)
        else:
            ## Try to see if there is a hit with the same start and end gene names and mge_length +- 50 bp or gene +- 1
            out_hit, missing_df = gene_name_tryer(prospective_csv, library_csv, out_hit, missing_df, mergio)




    return(out_hit, missing_df)

def gene_name_finder(flanks_csv, back_it_up):
    ## Function to find the first row in a flanks region that has a gene name and
    ## return that gene name. Removes the paralog suffix attached.
    ## Input: flanks_csv =  Either before or after 5k region from library pros
    ##        back_it_up = Whether to work backwards throgh the csv to get the closest
    ##                  hits

    new_flanks = flanks_csv.reset_index(drop=True)
    search_res = ""

    if back_it_up:
        for k in reversed(range(len(new_flanks.index))):
            search_res = re.findall('gene=.*?;', new_flanks.iloc[k, 8])
            if search_res == []:
                continue
            else:
                search_res = search_res[0]
                search_res = re.findall('=.*?;', search_res)
                search_res = search_res[0]

                search_res = re.sub("=","",search_res)
                search_res = re.sub(";","",search_res)
                #search_res = search_res[3:-3]

                break
    else:
        for k in range(len(new_flanks.index)):
            search_res = re.findall('gene=.*?;', new_flanks.iloc[k,8])
            if search_res == []:
                continue
            else:
                search_res = search_res[0]
                search_res = re.findall('=.*?;', search_res)
                search_res = search_res[0]

                search_res = re.sub("=", "", search_res)
                search_res = re.sub(";", "", search_res)

                break

    if search_res == []:
        search_res = "NONE"

    if re.search("_[0-9]$",search_res):
        search_res = re.sub("_[0-9]$", "", search_res)


    return  search_res

def library_trim(library_csv):
    ## Function to run through the library csv and removes any that look like they could be duplicates
    ## This will be based solely on the before and after genes

    before_after_combo = library_csv['before_gene_name'].str.cat(library_csv['after_gene_name'], sep = "£")
    before_after_unique = before_after_combo.unique()
    indies_to_remove = []
    for k in range(len(before_after_unique)):
        combo = before_after_unique[k]
        before = re.split("£", combo, 2)[0]
        after = re.split("£", combo, 2)[1]

        bef_and_aft = library_csv[(library_csv['before_gene_name'] == before) &\
                                  (library_csv['after_gene_name'] == after)]

        if len(bef_and_aft.index) == 1:
            continue
        else:
            ## Work through each hit to see if they match closely to the others and then
            ## remove them.
            ids_to_drop = []
            for hit in bef_and_aft.index:

                current_hit = bef_and_aft.loc[hit]

                if current_hit['id'] in ids_to_drop:
                    continue
                ## hits based on insert length
                current_hit_length = bef_and_aft[(bef_and_aft['insert_length'] >= (current_hit['insert_length'] - 100)) &\
                                                 (bef_and_aft['insert_length'] <= (current_hit['insert_length'] + 100))]
                current_hit_gene = bef_and_aft[(bef_and_aft['insert_genes'] >= (current_hit['insert_genes'] - 1)) &\
                                               (bef_and_aft['insert_genes'] <= (current_hit['insert_genes'] + 1))]

                tot_hit = pandas.concat([current_hit_length, current_hit_gene], sort=False, ignore_index=True)
                tot_hit = tot_hit.drop_duplicates()



                if tot_hit.empty:
                    continue
                else:

                    ## choose max flank lengths as hit

                    max_row = pandas.to_numeric(tot_hit['flanks_length']).idxmax()
                    remove_rows = tot_hit.drop(max_row)
                    ids_to_lose = remove_rows['id'].tolist()
                    ids_to_drop.extend(ids_to_lose)

            ids = list(library_csv['id'][library_csv['id'].isin(ids_to_drop)].index)
            indies_to_remove.extend(ids)


    out_lib = library_csv.drop(indies_to_remove)
    print("Removing this many duplicates from library: %s, now %s hits" % (len(indies_to_remove), len(out_lib.index)))
    return  out_lib

def subject_checker(hit, ori, compo_csv, mge_bounds, mge_ori):
    ## Function to take in hits that don't have close hits in the query strand
    ## from the hit_mover function and then to check if these hits are reasonably
    ## close in the subject strand to allow for merging.
    ## Input: hit: the initial hit input into hit_mover
    ##        ori: The location of the hit before or after
    ##        compo_csv: The poss_hits csv that is narrowed to those with hits on the same contig and proximate to hit
    ##        mge_bounds: The contig bounds of the contig the mge is found on
    ##        mge_ori: The orientation of the MGE in question.

    if hit['sstart'] < hit['send']:
        hit_subject_ori = "forward"
    elif hit['send'] < hit['sstart']:
        hit_subject_ori = "reverse"

        #sys.exit()

    hit_new = hit.copy()

    if ori == "before":
        if mge_ori == "forward":
            if hit_subject_ori == "forward":
                narrowed_poss_hits = compo_csv[(compo_csv['send'] >= (hit['sstart'] - 1000)) & \
                                               (compo_csv['send'] <= (hit['sstart'] + 50))]
                narrowed_poss_hits = narrowed_poss_hits.sort_values(by = ['send'], ascending = False)
                if not narrowed_poss_hits.empty:
                    new_align = hit['qend'] - narrowed_poss_hits['qstart'].iloc[0]
                    if new_align >= 5000:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        altered = "merge"
                    else:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        narrowed_poss_hits_2 = compo_csv[(compo_csv['send'] >= (narrowed_poss_hits['sstart'].iloc[0] - 1000)) & \
                                                         (compo_csv['send'] <= (narrowed_poss_hits['sstart'].iloc[0] + 50))]
                        narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        while not narrowed_poss_hits_2.empty and new_align < 5000:
                            new_align = abs(hit['qend'] - narrowed_poss_hits_2['qstart'].iloc[0])
                            extra_hit = pandas.DataFrame(narrowed_poss_hits_2.iloc[0]).transpose()

                            hit_new = pandas.concat([hit_new, extra_hit.iloc[[0]]], sort=False, ignore_index=False)
                            hit_new = hit_new.reset_index(drop=True)
                            if new_align < 5000:
                                narrowed_poss_hits_2 = compo_csv[(compo_csv['send'] >= (narrowed_poss_hits['sstart'].iloc[0] - 1000)) & \
                                                                 (compo_csv['send'] <= (narrowed_poss_hits['sstart'].iloc[0] + 50))]
                                narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        if new_align >= 5000:
                            altered = "merge"
                else:
                    new_align = 0
            elif hit_subject_ori == "reverse":
                narrowed_poss_hits = compo_csv[(compo_csv['send'] >= (hit['sstart'] - 50)) & \
                                               (compo_csv['send'] <= (hit['sstart'] + 1000))]
                narrowed_poss_hits = narrowed_poss_hits.sort_values(by=['send'], ascending=True)
                if not narrowed_poss_hits.empty:
                    new_align = hit['qend'] - narrowed_poss_hits['qstart'].iloc[0]
                    if new_align >= 5000:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        altered = "merge"
                    else:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        narrowed_poss_hits_2 = compo_csv[(compo_csv['send'] >= (narrowed_poss_hits['sstart'].iloc[0] - 50)) & \
                                                         (compo_csv['send'] <= (narrowed_poss_hits['sstart'].iloc[0] + 1000))]
                        narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        while not narrowed_poss_hits_2.empty and new_align < 5000:
                            new_align = abs(hit['qend'] - narrowed_poss_hits_2['qstart'].iloc[0])
                            extra_hit = pandas.DataFrame(narrowed_poss_hits_2.iloc[0]).transpose()

                            hit_new = pandas.concat([hit_new, extra_hit.iloc[[0]]], sort=False, ignore_index=False)
                            hit_new = hit_new.reset_index(drop=True)
                            if new_align < 5000:
                                narrowed_poss_hits_2 = compo_csv[(compo_csv['send'] >= (narrowed_poss_hits['sstart'].iloc[0] - 50)) & \
                                                                 (compo_csv['send'] <= (narrowed_poss_hits['sstart'].iloc[0] + 1000))]
                                narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        if new_align >= 5000:
                            altered = "merge"
                else:
                    new_align = 0
        elif mge_ori == "reverse":
            if hit_subject_ori == "forward":

                narrowed_poss_hits = compo_csv[(compo_csv['sstart'] <= (hit['send'] + 1000))  & \
                                               (compo_csv['sstart'] >= (hit['send'] - 50))]
                narrowed_poss_hits = narrowed_poss_hits.sort_values(by=['sstart'], ascending=True)
                if not narrowed_poss_hits.empty:
                    new_align = abs(hit['qstart'] - narrowed_poss_hits['qend'].iloc[0])
                    if new_align >= 5000:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        altered = "merge"
                    else:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        narrowed_poss_hits_2 = compo_csv[
                            (compo_csv['sstart'] >= (narrowed_poss_hits['send'].iloc[0] - 50)) & \
                            (compo_csv['sstart'] <= (narrowed_poss_hits['send'].iloc[0] + 1000))]
                        narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        while not narrowed_poss_hits_2.empty and new_align < 5000:
                            new_align = abs(hit['qstart'] - narrowed_poss_hits_2['qend'].iloc[0])
                            extra_hit = pandas.DataFrame(narrowed_poss_hits_2.iloc[0]).transpose()

                            hit_new = pandas.concat([hit_new, extra_hit.iloc[[0]]], sort=False, ignore_index=False)
                            hit_new = hit_new.reset_index(drop=True)
                            if new_align < 5000:
                                narrowed_poss_hits_2 = compo_csv[
                                    (compo_csv['sstart'] >= (narrowed_poss_hits['send'].iloc[0] - 50)) & \
                                    (compo_csv['sstart'] <= (narrowed_poss_hits['send'].iloc[0] + 1000))]
                                narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        if new_align >= 5000:
                            altered = "merge"
                else:
                    new_align = 0
            elif hit_subject_ori == "reverse":
                narrowed_poss_hits = compo_csv[(compo_csv['sstart'] <= (hit['send'] + 50)) & \
                                               (compo_csv['sstart'] >= (hit['send'] - 1000))]
                narrowed_poss_hits = narrowed_poss_hits.sort_values(by=['sstart'], ascending=False)
                if not narrowed_poss_hits.empty:
                    new_align = abs(hit['qstart'] - narrowed_poss_hits['qend'].iloc[0])
                    if new_align >= 5000:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        altered = "merge"
                    else:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        narrowed_poss_hits_2 = compo_csv[(compo_csv['sstart'] >= (narrowed_poss_hits['send'].iloc[0] - 1000)) & \
                                                         (compo_csv['sstart'] <= (narrowed_poss_hits['send'].iloc[0] + 50))]
                        narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        while not narrowed_poss_hits_2.empty and new_align < 5000:
                            new_align = abs(hit['qstart'] - narrowed_poss_hits_2['qend'].iloc[0])
                            extra_hit = pandas.DataFrame(narrowed_poss_hits_2.iloc[0]).transpose()

                            hit_new = pandas.concat([hit_new, extra_hit.iloc[[0]]], sort=False, ignore_index=False)
                            hit_new = hit_new.reset_index(drop=True)
                            if new_align < 5000:
                                narrowed_poss_hits_2 = compo_csv[(compo_csv['sstart'] >= (narrowed_poss_hits['send'].iloc[0] - 1000)) & \
                                                                 (compo_csv['sstart'] <= (narrowed_poss_hits['send'].iloc[0] + 50))]
                                narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        if new_align >= 5000:
                            altered = "merge"

                else:
                    new_align = 0
    elif ori == "after":
        if mge_ori == "forward":
            if hit_subject_ori == "forward":
                narrowed_poss_hits = compo_csv[(compo_csv['sstart'] <= (hit['send'] + 1000)) &\
                                               (compo_csv['sstart'] >= (hit['send'] - 50))]
                narrowed_poss_hits = narrowed_poss_hits.sort_values(by=['sstart'], ascending=True)
                if not narrowed_poss_hits.empty:
                    new_align = abs(hit['qstart'] - narrowed_poss_hits['qend'].iloc[0])
                    if new_align >= 5000:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        altered = "merge"
                    else:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        narrowed_poss_hits_2 = compo_csv[
                            (compo_csv['sstart'] >= (narrowed_poss_hits['send'].iloc[0] - 50)) & \
                            (compo_csv['sstart'] <= (narrowed_poss_hits['send'].iloc[0] + 1000))]
                        narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        while not narrowed_poss_hits_2.empty and new_align < 5000:
                            new_align = abs(hit['qstart'] - narrowed_poss_hits_2['qend'].iloc[0])
                            extra_hit = pandas.DataFrame(narrowed_poss_hits_2.iloc[0]).transpose()

                            hit_new = pandas.concat([hit_new, extra_hit.iloc[[0]]], sort=False, ignore_index=False)
                            hit_new = hit_new.reset_index(drop=True)
                            if new_align < 5000:
                                narrowed_poss_hits_2 = compo_csv[
                                    (compo_csv['sstart'] >= (narrowed_poss_hits['send'].iloc[0] - 50)) & \
                                    (compo_csv['sstart'] <= (narrowed_poss_hits['send'].iloc[0] + 1000))]
                                narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        if new_align >= 5000:
                            altered = "merge"
                else:
                    new_align = 0
            elif hit_subject_ori == "reverse":
                narrowed_poss_hits = compo_csv[(compo_csv['sstart'] <= (hit['send'] + 50)) & \
                                               (compo_csv['sstart'] >= (hit['send'] - 1000))]
                narrowed_poss_hits = narrowed_poss_hits.sort_values(by=['sstart'], ascending=False)
                if not narrowed_poss_hits.empty:
                    new_align = abs(hit['qstart'] - narrowed_poss_hits['qend'].iloc[0])
                    if new_align >= 5000:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        altered = "merge"
                    else:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        narrowed_poss_hits_2 = compo_csv[
                            (compo_csv['sstart'] >= (narrowed_poss_hits['send'].iloc[0] - 50)) & \
                            (compo_csv['sstart'] <= (narrowed_poss_hits['send'].iloc[0] + 50))]
                        narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        while not narrowed_poss_hits_2.empty and new_align < 5000:
                            new_align = abs(hit['qstart'] - narrowed_poss_hits_2['qend'].iloc[0])
                            extra_hit = pandas.DataFrame(narrowed_poss_hits_2.iloc[0]).transpose()

                            hit_new = pandas.concat([hit_new, extra_hit.iloc[[0]]], sort=False, ignore_index=False)
                            hit_new = hit_new.reset_index(drop=True)
                            if new_align < 5000:
                                narrowed_poss_hits_2 = compo_csv[
                                    (compo_csv['sstart'] >= (narrowed_poss_hits['send'].iloc[0] - 50)) & \
                                    (compo_csv['sstart'] <= (narrowed_poss_hits['send'].iloc[0] + 50))]
                                narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        if new_align >= 5000:
                            altered = "merge"
                        altered = "no"
                else:
                    new_align = 0
        elif mge_ori == "reverse":
            if hit_subject_ori == "forward":
                narrowed_poss_hits = compo_csv[(compo_csv['send'] >= (hit['sstart'] - 1000)) & \
                                               (compo_csv['send'] <= (hit['sstart'] + 50))]
                narrowed_poss_hits = narrowed_poss_hits.sort_values(by = ['send'], ascending = False)

                if not narrowed_poss_hits.empty:
                    new_align = abs(hit['qend'] - narrowed_poss_hits['qstart'].iloc[0])

                    if new_align >= 5000:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        altered = "merge"
                    else:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        narrowed_poss_hits_2 = compo_csv[(compo_csv['send'] >= (narrowed_poss_hits['sstart'].iloc[0] - 1000)) & \
                                                       (compo_csv['send'] <= (narrowed_poss_hits['sstart'].iloc[0] + 50))]
                        narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        while not narrowed_poss_hits_2.empty and new_align < 5000:
                            new_align = abs(hit['qend'] - narrowed_poss_hits_2['qstart'].iloc[0])
                            extra_hit = pandas.DataFrame(narrowed_poss_hits_2.iloc[0]).transpose()

                            hit_new = pandas.concat([hit_new, extra_hit.iloc[[0]]], sort=False, ignore_index=False)
                            hit_new = hit_new.reset_index(drop=True)
                            if new_align < 5000:
                                narrowed_poss_hits_2 = compo_csv[(compo_csv['send'] >= (narrowed_poss_hits['sstart'].iloc[0] - 1000)) & \
                                                                 (compo_csv['send'] <= (narrowed_poss_hits['sstart'].iloc[0] + 50))]
                                narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        if new_align >= 5000:
                            altered = "merge"


                else:
                    new_align = 0
            elif hit_subject_ori == "reverse":
                narrowed_poss_hits = compo_csv[(compo_csv['sstart'] >= (hit['send'] - 50) )& \
                                               (compo_csv['sstart'] <= (hit['send'] + 1000))]
                narrowed_poss_hits = narrowed_poss_hits.sort_values(by=['sstart'], ascending=True)

                if not narrowed_poss_hits.empty:
                    new_align = abs(hit['qend'] - narrowed_poss_hits['qstart'].iloc[0])


                    if new_align >= 5000:
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        altered = "merge"
                        hit_new = pandas.DataFrame(hit_new).transpose()
                        current_hit = pandas.DataFrame(narrowed_poss_hits.iloc[0]).transpose()
                        hit_new = pandas.concat([hit_new, current_hit], sort=False, ignore_index=False)
                        hit_new = hit_new.reset_index(drop=True)
                        narrowed_poss_hits_2 = compo_csv[(compo_csv['sstart'] >= (narrowed_poss_hits['send'].iloc[0] - 50)) & \
                                                         (compo_csv['sstart'] <= (narrowed_poss_hits['send'].iloc[0] + 1000))]
                        narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        while not narrowed_poss_hits_2.empty and new_align < 5000:
                            new_align = abs(hit['qend'] - narrowed_poss_hits_2['qstart'].iloc[0])
                            extra_hit = pandas.DataFrame(narrowed_poss_hits_2.iloc[0]).transpose()

                            hit_new = pandas.concat([hit_new, extra_hit.iloc[[0]]], sort=False, ignore_index=False)
                            hit_new = hit_new.reset_index(drop=True)
                            if new_align < 5000:
                                narrowed_poss_hits_2 = compo_csv[
                                    (compo_csv['sstart'] >= (narrowed_poss_hits['send'].iloc[0] - 50)) & \
                                    (compo_csv['sstart'] <= (narrowed_poss_hits['send'].iloc[0] + 1000))]
                                narrowed_poss_hits_2 = narrowed_poss_hits_2.sort_values(by='align', ascending=False)

                        if new_align >= 5000:
                            altered = "merge"
                else:
                    new_align = 0

    return hit_new, new_align

def multi_hit_merger(hit_new_bef):
    ## Function to merge the many act compos found close to a hit into one hit
    ## on the compo scale.
    ## Input: hit_new_bef: Pandas DF of multiple act comparisons to merge together
    max_row = pandas.to_numeric(hit_new_bef['qend']).idxmax()
    min_row = pandas.to_numeric(hit_new_bef['qstart']).idxmin()
    length_bef = hit_new_bef.iloc[max_row, 7] - hit_new_bef.iloc[min_row, 6]
    hit_bef_out = pandas.Series({'query': hit_new_bef['query'].iloc[0],
                                 'subject': hit_new_bef['subject'].iloc[0],
                                 'pid': hit_new_bef['pid'].iloc[0],
                                 'align': length_bef,
                                 'gap': hit_new_bef['gap'].iloc[0],
                                 'mismatch': hit_new_bef['mismatch'].iloc[0],
                                 'qstart': hit_new_bef.iloc[min_row, 6],
                                 'qend': hit_new_bef.iloc[max_row, 7],
                                 'sstart': hit_new_bef.iloc[min_row, 8],
                                 'send': hit_new_bef.iloc[max_row, 9],
                                 'eval': hit_new_bef['eval'].iloc[0],
                                 'bitscore': hit_new_bef['bitscore'].iloc[0]})

    return hit_bef_out

def hit_mover(hit_before, hit_after, compo_csv, isolate_id, mge_bounds, mge_ori):
    ## Function to include hits after the current if its less than 2k and the new hit
    ## is very close. Finds closest hits, if they are within 50bp of the next nearest.

    hit_before_length = abs(hit_before.iloc[7] - hit_before.iloc[6])
    hit_after_length = abs(hit_after.iloc[7] - hit_after.iloc[6])

    hit_new_bef = hit_before.copy()
    hit_new_aft = hit_after.copy()



    before_pass = False
    after_pass = False

    if hit_before_length < 5000:
        if mge_ori == "forward":
            ## looking for hits previous to the current before
            poss_hit_before = compo_csv[(compo_csv['qend'] <= (hit_before['qstart'] + 50)) &\
                                         (compo_csv['qstart'] >= (mge_bounds[0] - 10)) ]
            poss_hit_before = poss_hit_before.sort_values(by = 'qend', ascending=False)
            new_align = hit_before_length
            if not poss_hit_before.empty:
                current_hit = poss_hit_before[poss_hit_before['qend'] >= (hit_new_bef.iloc[6] - 50)]
                if not current_hit.empty:
                    new_align = hit_new_bef.iloc[7] - current_hit['qstart'].iloc[0]
                    if new_align >= 5000:
                        hit_new_bef = pandas.DataFrame(hit_new_bef).transpose()
                        current_hit = pandas.DataFrame(current_hit.iloc[0]).transpose()
                        hit_new_bef = pandas.concat([hit_new_bef, current_hit], sort=False, ignore_index=False)
                        hit_new_bef = hit_new_bef.reset_index(drop=True)
                    else:

                        hit_new_bef = pandas.DataFrame(hit_new_bef).transpose()
                        current_hit = pandas.DataFrame(current_hit.iloc[0]).transpose()
                        hit_new_bef = pandas.concat([hit_new_bef, current_hit], sort=False, ignore_index=False)
                        hit_new_bef = hit_new_bef.reset_index(drop=True)
                        while_n = 0
                        while not current_hit.empty and new_align < 5000:
                            while_n += 1
                            if while_n != 1:
                                hit_new_bef = pandas.concat([hit_new_bef, current_hit.iloc[[0]]], sort=False,ignore_index=False)
                                hit_new_bef = hit_new_bef.reset_index(drop=True)
                            end_indy = max(hit_new_bef.index.values) - 1
                            current_hit = poss_hit_before[(poss_hit_before['qend'] >= (hit_new_bef.iloc[end_indy, 6] - 50))]


                            if not current_hit.empty:
                                new_align = abs(hit_new_bef.iloc[0, 7] - current_hit.iloc[0, 6])
                                poss_hit_before = poss_hit_before.drop(current_hit.index[0])
                        if new_align < 5000:
                            hit_new_bef = multi_hit_merger(hit_new_bef)
                            hit_new_bef, new_align = subject_checker(hit_new_bef, "before",poss_hit_before, mge_bounds, "forward")


                else:

                    hit_new_bef, new_align = subject_checker(hit_before, "before", poss_hit_before, mge_bounds, "forward")

            if isinstance(hit_new_bef, pandas.DataFrame):
                before_process = "merge"
            else:
                before_process = "no"

        elif mge_ori == "reverse":
            poss_hit_before = compo_csv[(compo_csv['qstart'] >= (hit_before['qend'] - 50)) & \
                                        (compo_csv['qend'] <= (mge_bounds[1] + 10)) ]
            poss_hit_before = poss_hit_before.sort_values(by='qstart', ascending=True)
            new_align = hit_before_length
            if not poss_hit_before.empty:
                current_hit = poss_hit_before[poss_hit_before['qstart'] <= (hit_new_bef.iloc[7] + 50)]
                if not current_hit.empty:
                    new_align = abs(hit_new_bef.iloc[6] - current_hit['qend'].iloc[0])

                    if new_align >= 5000:
                        hit_new_bef = pandas.DataFrame(hit_new_bef).transpose()
                        current_hit = pandas.DataFrame(current_hit.iloc[0]).transpose()
                        hit_new_bef = pandas.concat([hit_new_bef, current_hit], sort=False, ignore_index=False)
                        hit_new_bef = hit_new_bef.reset_index(drop=True)

                    else:
                        hit_new_bef = pandas.DataFrame(hit_new_bef).transpose()
                        current_hit = pandas.DataFrame(current_hit.iloc[0]).transpose()
                        hit_new_bef = pandas.concat([hit_new_bef, current_hit], sort=False, ignore_index=False)
                        hit_new_bef = hit_new_bef.reset_index(drop=True)
                        while_n = 0
                        while not current_hit.empty and new_align < 5000:
                            while_n += 1
                            if while_n != 1:
                                hit_new_bef = pandas.concat([hit_new_bef, current_hit.iloc[[0]]], sort=False, ignore_index=False)
                                hit_new_bef = hit_new_bef.reset_index(drop=True)
                            end_indy = max(hit_new_bef.index.values) - 1
                            current_hit = poss_hit_before[poss_hit_before['qstart'] <= (hit_new_bef.iloc[end_indy, 7] + 50)]
                            if not current_hit.empty:
                                new_align = abs(hit_new_bef.iloc[0, 6] - current_hit.iloc[0, 7])
                                poss_hit_before = poss_hit_before.drop(current_hit.index[0])
                        if new_align < 5000:
                            hit_new_bef = multi_hit_merger(hit_new_bef)
                            hit_new_bef, new_align = subject_checker(hit_new_bef, "before",poss_hit_before, mge_bounds, "reverse")
                else:

                    hit_new_bef, new_align = subject_checker(hit_before, "before", poss_hit_before, mge_bounds, "reverse")

            if isinstance(hit_new_bef, pandas.DataFrame):
                before_process = "merge"
            else:
                before_process = "no"


    else:
        before_pass = True
        before_process = "No"

    if hit_after_length < 5000:
        if mge_ori == "forward":
            poss_hit_after = compo_csv[(compo_csv['qstart'] >= (hit_after['qend'] - 50)) & \
                                        (compo_csv['qend'] <= (mge_bounds[1] + 10)) ]
            poss_hit_after = poss_hit_after.sort_values(by='qstart', ascending=True)
            new_align = hit_after_length
            if not poss_hit_after.empty:
                current_hit = poss_hit_after[poss_hit_after['qstart'] <= (hit_new_aft.iloc[7] + 50)]
                if not current_hit.empty:
                    new_align = abs(hit_new_aft.iloc[6] - current_hit['qend'].iloc[0])
                    if new_align >= 5000:
                        hit_new_aft = pandas.DataFrame(hit_new_aft).transpose()
                        current_hit = pandas.DataFrame(current_hit.iloc[0]).transpose()
                        hit_new_aft = pandas.concat([hit_new_aft, current_hit], sort=False, ignore_index=False)
                        hit_new_aft = hit_new_aft.reset_index(drop=True)

                    else:
                        hit_new_aft = pandas.DataFrame(hit_new_aft).transpose()
                        current_hit = pandas.DataFrame(current_hit.iloc[0]).transpose()
                        hit_new_aft = pandas.concat([hit_new_aft, current_hit], sort=False, ignore_index=False)
                        hit_new_aft = hit_new_aft.reset_index(drop=True)
                        while_n = 0
                        while not current_hit.empty and new_align < 5000:
                            while_n += 1
                            if while_n != 1:
                                hit_new_aft = pandas.concat([hit_new_aft, current_hit.iloc[[0]]], sort=False, ignore_index=False)
                                hit_new_aft = hit_new_aft.reset_index(drop=True)
                            end_indy = max(hit_new_aft.index.values) - 1
                            current_hit = poss_hit_after[poss_hit_after['qstart'] <= (hit_new_aft.iloc[end_indy, 7] + 50)]
                            if not current_hit.empty:
                                new_align = abs(hit_new_aft.iloc[0, 6] - current_hit.iloc[0, 7])
                                poss_hit_after = poss_hit_after.drop(current_hit.index[0])
                        if new_align < 5000:
                            hit_new_aft = multi_hit_merger(hit_new_aft)
                            hit_new_aft, new_align = subject_checker(hit_new_aft, "after",poss_hit_after, mge_bounds, "forward")
            else:
                hit_new_aft, new_align = subject_checker(hit_after, "after", poss_hit_after, mge_bounds, "forward")

            if isinstance(hit_new_aft, pandas.DataFrame):
                after_process = "merge"
            else:
                after_process = "no"
        elif mge_ori == "reverse":
            ## looking for hits previous to the current before
            poss_hit_after = compo_csv[(compo_csv['qend'] <= (hit_after['qstart'] + 50)) & \
                                        (compo_csv['qstart'] >= (mge_bounds[0] - 10)) ]
            poss_hit_after = poss_hit_after.sort_values(by='qend', ascending=False)
            new_align = hit_after_length
            if not poss_hit_after.empty:
                current_hit = poss_hit_after[poss_hit_after['qend'] >= (hit_new_aft.iloc[ 6] - 50)]
                if not current_hit.empty:

                    new_align = hit_new_aft.iloc[7] - current_hit['qstart'].iloc[0]
                    if new_align >= 5000:
                        hit_new_aft = pandas.DataFrame(hit_new_aft).transpose()
                        current_hit = pandas.DataFrame(current_hit.iloc[0]).transpose()
                        hit_new_aft = pandas.concat([hit_new_aft, current_hit], sort=False, ignore_index=False)
                        hit_new_aft = hit_new_aft.reset_index(drop=True)

                    else:
                        hit_new_aft = pandas.DataFrame(hit_new_aft).transpose()
                        current_hit = pandas.DataFrame(current_hit.iloc[0]).transpose()
                        hit_new_aft = pandas.concat([hit_new_aft, current_hit], sort=False, ignore_index=False)
                        hit_new_aft = hit_new_aft.reset_index(drop=True)
                        while_n = 0
                        while not current_hit.empty and new_align < 5000:
                            while_n += 1
                            if while_n != 1:
                                hit_new_aft = pandas.concat([hit_new_aft, current_hit.iloc[[0]]], sort=False, ignore_index=False)
                                hit_new_aft = hit_new_aft.reset_index(drop=True)
                            end_indy = max(hit_new_aft.index.values) - 1
                            current_hit = poss_hit_after[poss_hit_after['qend'] >= (hit_new_aft.iloc[end_indy, 6] - 50)]
                            if not current_hit.empty:
                                new_align = hit_new_aft.iloc[0, 7] - current_hit.iloc[0, 6]
                                poss_hit_after = poss_hit_after.drop(current_hit.index[0])
                        if new_align < 5000:
                            hit_new_aft = multi_hit_merger(hit_new_aft)
                            hit_new_aft, new_align = subject_checker(hit_new_aft, "after",poss_hit_after, mge_bounds, "forward")


                else:
                    hit_new_aft, new_align = subject_checker(hit_after, "after", poss_hit_after, mge_bounds, "reverse")

            if isinstance(hit_new_aft, pandas.DataFrame):

                after_process = "merge"
            else:
                after_process = "no"
    else:
        after_pass = True
        after_process = "No"


    if before_process == "merge":
        max_row = pandas.to_numeric(hit_new_bef['qend']).idxmax()
        min_row = pandas.to_numeric(hit_new_bef['qstart']).idxmin()
        length_bef = hit_new_bef.iloc[max_row, 7] - hit_new_bef.iloc[min_row, 6]
        hit_bef_out = pandas.Series({'query': hit_before['query'],
                                 'subject': hit_before['subject'],
                                 'pid': hit_before['pid'],
                                 'align': length_bef,
                                 'gap': hit_before['gap'],
                                 'mismatch': hit_before['mismatch'],
                                 'qstart': hit_new_bef.iloc[min_row, 6],
                                 'qend': hit_new_bef.iloc[max_row, 7],
                                 'sstart': hit_new_bef.iloc[min_row,8],
                                 'send': hit_new_bef.iloc[max_row, 9],
                                 'eval': hit_before['eval'],
                                 'bitscore': hit_before['bitscore']})

        before_pass = length_bef >= 5000
    else:
        hit_bef_out = hit_before.copy()

    if after_process == "merge":
        max_row = pandas.to_numeric(hit_new_aft['qend']).idxmax()
        min_row = pandas.to_numeric(hit_new_aft['qstart']).idxmin()
        length_aft = hit_new_aft.iloc[max_row, 7] - hit_new_aft.iloc[min_row, 6]
        hit_aft_out = pandas.Series({'query': hit_after['query'],
                                     'subject': hit_after['subject'],
                                     'pid': hit_after['pid'],
                                     'align': length_aft,
                                     'gap': hit_after['gap'],
                                     'mismatch': hit_after['mismatch'],
                                     'qstart': hit_new_aft.iloc[min_row, 6],
                                     'qend': hit_new_aft.iloc[max_row, 7],
                                     'sstart': hit_new_aft.iloc[min_row, 8],
                                     'send': hit_new_aft.iloc[max_row, 9],
                                     'eval': hit_after['eval'],
                                     'bitscore': hit_after['bitscore']})
        after_pass = length_aft >= 5000
    else:
        hit_aft_out = hit_after.copy()


    tot_pass = before_pass and after_pass


    return hit_bef_out, hit_aft_out, tot_pass

def final_acceptor(hit_before, hit_after, isolate_id, mge_bounds, mge_ori):
    ## Function to check whether the merged bounds are above 2500 and then
    ## to just take the 5000 from here if it's within the contig bounds
    ## input: hit_before: Merged hit before from hit_mover
    ##        hit_after: Merged hit after from hit_mover
    ##        compo_csv: The act comparison csv
    ##        isolate_id: The current isolate_id
    ##        mge_bounds: The contig bounds for the mge
    ##        mge_ori: The orientation of the mge

    hit_before_length = abs(hit_before.iloc[7] - hit_before.iloc[6])
    hit_after_length = abs(hit_after.iloc[7] - hit_after.iloc[6])

    hit_new_bef = hit_before.copy()
    hit_new_aft = hit_after.copy()

    if hit_before_length < 5000:
        if mge_ori == "forward":
            new_end = hit_before.iloc[7] - 5000
            if new_end > mge_bounds[0]:
                max_send = max([hit_before['sstart'], hit_before['send']])

                hit_bef_out = pandas.Series({'query': hit_before['query'],
                                             'subject': hit_before['subject'],
                                             'pid': hit_before['pid'],
                                             'align': 5000,
                                             'gap': hit_before['gap'],
                                             'mismatch': hit_before['mismatch'],
                                             'qstart': new_end,
                                             'qend': hit_new_bef['qend'],
                                             'sstart': max_send - 5000,
                                             'send': max_send,
                                             'eval': hit_before['eval'],
                                             'bitscore': hit_before['bitscore']})
                hit_before_length = 5000
            else:
                hit_bef_out = hit_new_bef
        elif mge_ori == "reverse":
            new_end = hit_before.iloc[6] + 5000
            if new_end > mge_bounds[0]:
                min_send = min([hit_before['sstart'], hit_before['send']])
                hit_bef_out = pandas.Series({'query': hit_before['query'],
                                             'subject': hit_before['subject'],
                                             'pid': hit_before['pid'],
                                             'align': 5000,
                                             'gap': hit_before['gap'],
                                             'mismatch': hit_before['mismatch'],
                                             'qstart': hit_before['qstart'],
                                             'qend': new_end,
                                             'sstart': min_send,
                                             'send': min_send + 5000,
                                             'eval': hit_before['eval'],
                                             'bitscore': hit_before['bitscore']})
                hit_before_length = 5000
            else:
                hit_bef_out = hit_new_bef
    else:
        hit_bef_out = hit_new_bef

    ## Now for the after hits
    if hit_after_length < 5000:
        if mge_ori == "reverse":
            new_end = hit_after.iloc[6] + 5000
            if new_end > mge_bounds[0]:
                min_send = min([hit_after['sstart'], hit_after['send']])
                hit_aft_out = pandas.Series({'query': hit_after['query'],
                                             'subject': hit_after['subject'],
                                             'pid': hit_after['pid'],
                                             'align': 5000,
                                             'gap': hit_after['gap'],
                                             'mismatch': hit_after['mismatch'],
                                             'qstart': hit_after['qstart'],
                                             'qend': new_end,
                                             'sstart': min_send,
                                             'send': min_send + 5000,
                                             'eval': hit_after['eval'],
                                             'bitscore': hit_after['bitscore']})
                hit_after_length = 5000
            else:
                hit_aft_out = hit_new_aft
        elif mge_ori == "forward":
            new_end = hit_after.iloc[7] - 5000
            if new_end > mge_bounds[0]:
                max_send = max([hit_after['sstart'], hit_after['send']])
                hit_aft_out = pandas.Series({'query': hit_after['query'],
                                             'subject': hit_after['subject'],
                                             'pid': hit_after['pid'],
                                             'align': 5000,
                                             'gap': hit_after['gap'],
                                             'mismatch': hit_after['mismatch'],
                                             'qstart': new_end,
                                             'qend': hit_new_aft['qend'],
                                             'sstart': max_send - 5000,
                                             'send': max_send,
                                             'eval': hit_after['eval'],
                                             'bitscore': hit_after['bitscore']})
                hit_after_length = 5000
            else:
                hit_aft_out = hit_new_aft
    else:
        hit_aft_out = hit_new_aft

    if hit_before_length >= 5000 and hit_after_length >= 5000:
        both_tigs = True
    else:
        both_tigs = False

    return  hit_bef_out, hit_aft_out, both_tigs

def before_and_after_check(hit_to_check, mge_locs, compo_table, hit_loc, other_hit, isolate_id):
    ## This function will take input from the hit distance checker. We're looking for
    ## any hits that might have been missed by our initial search in before and after
    ## hits function, i.e they have an align < 2000. They have to be adjacent to the
    ## mge and adjacent to the hit_to_check in the reference genome.

    hit_check_query = hit_to_check[[6, 7]]
    hit_check_subject = hit_to_check[[8, 9]]
    other_hit_subject = other_hit[[8,9]]

    if mge_locs[0] < mge_locs[1]:
        mge_ori = "forward"
    elif mge_locs[0] > mge_locs[1]:
        mge_ori = "reverse"

    if hit_check_subject[0] < hit_check_subject[1]:
        check_ori = "forward"
    elif hit_check_subject[0] > hit_check_subject[1]:
        check_ori = "reverse"

    added_hits = pandas.DataFrame()

    if hit_loc == "before":
        if mge_ori == "forward":
            if check_ori == "forward":
                poss_hits = compo_table[(compo_table['qend'] < (mge_locs[0] + 50)) &\
                                        (compo_table['qend'] > hit_check_query[0]) &\
                                        (compo_table['send'] > hit_check_subject[1]) &\
                                        (compo_table['sstart'] < other_hit_subject[1])]
                poss_hits = poss_hits.sort_values(by=['qstart'], ascending=True)
                if not poss_hits.empty:
                    hit_orig = pandas.DataFrame(hit_to_check).transpose()
                    poss_hits = poss_hits.reset_index(drop=True)
                    current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                    added_hits = pandas.concat([added_hits,hit_orig, current_hit], ignore_index=True, sort=False)
                    poss_hits = poss_hits.drop(0)
                    new_row = 0
                    while not poss_hits.empty:
                        new_row += 1
                        current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                        added_hits = pandas.concat([added_hits, current_hit], ignore_index=True, sort=False)
                        poss_hits = poss_hits.drop(new_row)



            elif check_ori == "reverse":
                poss_hits = compo_table[(compo_table['qend'] < (mge_locs[0] + 50)) & \
                                        (compo_table['qend'] > hit_check_query[0]) & \
                                        (compo_table['send'] < hit_check_subject[1]) &\
                                        (compo_table['sstart'] > other_hit_subject[1])]
                poss_hits = poss_hits.sort_values(by=['qend'], ascending=True)
                if not poss_hits.empty:
                    hit_orig = pandas.DataFrame(hit_to_check).transpose()
                    poss_hits = poss_hits.reset_index(drop=True)
                    current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                    added_hits = pandas.concat([added_hits,hit_orig, current_hit], ignore_index=True, sort=False)
                    poss_hits = poss_hits.drop(0)
                    new_row = 0
                    while not poss_hits.empty:
                        new_row += 1
                        current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                        added_hits = pandas.concat([added_hits, current_hit], ignore_index=True, sort=False)
                        poss_hits = poss_hits.drop(new_row)

        elif mge_ori == "reverse":
            if check_ori == "forward":
                poss_hits = compo_table[(compo_table['qstart'] > (mge_locs[1] - 50)) & \
                                        (compo_table['qstart'] < hit_check_query[1]) & \
                                        (compo_table['sstart'] < hit_check_subject[0]) & \
                                        (compo_table['send'] > other_hit_subject[0])]
                poss_hits = poss_hits.sort_values(by=['qstart'], ascending=False)
                if not poss_hits.empty:
                    hit_orig = pandas.DataFrame(hit_to_check).transpose()
                    poss_hits = poss_hits.reset_index(drop=True)
                    current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                    added_hits = pandas.concat([added_hits, hit_orig, current_hit], ignore_index=True, sort=False)
                    poss_hits = poss_hits.drop(0)
                    new_row = 0
                    while not poss_hits.empty:
                        new_row += 1
                        current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                        added_hits = pandas.concat([added_hits, current_hit], ignore_index=True, sort=False)
                        poss_hits = poss_hits.drop(new_row)

            elif check_ori == "reverse":
                poss_hits = compo_table[(compo_table['qstart'] > (mge_locs[1] - 50)) & \
                                        (compo_table['qstart'] < hit_check_query[1]) & \
                                        (compo_table['sstart'] > hit_check_subject[0]) &\
                                        (compo_table['send'] < other_hit_subject[0])]
                poss_hits = poss_hits.sort_values(by=['qstart'], ascending=False)
                if not poss_hits.empty:
                    hit_orig = pandas.DataFrame(hit_to_check).transpose()
                    poss_hits = poss_hits.reset_index(drop=True)
                    current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                    added_hits = pandas.concat([added_hits, hit_orig, current_hit], ignore_index=True, sort=False)
                    poss_hits = poss_hits.drop(0)
                    new_row = 0
                    while not poss_hits.empty:
                        new_row += 1
                        current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                        added_hits = pandas.concat([added_hits, current_hit], ignore_index=True, sort=False)
                        poss_hits = poss_hits.drop(new_row)

    elif hit_loc == "after":
        if mge_ori == "forward":
            if check_ori == "forward":
                poss_hits = compo_table[(compo_table['qstart'] > (mge_locs[1] - 50)) & \
                                        (compo_table['qstart'] < hit_check_query[1]) & \
                                        (compo_table['sstart'] < hit_check_subject[0])&\
                                        (compo_table['send'] > other_hit_subject[0])]
                poss_hits = poss_hits.sort_values(by=['qstart'], ascending=False)
                if not poss_hits.empty:
                    hit_orig = pandas.DataFrame(hit_to_check).transpose()
                    poss_hits = poss_hits.reset_index(drop=True)
                    current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                    added_hits = pandas.concat([added_hits, hit_orig, current_hit], ignore_index=True, sort=False)
                    poss_hits = poss_hits.drop(0)
                    new_row = 0
                    while not poss_hits.empty:
                        new_row += 1
                        current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                        added_hits = pandas.concat([added_hits, current_hit], ignore_index=True, sort=False)
                        poss_hits = poss_hits.drop(new_row)


            elif check_ori == "reverse":
                poss_hits = compo_table[(compo_table['qstart'] > (mge_locs[1] - 50)) & \
                                        (compo_table['qstart'] < hit_check_query[1]) & \
                                        (compo_table['sstart'] > hit_check_subject[0]) &\
                                        (compo_table['send'] < other_hit_subject[0])]
                poss_hits = poss_hits.sort_values(by=['qstart'], ascending=False)
                if not poss_hits.empty:
                    hit_orig = pandas.DataFrame(hit_to_check).transpose()
                    poss_hits = poss_hits.reset_index(drop=True)
                    current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                    added_hits = pandas.concat([added_hits, hit_orig, current_hit], ignore_index=True, sort=False)
                    poss_hits = poss_hits.drop(0)
                    new_row = 0
                    while not poss_hits.empty:
                        new_row += 1
                        current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                        added_hits = pandas.concat([added_hits, current_hit], ignore_index=True, sort=False)
                        poss_hits = poss_hits.drop(new_row)

        elif mge_ori == "reverse":
            if check_ori == "forward":
                poss_hits = compo_table[(compo_table['qend'] < (mge_locs[0] + 50)) & \
                                        (compo_table['qend'] > hit_check_query[0]) & \
                                        (compo_table['send'] > hit_check_subject[1]) &\
                                        (compo_table['sstart'] < other_hit_subject[1])]
                poss_hits = poss_hits.sort_values(by=['qend'], ascending=True)
                if not poss_hits.empty:
                    hit_orig = pandas.DataFrame(hit_to_check).transpose()
                    poss_hits = poss_hits.reset_index(drop=True)
                    current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                    added_hits = pandas.concat([added_hits, hit_orig, current_hit], ignore_index=True, sort=False)
                    poss_hits = poss_hits.drop(0)
                    new_row = 0
                    while not poss_hits.empty:
                        new_row += 1
                        current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                        added_hits = pandas.concat([added_hits, current_hit], ignore_index=True, sort=False)
                        poss_hits = poss_hits.drop(new_row)


            elif check_ori == "reverse":
                poss_hits = compo_table[(compo_table['qend'] < (mge_locs[0] + 50)) & \
                                        (compo_table['qend'] > hit_check_query[0]) & \
                                        (compo_table['send'] < hit_check_subject[1]) &\
                                        (compo_table['sstart'] > other_hit_subject[1])]
                poss_hits = poss_hits.sort_values(by=['qend'], ascending=True)
                if not poss_hits.empty:
                    hit_orig = pandas.DataFrame(hit_to_check).transpose()
                    poss_hits = poss_hits.reset_index(drop=True)
                    current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                    added_hits = pandas.concat([added_hits, hit_orig, current_hit], ignore_index=True, sort=False)
                    poss_hits = poss_hits.drop(0)
                    new_row = 0
                    while not poss_hits.empty:
                        new_row += 1
                        current_hit = pandas.DataFrame(poss_hits.iloc[0]).transpose()
                        added_hits = pandas.concat([added_hits, current_hit], ignore_index=True, sort=False)
                        poss_hits = poss_hits.drop(new_row)


    return added_hits

def isolate_name_producer(names_series):
    ## Function to turn a series of isolate gff/fasta locations into a list of isolates:
    ## Input: names_series: A series object of fasta locations to work through

    ids = []

    for k in range(len(names_series.index)):
        current_id = names_series.iloc[k]
        current_base = os.path.basename(current_id)
        current_isolate = re.sub("\..*[a-z,A-Z].*$","", current_base)
        ids.append(current_isolate)

    return  ids

def all_presence_checker(id_list, series_ids):
    ## Function to check whether the total list of ids is within the blast csv
    ## Input: id_list: List of ids to check for presence all must be within the series ids
    ##        series_ids: the series of ids to check through
    skip = False
    num_in = 0
    for k in range(len(id_list)):
        if any(series_ids.isin([id_list[k]])):
            num_in += 1

    if num_in == len(id_list):
        skip = True

    return skip

def nice_proper_hits_ids(ids_list):
    ## Function to turn all the double underscore ids into standard isolate names for the act compos:
    ## Input: ids_list: list format of the first column of the proper hits
    nicer = []
    for id in ids_list:
        if id.count('_') == 2:
            last_occy = id.rfind('_')
            isolate_id = id[0:last_occy]
            nicer.append(isolate_id)
        else:
            nicer.append(id)

    return nicer

def act_mapper(hit_before, hit_after, act_loc, current_insert_locs):
    ## Function to take in the before and after hits generated from a new reference and convert
    ## the before and after hits sstart and send back into hits within the gubbins reference.
    ## Input: hit_before: The calculated hit before in the new reference
    ##        hit_after: The calculate hit after in the new reference
    ##        act_loc: The location of the act_comparison between the old reference and the new reference
    ##        current_insert_locs: The locs in the new reference with the insert in between.
    ## So I think we'll just calculate the location of the two closest points in both (so representative of either
    ## side of the insertion) then take the hit the other two as +- 5000 bp.

    ## Load up the act compo
    compo_names = ['query', 'subject', 'pid', 'align', 'gap', 'mismatch', 'qstart',
                   'qend', 'sstart', 'send', 'eval', 'bitscore']
    compo_table = pandas.read_csv(act_loc, sep='\t', names=compo_names)

    nuevo_before = hit_before.copy()
    nuevo_after = hit_after.copy()
    ## first test if both are within a single hit likely from an insertion away form where mge is within reference

    single_hit = compo_table[(compo_table['qstart'] <= current_insert_locs[0]) & (compo_table['qend'] >= current_insert_locs[1])]

    if single_hit.empty:
        ## Maybe this represent the insertion of the mge in the element too so lets look for hits either side
        hit_1 = compo_table[(compo_table['qend'] >= (current_insert_locs[0] - 50)) & (compo_table['qstart'] <= (current_insert_locs[1] + 50))]
        hit_2 = compo_table[(compo_table['qstart'] >= (current_insert_locs[0] - 50)) & (compo_table['qend'] <= (current_insert_locs[1] + 50))]
        hit_1 = hit_1.sort_values(by='align', ascending=False)
        hit_2 = hit_2.sort_values(by='align', ascending=False)

        if hit_1.empty or hit_2.empty:
            return  hit_before, hit_after, False
        else:
            if isinstance(hit_1, pandas.Series):
                hit_1 = hit_1.to_frame.transpose()
            if isinstance(hit_2, pandas.Series):
                hit_2 = hit_2.to_frame.transpose()
            send_gap = hit_1['qend'].iloc[0] - current_insert_locs[0]
            sstart_gap = current_insert_locs[1] - hit_2['qstart'].iloc[0]

            if hit_1['sstart'].iloc[0] < hit_1['send'].iloc[0]:
                new_sstart = hit_1['send'].iloc[0] - send_gap
                bef_adder = -5000
            else:
                new_sstart = hit_1['send'].iloc[0] + send_gap
                bef_adder = 5000
            if hit_2['sstart'].iloc[0] < hit_2['send'].iloc[0]:
                new_send = hit_2['sstart'].iloc[0] + sstart_gap
                aft_adder = 5000
            else:
                new_send = hit_2['sstart'].iloc[0] - sstart_gap
                aft_adder = -5000

            nuevo_before['send'] = new_sstart
            nuevo_before['sstart'] = new_sstart + bef_adder
            nuevo_after['sstart'] = new_send
            nuevo_after['send'] = new_send + aft_adder

            mapped = True
    else:
        single_hit = single_hit.sort_values(by='align', ascending=False)
        if isinstance(single_hit, pandas.Series):
            single_hit = single_hit.to_frame().transpose()

        send_gap = single_hit['qend'].iloc[0] - current_insert_locs[1]
        sstart_gap = current_insert_locs[0] - single_hit['qstart'].iloc[0]

        if single_hit['sstart'].iloc[0] < single_hit['send'].iloc[0]:
            new_sstart = single_hit['sstart'].iloc[0] + sstart_gap
            new_send = single_hit['send'].iloc[0] - send_gap
            bef_adder = -5000
            aft_adder = 5000
        else:
            new_sstart = single_hit['sstart'].iloc[0] - sstart_gap
            new_send = single_hit['send'].iloc[0] + send_gap
            bef_adder = 5000
            aft_adder = -5000

        nuevo_before['send'] = new_sstart
        nuevo_before['sstart'] = new_sstart + bef_adder
        nuevo_after['sstart'] = new_send
        nuevo_after['send'] = new_send + aft_adder
        mapped = True

    return  nuevo_before, nuevo_after, mapped


if __name__ == '__main__':
    tab = str.maketrans("ACTG", "TGAC")
    tic = time.perf_counter()
    files_for_input = get_options()

    python_dir = os.path.dirname(os.path.abspath(__file__))
    perl_dir = re.sub("python$", "perl", python_dir)

    pandas.set_option('display.max_columns', 500)

    contig_file_abs_path = files_for_input.contig_loc
    absolute_act_path = files_for_input.act_loc
    fasta_csv = files_for_input.fasta_csv
    fasta_pandas = pandas.read_csv(fasta_csv)
    fasta_pandas.columns = ['isolate', 'reference']




    proper_hits = pandas.read_csv(files_for_input.blast_csv)
    merged_csv, merged_locs = merged_contig_checker(proper_hits, contig_file_abs_path, absolute_act_path)
    is_2k = merged_csv['align'] >= int(files_for_input.align)

    proper_hits = merged_csv[is_2k]
    proper_hits = proper_hits.reset_index(drop=True)

    nice_ids = nice_proper_hits_ids(proper_hits.iloc[:, 0].tolist())
    proper_hits['nicest_ids'] = pandas.Series(nice_ids, index=proper_hits.index)
    ## Now lets load up the csv with the isolate names and their reference location

    isolate_ref_gff = pandas.read_csv(files_for_input.reference_csv)

    library_dat = pandas.read_csv(files_for_input.lib_csv)

    ## Now lets remove the library ids from the proper hits index

    narrowed_prop = proper_hits[proper_hits['nicest_ids'].isin(library_dat['id']) == False]

    ## Ok so narrowed prop now contains all the hits we have to loop through and assign as hits in the
    ## library or not


    library_dat['insert_name'] = pandas.Series(range(1,(len(library_dat.index) + 1)))

    refs_to_alter = []
    clusters_to_skip = []
    new_refs = []
    ## Set up the csv for isolates that were missed due to incomplete hits
    missing_colnames = ["id","mge_start","mge_end","mge_length","ref_name","cluster_name","reason"]
    missing_isolate = pandas.DataFrame(columns=missing_colnames)

    ## Set up the csv for the isolates with matching hits

    hit_col_names = ["id", "mge_start", "mge_end", "insert_start", "insert_end","before_sstart","before_send","after_sstart",
                     "after_send", "mge_length",
                     "insert_length", "insert_genes", "mge_genes", "flank_genes",
                     "mean_flank_gene_length", 'flanks_length', 'before_flank_gene', 'after_flank_gene',
                     'before_flank_avg',
                     'after_flank_avg', 'before_gene_name', 'after_gene_name', 'ref_name', 'cluster_name',
                     'insert_name']

    hit_df = pandas.DataFrame(columns=hit_col_names)
    tic_1 = time.perf_counter()
    print("")
    print("This many hits to get through: %s" % len(proper_hits.index))
    print("")

    ## Now loop through the blast results ##

    for k in range(len(narrowed_prop.index)):
        ## First we'll get the hit locs in place for each of the hits
        start_len = len(missing_isolate.index) + len(hit_df.index)
        print("On isolate: ", k, end='\r', flush=True)

        current_row = narrowed_prop.iloc[[k]]
        hitters = (list(current_row.iloc[0, [5, 6]]))
        if hitters[0] < hitters[1]:
            mge_ori = "forward"
        elif hitters[0] > hitters[1]:
            mge_ori = "reverse"

        ## Now we'll go into the tab files and get the compo files
        isolate_id_z = narrowed_prop.iloc[k, 0]
        if isolate_id_z.count('_') == 2:
            last_occy = isolate_id_z.rfind('_')
            isolate_id = isolate_id_z[0:last_occy]
        else:
            isolate_id = isolate_id_z

        current_mge_length = abs(hitters[1] - hitters[0])
        # cluster_2 = ["10900_6#9", "15608_3#45", "15608_3#49", "15608_3#55", "15608_3#71", "15608_3#75", "15608_4#18", \
        #            "15608_4#30", "15608_4#2", "15608_4#24"]
        # cluster_2 = ["6678_3#5","6569_4#15"]

        # cluster_2 = ["12291_5#65","15682_1#29","15682_1#38","15682_1#51","15682_1#65","15682_2#52","20925_3#47","21127_1#10",\
        #             "21127_1#104","21127_1#108","21127_1#35","21127_1#44","21127_1#67","21127_1#70","21127_1#83","21127_1#97",\
        #             "22841_4#117","6569_4#16","6569_4#17","6569_4#18"]
        # cluster_2 = ["15608_5#88","15682_1#66", "19084_7#62", "19084_7#68","20925_3#60","20925_3#68","6187_5#9"]

        # if isolate_id not in cluster_2:
        #    continue

        current_gff_loc, ref_loc, cluster_name = gff_finder(isolate_ref_gff, isolate_id, True)
        ref_name = os.path.basename(ref_loc.iloc[0])
        ref_name = re.sub("\..*[a-zA-Z]*$", "", ref_name)
        cluster_name = cluster_name.iloc[0]

        act_map = False
        if cluster_name in clusters_to_skip:
            missing_current = pandas.DataFrame()
            missing_current['id'] = pandas.Series(isolate_id)
            missing_current['mge_start'] = pandas.Series(hitters[0], index=missing_current.index)
            missing_current['mge_end'] = pandas.Series(hitters[1], index=missing_current.index)
            missing_current['ref_name'] = pandas.Series(ref_name, index=missing_current.index)
            missing_current['cluster_name'] = pandas.Series(cluster_name, index=missing_current.index)
            missing_current['mge_length'] = pandas.Series(current_mge_length, index=missing_current.index)
            missing_current['reason'] = pandas.Series(["In clusters to skip"], index=missing_current.index)
            missing_isolate = missing_isolate.append(missing_current, sort = False)

            continue

        if ref_name not in refs_to_alter:

            if ref_name in nice_ids:
                ## Need to re-run act for this cluster with a different id. First check if all the cluster
                ## are in the isolate_ids. If so, we have to skip the whole cluster
                current_cluster_vals = isolate_ref_gff[isolate_ref_gff['reference'] == ref_loc.iloc[0]]

                current_ids = isolate_name_producer(current_cluster_vals['isolate'])
                skip = all_presence_checker(current_ids, narrowed_prop.iloc[:, 0])
                if skip:
                    missing_current = pandas.DataFrame()
                    missing_current['id'] = pandas.Series(isolate_id)
                    missing_current['mge_start'] = pandas.Series(hitters[0], index=missing_current.index)
                    missing_current['mge_end'] = pandas.Series(hitters[1], index=missing_current.index)
                    missing_current['ref_name'] = pandas.Series(ref_name,  index=missing_current.index)
                    missing_current['cluster_name'] = pandas.Series(cluster_name,  index=missing_current.index)
                    missing_current['mge_length'] = pandas.Series(current_mge_length,  index=missing_current.index)
                    missing_current['reason'] = pandas.Series(["In clusters to skip"],  index=missing_current.index)
                    missing_isolate = missing_isolate.append(missing_current)
                    clusters_to_skip.append(cluster_name, sort = False)
                    continue
                else:
                    ## Find n50 of those not in the blast file  re-run act compos
                    all_ids = nice_ids
                    no_element_ids = numpy.setdiff1d(current_ids, all_ids).tolist()

                    element_ids = list(set(current_ids) & set(all_ids))
                    ## no element ids contains the isolates without the element in this cluster
                    ## Now we need to find their fastas and then run through the act comparison
                    fastas_to_run = []
                    for fasta_fn in no_element_ids:
                        current_loc, ref_loc, ignore_me = gff_finder(fasta_pandas, fasta_fn, False)
                        fastas_to_run.append(current_loc.iloc[0])

                    fastas_to_act = []
                    for fasta_fn in element_ids:
                        current_loc, ref_loc, ignore_me = gff_finder(fasta_pandas, fasta_fn, False)
                        fastas_to_act.append(current_loc.iloc[0])

                    ## Now we've got the fasta list lets get the best n50 loc

                    new_ref = n50_calc(fastas_to_run, ref_name)
                    new_ref_name = os.path.basename(re.sub("\..*[a-z,A-Z].*$", "", new_ref))
                    ## Now create the new fastas df to input to the act compo script
                    new_act_df = pandas.DataFrame()
                    new_act_df['isolate'] = pandas.Series(fastas_to_act)
                    new_act_df['reference'] = pandas.Series(numpy.repeat(new_ref, len(fastas_to_act)).tolist(),index=new_act_df.index)

                    df_loc = "./" + cluster_name + "_altered_fasta_for_ACT.csv"
                    new_act_df.to_csv(path_or_buf=df_loc, index=False)
                    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                    print("Rerunning act comparisons for %s isolates in cluster %s with new ref %s" % (
                    len(element_ids), cluster_name, new_ref_name))
                    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                    ## Now lets run these new acts to replace the current ones with the reference
                    act_command = "python " + python_dir + "/running_act_comparisons.py" + \
                                  " --csv " + df_loc + " --perl_dir " + perl_dir + "/" + " --act_dir ./act_compos/"

                    subprocess.call(act_command, shell=True)  # , stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                    refs_to_alter.append(ref_name)
                    new_refs.append(new_ref_name)
                    act_map = True
        else:
            act_map = True
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
        ref_contig = re.sub("#", "_", ref_name)
        ref_contig_file = contig_file_abs_path + ref_contig + contig_suffix
        contig_tab = pandas.read_csv(contig_file_path)
        ref_contig_tab = pandas.read_csv(ref_contig_file)

        contig_mge = contig_checker(contig_tab, hitters)

        mge_bounds = bounds_of_contig(contig_tab, contig_mge)

        ###########################################################################
        ## Now we'll do a check on the reference to see if it has the same ########
        ## identical hits as the current isolate in question ######################
        ###########################################################################

        hit_before, hit_after, overlap = ref_contains_hit(compo_table, hitters, mge_bounds, isolate_id)

        if overlap == "No":

            hit_before, hit_after, which_hit = before_and_after_hits(hitters, compo_table, mge_bounds,
                                                                     hits_to_search="both")

        else:
            which_hit = "both"

        if hit_before[0] == 0:
            contig_before = None
            hit_before_length = 0
        else:
            hit_before_loc = hit_before.iloc[[6, 7]]
            hit_before_sub_loc = hit_before.iloc[[8, 9]]
            hit_before_length = abs(hit_before_loc[1] - hit_before_loc[0])
            contig_before = contig_checker(contig_tab, hit_before_loc)
            hit_before_sub_loc = hit_before_sub_loc.sort_values(ascending=True)
            contig_bef_ref = contig_checker(ref_contig_tab, hit_before_sub_loc)

        if hit_after[0] == 0:
            contig_after = None
            hit_after_length = 0
        else:
            hit_after_loc = hit_after.iloc[[6, 7]]
            hit_after_sub_loc = hit_after.iloc[[8, 9]]
            hit_after_length = abs(hit_after_loc[1] - hit_after_loc[0])
            contig_after = contig_checker(contig_tab, hit_after_loc)
            hit_after_sub_loc = hit_after_sub_loc.sort_values(ascending=True)
            contig_aft_ref = contig_checker(ref_contig_tab, hit_after_sub_loc)

        all_one_tig = contig_before == contig_mge and contig_mge == contig_after and which_hit == "both"

        if all_one_tig and contig_bef_ref == contig_aft_ref and overlap == "No":

            hit_before_check = before_and_after_check(hit_before, hitters, compo_table, "before", hit_after, isolate_id)
            hit_after_check = before_and_after_check(hit_after, hitters, compo_table, "after", hit_before, isolate_id)

            if not hit_before_check.empty:
                hit_before = multi_hit_merger(hit_before_check)
                hit_before_loc = hit_before.iloc[[6, 7]]
                hit_before_length = abs(hit_before_loc[1] - hit_before_loc[0])

            if not hit_after_check.empty:
                hit_after = multi_hit_merger(hit_after_check)
                hit_after_loc = hit_after.iloc[[6, 7]]
                hit_after_length = abs(hit_after_loc[1] - hit_after_loc[0])

        all_one_tig_5k = hit_before_length >= 5000 and hit_after_length >= 5000 and all_one_tig

        if all_one_tig and not all_one_tig_5k:
            hit_before, hit_after, all_one_tig_5k = hit_mover(hit_before, hit_after, compo_table, isolate_id,
                                                              mge_bounds, mge_ori)
            hit_before_loc = hit_before.iloc[[6, 7]]
            hit_after_loc = hit_after.iloc[[6, 7]]

            if not all_one_tig_5k:

               hit_before, hit_after, all_one_tig_5k = final_acceptor(hit_before, hit_after, isolate_id, mge_bounds, mge_ori)
               hit_before_loc = hit_before.iloc[[6, 7]]
               hit_after_loc = hit_after.iloc[[6, 7]]



        if all_one_tig_5k:

            ## If the before and after hits to the isolate all line up on one contig we can assume this is
            ## a good enough hit to be used in library creation. Also needs to have at least 5k bp either side
            ## in this hit

            current_gff = pandas.read_csv(current_gff_loc.iloc[0], sep='\t',
                                          names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                                 'attributes'],
                                          header=None)

            if len(contig_tab.index) > 1:

                current_gff = gff_to_dna(current_gff, contig_tab, isolate_id, input_k=k)

            if mge_ori == "forward":
                ## Lets get the element length
                ## The insertion length
                ## Number of genes in the inserton
                ## number genes & average length 500bp either side.
                if current_row['merged'].values == "Yes":
                    current_mge_length = current_row['align']
                    mergio = True
                else:
                    current_mge_length = hitters[1] - hitters[0]
                    mergio = False


                current_insert_locs = [hit_before_loc[1], hit_after_loc[0]]
                current_insert_length = int(hit_after_loc[0]) - int(hit_before_loc[1])
                current_insert_s_locs = [hit_before.iloc[9], hit_after.iloc[8]]
                genes_insert = current_gff[
                    (current_gff['start'] >= current_insert_locs[0]) & (current_gff['end'] <= current_insert_locs[1])]
                genes_mge = current_gff[(current_gff['start'] >= hitters[0]) & (current_gff['end'] <= hitters[1])]
                genes_mge_num = len(genes_mge.index)
                gene_insert_num = len(genes_insert.index)

                ## 5000 bp regions.
                ## Include any genes overlapping region, so just based on end

                before_loc_gens = current_gff[(current_gff['end'] > (hit_before_loc[1] - 5000)) & (
                            current_gff['end'] <= (hit_before_loc[1] + 100))]
                num_genes_before = len(before_loc_gens.index)
                mean_length_before = [0]
                if not before_loc_gens.empty:
                    before_gene_lengths = []
                    for l in range(len(before_loc_gens.index)):
                        current_length = before_loc_gens.iloc[l, 4] - before_loc_gens.iloc[l, 3]
                        before_gene_lengths.append([current_length])
                    mean_length_before = numpy.mean(before_gene_lengths)

                after_loc_gens = current_gff[(current_gff['start'] >= (hit_after_loc[0] - 100)) & (
                            current_gff['start'] < (hit_after_loc[0] + 5000))]
                num_genes_after = len(after_loc_gens.index)
                mean_length_after = [0]
                if not after_loc_gens.empty:
                    after_loc_lengths = []
                    for l in range(len(after_loc_gens.index)):
                        current_length = after_loc_gens.iloc[l, 4] - after_loc_gens.iloc[l, 3]
                        after_loc_lengths.append([current_length])
                    mean_length_after = numpy.mean(after_loc_lengths)

                before_gene = gene_name_finder(before_loc_gens, back_it_up=True)
                after_gene = gene_name_finder(after_loc_gens, back_it_up=False)

                flanks_genes = num_genes_before + num_genes_after

                library_flank_gene_length = numpy.mean(before_gene_lengths + after_loc_lengths)
                # library_flank_gene_length = numpy.mean([mean_length_before, mean_length_after])
                tot_flanks_length = hit_before_length + hit_after_length

                if act_map:
                    ## Get the altered ref to use, then check if we can actually map this back to the reference in
                    ## the act_mapper function, if not skip hit
                    altered_index = [i for i, x in enumerate(refs_to_alter) if x == ref_name]
                    current_new = new_refs[altered_index[0]]
                    new_act_loc = absolute_act_path + current_new + ".crunch.gz"
                    if isolate_id != ref_name:
                        hit_before, hit_after, mapped = act_mapper(hit_before, hit_after, new_act_loc,
                                                                   current_insert_s_locs)
                        if not mapped:
                            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                            print("Can't map")
                            print(isolate_id)
                            print(hit_before_loc)
                            print(hit_after_loc)
                            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                            continue
                    else:
                        hit_before['sstart'] = hit_before['qstart']
                        hit_before['send'] = hit_before['qend']
                        hit_after['sstart'] = hit_after['qstart']
                        hit_after['send'] = hit_after['qend']
                library_pros = pandas.DataFrame()
                library_pros['id'] = pandas.Series(isolate_id)
                library_pros['mge_start'] = pandas.Series(hitters[0], index=library_pros.index)
                library_pros['mge_end'] = pandas.Series(hitters[1], index=library_pros.index)
                library_pros['insert_start'] = pandas.Series(current_insert_locs[0], index=library_pros.index)
                library_pros['insert_end'] = pandas.Series(current_insert_locs[1], index=library_pros.index)
                library_pros['before_sstart'] = pandas.Series(hit_before['sstart'], index=library_pros.index)
                library_pros['before_send'] = pandas.Series(hit_before['send'], index=library_pros.index)
                library_pros['after_sstart'] = pandas.Series(hit_after['sstart'], index=library_pros.index)
                library_pros['after_send'] = pandas.Series(hit_after['send'], index=library_pros.index)
                library_pros['mge_length'] = pandas.Series(current_mge_length, index=library_pros.index)
                library_pros['insert_length'] = pandas.Series(current_insert_length, index=library_pros.index)
                library_pros['insert_genes'] = pandas.Series(gene_insert_num, index=library_pros.index)
                library_pros['mge_genes'] = pandas.Series(genes_mge_num, index=library_pros.index)
                library_pros['flank_genes'] = pandas.Series(flanks_genes, index=library_pros.index)
                library_pros['mean_flank_gene_length'] = pandas.Series(library_flank_gene_length,
                                                                       index=library_pros.index)
                library_pros['flanks_length'] = pandas.Series(tot_flanks_length, index=library_pros.index)
                library_pros['before_flank_gene'] = pandas.Series(num_genes_before, index=library_pros.index)
                library_pros['after_flank_gene'] = pandas.Series(num_genes_after, index=library_pros.index)
                library_pros['before_flank_avg'] = pandas.Series(numpy.mean(before_gene_lengths),
                                                                 index=library_pros.index)
                library_pros['after_flank_avg'] = pandas.Series(numpy.mean(after_loc_lengths), index=library_pros.index)

                library_pros['before_gene_name'] = pandas.Series(before_gene, index=library_pros.index)
                library_pros['after_gene_name'] = pandas.Series(after_gene, index=library_pros.index)

                library_pros['ref_name'] = pandas.Series(ref_name, index=library_pros.index)
                library_pros['cluster_name'] = pandas.Series(cluster_name, index=library_pros.index)
                ## check if to add in to library csv

                if genes_mge_num <= gene_insert_num:
                    hit_df, missing_isolate = hit_detector(library_dat, library_pros, isolate_id, hit_df, missing_isolate, mergio)
                else:
                    missing_current = pandas.DataFrame()
                    missing_current['id'] = pandas.Series(isolate_id)
                    missing_current['mge_start'] = pandas.Series(hitters[0], index=missing_current.index)
                    missing_current['mge_end'] = pandas.Series(hitters[1], index=missing_current.index)
                    missing_current['ref_name'] = pandas.Series(ref_name, index=missing_current.index)
                    missing_current['cluster_name'] = pandas.Series(cluster_name, index=missing_current.index)
                    missing_current['mge_length'] = pandas.Series(current_mge_length, index=missing_current.index)
                    missing_current['reason'] = pandas.Series(["More mge than insert"], missing_current.index)
                    missing_isolate = missing_isolate.append(missing_current, sort=False)



            elif mge_ori == "reverse":

                if current_row['merged'].values == "Yes":
                    current_mge_length = current_row['align']
                    mergio = True
                else:
                    current_mge_length = hitters[0] - hitters[1]
                    mergio = False


                current_insert_locs = [hit_after_loc[1], hit_before_loc[0]]
                current_insert_length = int(hit_before_loc[0]) - int(hit_after_loc[1])
                current_insert_s_locs = [hit_after.iloc[9], hit_before.iloc[8]]
                genes_insert = current_gff[(current_gff['start'] >= current_insert_locs[0]) \
                                           & (current_gff['end'] <= current_insert_locs[1])]
                genes_mge = current_gff[(current_gff['start'] >= hitters[1]) & (current_gff['end'] <= hitters[0])]

                genes_mge_num = len(genes_mge.index)
                gene_insert_num = len(genes_insert.index)

                ## 5000 bp regions.
                ## Include any genes overlapping region, so just based on end
                before_loc_gens = current_gff[(current_gff['start'] > (hit_before_loc[0] - 100)) & \
                                              (current_gff['start'] <= (hit_before_loc[0] + 5000))]
                num_genes_before = len(before_loc_gens.index)
                mean_length_before = [0]
                if not before_loc_gens.empty:
                    before_gene_lengths = []
                    for l in range(len(before_loc_gens.index)):
                        current_length = before_loc_gens.iloc[l, 4] - before_loc_gens.iloc[l, 3]
                        before_gene_lengths.append([current_length])
                    mean_length_before = numpy.mean(before_gene_lengths)

                after_loc_gens = current_gff[(current_gff['end'] > (hit_after_loc[1] - 5000)) & \
                                             (current_gff['end'] <= (hit_after_loc[1] + 100))]
                num_genes_after = len(after_loc_gens.index)
                mean_length_after = [0]
                if not after_loc_gens.empty:
                    after_loc_lengths = []
                    for l in range(len(after_loc_gens.index)):
                        current_length = after_loc_gens.iloc[l, 4] - after_loc_gens.iloc[l, 3]
                        after_loc_lengths.append([current_length])
                    mean_length_after = numpy.mean(after_loc_lengths)

                flanks_genes = num_genes_before + num_genes_after

                before_gene = gene_name_finder(before_loc_gens, back_it_up=False)
                after_gene = gene_name_finder(after_loc_gens, back_it_up=True)

                library_flank_gene_length = numpy.mean([mean_length_before, mean_length_after])
                tot_flanks_length = hit_before_length + hit_after_length

                if act_map:
                    ## Get the altered ref to use, then check if we can actually map this back to the reference in
                    ## the act_mapper function, if not skip hit
                    altered_index = [i for i, x in enumerate(refs_to_alter) if x == ref_name]
                    current_new = new_refs[altered_index[0]]
                    new_act_loc = absolute_act_path + current_new + ".crunch.gz"
                    if isolate_id != ref_name:
                        hit_before, hit_after, mapped = act_mapper(hit_before, hit_after, new_act_loc,
                                                                   current_insert_s_locs)
                        if not mapped:
                            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                            print("Can't map")
                            print(isolate_id)
                            print(current_insert_s_locs)
                            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                            continue
                    else:
                        hit_before['sstart'] = hit_before['qstart']
                        hit_before['send'] = hit_before['qend']
                        hit_after['sstart'] = hit_after['qstart']
                        hit_after['send'] = hit_after['qend']
                library_pros = pandas.DataFrame()
                library_pros['id'] = pandas.Series(isolate_id)
                library_pros['mge_start'] = pandas.Series(hitters[0], index=library_pros.index)
                library_pros['mge_end'] = pandas.Series(hitters[1], index=library_pros.index)
                library_pros['insert_start'] = pandas.Series(current_insert_locs[0], index=library_pros.index)
                library_pros['insert_end'] = pandas.Series(current_insert_locs[1], index=library_pros.index)
                library_pros['before_sstart'] = pandas.Series(hit_before['sstart'], index=library_pros.index)
                library_pros['before_send'] = pandas.Series(hit_before['send'], index=library_pros.index)
                library_pros['after_sstart'] = pandas.Series(hit_after['sstart'], index=library_pros.index)
                library_pros['after_send'] = pandas.Series(hit_after['send'], index=library_pros.index)
                library_pros['mge_length'] = pandas.Series(current_mge_length, index=library_pros.index)
                library_pros['insert_length'] = pandas.Series(current_insert_length, index=library_pros.index)
                library_pros['insert_genes'] = pandas.Series(gene_insert_num, index=library_pros.index)
                library_pros['mge_genes'] = pandas.Series(genes_mge_num, index=library_pros.index)
                library_pros['flank_genes'] = pandas.Series(flanks_genes, index=library_pros.index)
                library_pros['mean_flank_gene_length'] = pandas.Series(library_flank_gene_length,
                                                                       index=library_pros.index)
                library_pros['flanks_length'] = pandas.Series(tot_flanks_length, index=library_pros.index)
                library_pros['before_flank_gene'] = pandas.Series(num_genes_before, index=library_pros.index)
                library_pros['after_flank_gene'] = pandas.Series(num_genes_after, index=library_pros.index)
                library_pros['before_flank_avg'] = pandas.Series(numpy.mean(before_gene_lengths),
                                                                 index=library_pros.index)
                library_pros['after_flank_avg'] = pandas.Series(numpy.mean(after_loc_lengths), index=library_pros.index)

                library_pros['before_gene_name'] = pandas.Series(before_gene, index=library_pros.index)
                library_pros['after_gene_name'] = pandas.Series(after_gene, index=library_pros.index)
                library_pros['ref_name'] = pandas.Series(ref_name, index=library_pros.index)
                library_pros['cluster_name'] = pandas.Series(cluster_name, index=library_pros.index)

                ## check if to add in to library csv
                if genes_mge_num <= gene_insert_num:
                    hit_df, missing_isolate = hit_detector(library_dat, library_pros, isolate_id, hit_df, missing_isolate, mergio)
                else:
                    missing_current = pandas.DataFrame()
                    missing_current['id'] = pandas.Series(isolate_id)
                    missing_current['mge_start'] = pandas.Series(hitters[0], index=missing_current.index)
                    missing_current['mge_end'] = pandas.Series(hitters[1], index=missing_current.index)
                    missing_current['ref_name'] = pandas.Series(ref_name, index=missing_current.index)
                    missing_current['cluster_name'] = pandas.Series(cluster_name, index=missing_current.index)
                    missing_current['mge_length'] = pandas.Series(current_mge_length, index=missing_current.index)
                    missing_current['reason'] = pandas.Series(["More mge than insert"], missing_current.index)
                    missing_isolate = missing_isolate.append(missing_current, sort=False)
        else:
            missing_current = pandas.DataFrame()
            missing_current['id'] = pandas.Series(isolate_id)
            missing_current['mge_start'] = pandas.Series(hitters[0], index=missing_current.index)
            missing_current['mge_end'] = pandas.Series(hitters[1], index=missing_current.index)
            missing_current['ref_name'] = pandas.Series(ref_name, index=missing_current.index)
            missing_current['cluster_name'] = pandas.Series(cluster_name, index=missing_current.index)
            if current_row['merged'].values[0] == "Yes":
                missing_current['mge_length'] = pandas.Series(current_row['align'], index=missing_current.index)
                missing_current['reason'] = pandas.Series(["MERGED No good hits before and after"],missing_current.index)
            else:
                missing_current['mge_length'] = pandas.Series(current_mge_length, index=missing_current.index)
                missing_current['reason'] = pandas.Series(["No good hits before and after"],
                                                          missing_current.index)
            missing_isolate = missing_isolate.append(missing_current, sort = False)
        end_len = len(missing_isolate.index) + len(hit_df.index)
        if end_len <= start_len:
            print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
            print(isolate_id)
            print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")


    hit_out_name = files_for_input.output + "_hits_df.csv"
    hit_df = hit_df.append(library_dat)

    missing_out_name = files_for_input.output + "_missing_df.csv"
    hit_df.to_csv(path_or_buf=hit_out_name, index=False)
    missing_isolate.to_csv(path_or_buf=missing_out_name, index = False)
    toc1 = time.perf_counter()
    toc = time.perf_counter()
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("This many hits allocated: %s of %s" % (len(hit_df.index), len(narrowed_prop.index) + len(library_dat.index)))
    print("This many missing hits: %s (%s %s)" % (len(missing_isolate.index), round((len(missing_isolate.index) / (len(narrowed_prop.index)+ len(library_dat.index))) * 100), "%"))
    print("Took this long for hit allocation: %s" % (toc1 - tic_1))
    print("Took this long overall: %s" % (toc - tic))
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hit allocator finished ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


