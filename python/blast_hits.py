import os
import pandas
from itertools import product
import numpy
import sys
import re
import difflib
import dendropy
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import subprocess
pandas.set_option('display.max_columns', 500)
pandas.set_option('display.width', 1000)
import time

def get_options():
    purpose = ''' This is a script to identify regions within isolates to search
    for with BLAST to detect the originaly progenitor of the sequence. 
    Usage: python blast_hit_identifiers.py <tree> <reccy_csv> <pyt_csv> <mut_per_branch_csv> <embl_reccys_to_csv> <act_compos> <flank_length>
    <reference_id> <fasta_directory> <out_directory> <out_csv>'''

    parser = argparse.ArgumentParser(description=purpose, prog='blast_hits.py')

    parser.add_argument('--gubbins_res', required=True,  help='Directory where all cluster dirs of gubbins res are stored"', type=str)
    parser.add_argument('--reccy_hits', required=True, help='csv from reccy finder with hits within recombinations', type=str)
    parser.add_argument('--hit_csv', required=True, help='hits csv out from hit_allocator', type=str)
    parser.add_argument('--contig_bounds', required=True, help='bounds of contigs used in reconstruction', type = str)
    parser.add_argument('--proper_hits',required=True, help='proper hits csv output from library creator', type=str)
    parser.add_argument('--act_compos', required=True, help='Location of act comparisons, with .referoo.fasta. prefix', type=str)
    parser.add_argument('--flank_length', required=True, help='Length to extract from flanks', type=int)
    parser.add_argument('--dna_dir', required=True, help='location of dna files', type=str)
    parser.add_argument('--out_dir', required=True, help='directory to save extracted flanks', type=str)
    parser.add_argument('--out_name', required=True, help='Prefix to append to out out_files', type=str)

    args = parser.parse_args()

    return args

def progress_bar():
    print("Narrowing down the isolates for blast searching", end="", flush=True)

def leaf_tips(tree, example_id, current_node):
    tip_node = tree.find_node_with_taxon_label(label=example_id)
    parent_node = tip_node.parent_node
    node_names = []




    if parent_node.label == current_node:
        tip_numbers = parent_node.leaf_nodes()
    else:
        while parent_node.label != current_node:
            parent_node = parent_node.parent_node

            if parent_node.label == current_node:
                tip_numbers = parent_node.leaf_nodes()

    for klose in tip_numbers:
        nodey_name = str(klose.taxon)
        noder_name = re.sub("\'","",nodey_name)
        node_name = re.sub("\'", "", noder_name)

        node_names.append(node_name)



    return node_names

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

def reference_present(insertion_node, tree, reference, isolate_examp):
    leaf_tippies = leaf_tips(tree,isolate_examp, insertion_node)

    if reference not in leaf_tippies:
        present = "No"
    else:
        present = "Yes"

    return  present

def all_presence_checker(id_list, series_ids):
    ## Function to check whether the total list of ids is within the blast csv
    ## Input: id_list: List of ids to check for presence all must be within the series ids
    ##        series_ids: the series of ids to check through
    skip = False
    num_in = 0
    for k in range(len(id_list)):
        if id_list[k] in series_ids:
            num_in += 1

    if num_in == len(id_list):
        skip = True


    return skip

def outside_control(insertion_node, tree, example_id, act_comp_dir, ref_insertion, nice_ids_tot):
    ## This function identifies control isolates outside of insertion node and the regions
    ## to look at

    outside_isolates, ultimate_node = chains_of_isolates_plus_one(tree,example_id, insertion_node)
    skip = all_presence_checker(outside_isolates, nice_ids_tot)

    while skip and ultimate_node != "internal_ROOT":
        outside_isolates, ultimate_node = chains_of_isolates_plus_one(tree,example_id,ultimate_node)
        skip = all_presence_checker(outside_isolates, nice_ids_tot)

    if skip == True:
        return 0, False, False



    lengths_of_izzys = {}

    for k in range(len(outside_isolates)):
        current_isolate = outside_isolates[k]
        chain_len = chains_of_isolates(tree, current_isolate, ultimate_node)
        lengths_of_izzys[current_isolate] = chain_len

    length_dict_sorted = {k: v for k, v in sorted(lengths_of_izzys.items(), key= lambda item: item[1])}
    finding_ref = True

    for ref in length_dict_sorted:
        if ref in nice_ids_tot:
            continue

        outside_iso = ref

        act_file = act_comp_dir + outside_iso + ".crunch.gz"
        compo_names = ['query', 'subject', 'pid', 'align', 'gap', 'mismatch', 'qstart',
                       'qend', 'sstart', 'send', 'eval', 'bitscore']
        compo_table = pandas.read_csv(act_file, sep='\t', names=compo_names)

        compo_table = compo_transformer(compo_table)



        forward_subset = compo_table[(compo_table['sstart'] <= ref_insertion[0]) & (compo_table['send'] >= ref_insertion[1])]
        if forward_subset.empty:
            reverse_subset = compo_table[(compo_table['send'] <= ref_insertion[0]) & (compo_table['sstart'] >= ref_insertion[1])]
            if reverse_subset.empty:
                before_forward = compo_table[(compo_table['sstart'] <= ref_insertion[0]) & (compo_table['send'] <= ref_insertion[0])]

                after_forward = compo_table[(compo_table['sstart'] >= ref_insertion[1]) & (compo_table['send'] >= ref_insertion[1])]

                ## INsert code here looking to arrange these values and then take the closest match to the start


                if after_forward.empty or before_forward.empty:
                    outside_iso = False
                    isolate_before_end = False
                    isolate_after_start = False

                    continue


                bef_forward_ordering = before_forward.sort_values(by=['send'], ascending=False)
                bef_reverse_ordering = before_forward.sort_values(by=['sstart'], ascending=False)
                bef_forward_send_row = bef_forward_ordering.iloc[0]
                bef_reverse_sstart_row = bef_reverse_ordering.iloc[0]
                bef_forward_send = bef_forward_send_row.iloc[9]
                bef_forward_sstart = bef_reverse_sstart_row.iloc[8]

                if bef_forward_send < bef_forward_sstart:
                    before_row = bef_reverse_sstart_row
                elif bef_forward_sstart < bef_forward_send:
                    before_row = bef_forward_send_row

                aft_reverse_ordering = after_forward.sort_values(by=['send'], ascending=True)
                aft_forward_ordering = after_forward.sort_values(by=['sstart'], ascending=True)
                aft_reverse_send_row = aft_reverse_ordering.iloc[0]
                aft_forward_sstart_row = aft_forward_ordering.iloc[0]
                aft_reverse_send = aft_reverse_send_row.iloc[9]
                aft_forward_sstart = aft_forward_sstart_row.iloc[8]

                if aft_reverse_send < aft_forward_sstart:
                    after_row = aft_reverse_send_row
                elif aft_forward_sstart < aft_reverse_send:
                    after_row = aft_forward_sstart_row




                isolate_before_end = before_row.iloc[7]
                isolate_after_start = after_row.iloc[6]


                break

            else:
                begin_diff = reverse_subset.iloc[0, 8] - ref_insertion[1]
                after_diff = ref_insertion[0] - reverse_subset.iloc[0, 9]

                isolate_before_end = forward_subset.iloc[0, 6] + begin_diff
                isolate_after_start = forward_subset.iloc[0, 7] - after_diff
                break



        else:

            begin_diff = ref_insertion[0] - forward_subset.iloc[0,8]
            after_diff = forward_subset.iloc[0,9] - ref_insertion[1]


            isolate_before_end = forward_subset.iloc[0, 6] + begin_diff
            isolate_after_start = forward_subset.iloc[0, 7] - after_diff
            break




    return outside_iso, isolate_before_end, isolate_after_start

def compo_transformer(comps_tab):
    for k in range(len(comps_tab.index)):
        current_sstart = comps_tab.iloc[k, 8]
        current_send = comps_tab.iloc[k, 9]
        if current_sstart > current_send:
            comps_tab.iloc[k, 8] = current_send
            comps_tab.iloc[k, 9] = current_sstart


    return comps_tab

def act_bounds(python_row):
    before_regions = python_row[['before_sstart','before_send']]
    after_regions = python_row[['after_sstart','after_send']]

    midpoint_bef = (before_regions.iloc[0,1] + before_regions.iloc[0,0]) / 2
    midpoint_aft = (after_regions.iloc[0,1] + after_regions.iloc[0,0]) / 2
    if midpoint_bef < midpoint_aft:
        closest_vals = [before_regions.min(axis=1).values[0], after_regions.max(axis=1).values[0]]
    elif midpoint_bef > midpoint_aft:
        closest_vals = [after_regions.min(axis=1).values[0], before_regions.max(axis=1).values[0]]


    return closest_vals

def chains_of_isolates(tree, example_id, end_node):
    tip_node = tree.find_node_with_taxon_label(label=example_id)
    parent_node = tip_node.parent_node
    node_names = []
    node_names.append(example_id)

    if parent_node.label == end_node:
        node_names.append(parent_node.label)

    else:
        while parent_node.label != end_node:
            node_names.append(parent_node.label)
            parent_node = parent_node.parent_node


    return node_names

def chains_of_isolates_plus_one(tree, example_id, end_node):
    tip_node = tree.find_node_with_taxon_label(label=example_id)
    parent_node = tip_node.parent_node
    node_names = []
    node_names.append(example_id)

    if parent_node.label == end_node:
        node_names.append(parent_node.label)
        ultimate_node = parent_node.parent_node

    else:
        while parent_node.label != end_node:
            node_names.append(parent_node.label)
            parent_node = parent_node.parent_node
            if parent_node.label  == end_node:
                ultimate_node = parent_node.parent_node

    insertion_node_tips = parent_node.leaf_nodes()
    ultimate_node_tips = ultimate_node.leaf_nodes()

    inny_nodes = Leaf_nodes_to_list(insertion_node_tips)
    ulty_nodes = Leaf_nodes_to_list(ultimate_node_tips)


    poss_nodes = numpy.setdiff1d(ulty_nodes, inny_nodes)



    return poss_nodes, ultimate_node.label

def Leaf_nodes_to_list(node_list):
    node_names = []
    for klose in node_list:
        nodey_name = str(klose.taxon)
        noder_name = re.sub("\'","",nodey_name)
        node_name = re.sub("\'", "", noder_name)

        node_names.append(node_name)

    return  node_names

def largest_reccombination_block(before_row, after_row):
    before_length = before_row.iloc[0,3]
    after_length = after_row.iloc[0,3]

    if before_length > after_length:
        start_return = before_row.iloc[0,6]
        end_return = before_row.iloc[0,7]
    elif after_length > before_length:
        start_return = after_row.iloc[0, 6]
        end_return = after_row.iloc[0, 7]

    return start_return, end_return

def compo_enlarger(start_act, ori, pos, act_compo, target_length, debug_id, pre_length):
    current_length = abs(start_act.iloc[0,6] - start_act.iloc[0, 7])
    current_gap =  current_length



    gap_list = [[start_act.iloc[0, 6],start_act.iloc[0,7]]]

    while current_gap < target_length:
        if ori == "forward":
            if pos == "bef":
                current_ref_pos = [start_act.iloc[0,8], start_act.iloc[0, 9]]
                if current_ref_pos[0] < current_ref_pos[1]:
                    next_act_event_send = act_compo[act_compo['send'] < current_ref_pos[0]]
                    next_act_event_send = next_act_event_send.sort_values(by = ['send'], ascending=False)
                    next_act_event_send_row = next_act_event_send.iloc[0]

                    ref_start_send = next_act_event_send_row.iloc[9]



                    next_act_event_sstart = act_compo[act_compo['sstart'] < current_ref_pos[0]]
                    next_act_event_sstart = next_act_event_sstart.sort_values(by = ['sstart'], ascending=False)
                    next_act_event_sstart_row = next_act_event_sstart.iloc[0]

                    ref_start_sstart = next_act_event_sstart_row.iloc[8]

                    if ref_start_sstart == ref_start_send:
                        alignos = [next_act_event_sstart_row.iloc[3], next_act_event_send_row.iloc[3]]
                        max_indy = alignos.index(max(alignos))
                        if max_indy == 0:
                            ref_start_sstart += 1
                        elif max_indy == 1:
                            ref_start_send += 1
                        else:
                            ref_start_send += 1


                    if ref_start_send >= ref_start_sstart:
                        distance_hits = abs(next_act_event_send_row.iloc[7] - next_act_event_send_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_send_row.iloc[6] + (current_gap - target_length), next_act_event_send_row.iloc[7]]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                            gap_list.append(gappers)
                    elif ref_start_send < ref_start_sstart:
                        distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_sstart_row.iloc[6] + (current_gap - target_length),
                                       next_act_event_sstart_row.iloc[7]]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                            gap_list.append(gappers)
                elif current_ref_pos[1] < current_ref_pos[0]:
                    next_act_event_send = act_compo[act_compo['send'] < current_ref_pos[1]]
                    next_act_event_send = next_act_event_send.sort_values(by=['send'], ascending=False)
                    next_act_event_send_row = next_act_event_send.iloc[0]

                    ref_start_send = next_act_event_send_row.iloc[9]

                    next_act_event_sstart = act_compo[act_compo['sstart'] < current_ref_pos[1]]
                    next_act_event_sstart = next_act_event_sstart.sort_values(by=['sstart'], ascending=False)
                    next_act_event_sstart_row = next_act_event_sstart.iloc[0]

                    ref_start_sstart = next_act_event_sstart_row.iloc[8]

                    if ref_start_sstart == ref_start_send:
                        alignos = [next_act_event_sstart_row.iloc[3], next_act_event_send_row.iloc[3]]
                        max_indy = alignos.index(max(alignos))
                        if max_indy == 0:
                            ref_start_sstart += 1
                        elif max_indy == 1:
                            ref_start_send += 1
                        else:
                            ref_start_send += 1


                    if ref_start_send >= ref_start_sstart:
                        distance_hits = abs(next_act_event_send_row.iloc[7] - next_act_event_send_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_send_row.iloc[6] + (current_gap - target_length),
                                       next_act_event_send_row.iloc[7]]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                            gap_list.append(gappers)
                    elif ref_start_send < ref_start_sstart:

                        distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_sstart_row.iloc[6] + (current_gap - target_length),
                                       next_act_event_sstart_row.iloc[7]]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                            gap_list.append(gappers)

            elif pos == "aft":
                current_ref_pos = [start_act.iloc[0, 8], start_act.iloc[0, 9]]
                if current_ref_pos[0] < current_ref_pos[1]:
                    next_act_event_send = act_compo[act_compo['send'] > current_ref_pos[1]]
                    next_act_event_send = next_act_event_send.sort_values(by=['send'], ascending=True)
                    next_act_event_send_row = next_act_event_send.iloc[0]

                    ref_start_send = next_act_event_send_row.iloc[9]

                    next_act_event_sstart = act_compo[act_compo['sstart'] > current_ref_pos[1]]
                    next_act_event_sstart = next_act_event_sstart.sort_values(by=['sstart'], ascending=True)
                    next_act_event_sstart_row = next_act_event_sstart.iloc[0]

                    ref_start_sstart = next_act_event_sstart_row.iloc[8]

                    if ref_start_sstart == ref_start_send:
                        alignos = [next_act_event_sstart_row.iloc[3], next_act_event_send_row.iloc[3]]
                        max_indy = alignos.index(max(alignos))
                        if max_indy == 0:
                            ref_start_sstart -= 1
                        elif max_indy == 1:
                            ref_start_send -= 1
                        else:
                            ref_start_send -= 1


                    if ref_start_send < ref_start_sstart:
                        distance_hits = abs(next_act_event_send_row.iloc[7] - next_act_event_send_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_send_row.iloc[6] ,
                                       next_act_event_send_row.iloc[7] - (current_gap - target_length)]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                            gap_list.append(gappers)
                    elif ref_start_send >= ref_start_sstart:
                        distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_sstart_row.iloc[6],
                                       next_act_event_sstart_row.iloc[7] - (current_gap - target_length)]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                            gap_list.append(gappers)
                elif current_ref_pos[1] < current_ref_pos[0]:
                    next_act_event_send = act_compo[act_compo['send'] > current_ref_pos[0]]
                    next_act_event_send = next_act_event_send.sort_values(by=['send'], ascending=True)
                    next_act_event_send_row = next_act_event_send.iloc[0]

                    ref_start_send = next_act_event_send_row.iloc[9]

                    next_act_event_sstart = act_compo[act_compo['sstart'] > current_ref_pos[0]]
                    next_act_event_sstart = next_act_event_sstart.sort_values(by=['sstart'], ascending=True)
                    next_act_event_sstart_row = next_act_event_sstart.iloc[0]

                    ref_start_sstart = next_act_event_sstart_row.iloc[8]

                    if ref_start_sstart == ref_start_send:
                        alignos = [next_act_event_sstart_row.iloc[3], next_act_event_send_row.iloc[3]]
                        max_indy = alignos.index(max(alignos))
                        if max_indy == 0:
                            ref_start_sstart -= 1
                        elif max_indy == 1:
                            ref_start_send -= 1
                        else:
                            ref_start_send -= 1


                    if ref_start_send < ref_start_sstart:
                        distance_hits = abs(next_act_event_send_row.iloc[7] - next_act_event_send_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_send_row.iloc[6] ,
                                       next_act_event_send_row.iloc[7] - (current_gap - target_length)]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                            gap_list.append(gappers)
                    elif ref_start_send >= ref_start_sstart:
                        distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_sstart_row.iloc[6] ,
                                       next_act_event_sstart_row.iloc[7] - (current_gap - target_length)]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                            gap_list.append(gappers)
        elif ori == "reverse":
            if pos == "bef":
                current_ref_pos = [start_act.iloc[0,8], start_act.iloc[0, 9]]
                if current_ref_pos[0] < current_ref_pos[1]:
                    next_act_event_send = act_compo[act_compo['send'] > current_ref_pos[1]]
                    next_act_event_send = next_act_event_send.sort_values(by = ['send'], ascending=True)
                    next_act_event_send_row = next_act_event_send.iloc[0]

                    ref_start_send = next_act_event_send_row.iloc[9]



                    next_act_event_sstart = act_compo[act_compo['sstart'] > current_ref_pos[1]]
                    next_act_event_sstart = next_act_event_sstart.sort_values(by = ['sstart'], ascending=True)
                    next_act_event_sstart_row = next_act_event_sstart.iloc[0]

                    ref_start_sstart = next_act_event_sstart_row.iloc[8]

                    if ref_start_sstart == ref_start_send:
                        alignos = [next_act_event_sstart_row.iloc[3], next_act_event_send_row.iloc[3]]
                        max_indy = alignos.index(max(alignos))
                        if max_indy == 0:
                            ref_start_sstart -= 1
                        elif max_indy == 1:
                            ref_start_send -= 1
                        else:
                            ref_start_send -= 1



                    if ref_start_send < ref_start_sstart:
                        distance_hits = abs(next_act_event_send_row.iloc[7] - next_act_event_send_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [ next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7] - (current_gap - target_length)]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                            gap_list.append(gappers)
                    elif ref_start_send >= ref_start_sstart:

                        distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_sstart_row.iloc[6] ,
                                       next_act_event_sstart_row.iloc[7] - (current_gap - target_length)]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                            gap_list.append(gappers)
                elif current_ref_pos[1] < current_ref_pos[0]:
                    next_act_event_send = act_compo[act_compo['send'] > current_ref_pos[0]]
                    next_act_event_send = next_act_event_send.sort_values(by=['send'], ascending=True)
                    next_act_event_send_row = next_act_event_send.iloc[0]

                    ref_start_send = next_act_event_send_row.iloc[9]

                    next_act_event_sstart = act_compo[act_compo['sstart'] > current_ref_pos[0]]
                    next_act_event_sstart = next_act_event_sstart.sort_values(by=['sstart'], ascending=True)
                    next_act_event_sstart_row = next_act_event_sstart.iloc[0]

                    ref_start_sstart = next_act_event_sstart_row.iloc[8]

                    if ref_start_sstart == ref_start_send:
                        alignos = [next_act_event_sstart_row.iloc[3], next_act_event_send_row.iloc[3]]
                        max_indy = alignos.index(max(alignos))
                        if max_indy == 0:
                            ref_start_sstart -= 1
                        elif max_indy == 1:
                            ref_start_send -= 1
                        else:
                            ref_start_send -= 1

                    if ref_start_send < ref_start_sstart:
                        distance_hits = abs(next_act_event_send_row.iloc[7] - next_act_event_send_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_send_row.iloc[6] ,
                                       next_act_event_send_row.iloc[7] - (current_gap - target_length)]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                            gap_list.append(gappers)
                    elif ref_start_send > ref_start_sstart:
                        distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_sstart_row.iloc[6] ,
                                       next_act_event_sstart_row.iloc[7] - (current_gap - target_length)]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                            gap_list.append(gappers)

            elif pos == "aft":
                current_ref_pos = [start_act.iloc[0, 8], start_act.iloc[0, 9]]




                if current_ref_pos[0] < current_ref_pos[1]:
                    next_act_event_send = act_compo[act_compo['send'] < current_ref_pos[0]]
                    next_act_event_send = next_act_event_send.sort_values(by=['send'], ascending=False)
                    next_act_event_send_row = next_act_event_send.iloc[0]

                    ref_start_send = next_act_event_send_row.iloc[9]

                    next_act_event_sstart = act_compo[act_compo['sstart'] < current_ref_pos[0]]
                    next_act_event_sstart = next_act_event_sstart.sort_values(by=['sstart'], ascending=False)
                    next_act_event_sstart_row = next_act_event_sstart.iloc[0]

                    ref_start_sstart = next_act_event_sstart_row.iloc[8]

                    if ref_start_sstart == ref_start_send:

                        alignos = [next_act_event_sstart_row.iloc[3], next_act_event_send_row.iloc[3]]
                        max_indy = alignos.index(max(alignos))
                        if max_indy == 0:
                            ref_start_sstart += 1
                        elif max_indy == 1:
                            ref_start_send += 1
                        else:
                            ref_start_send += 1

                    if ref_start_send > ref_start_sstart:
                        distance_hits = abs(next_act_event_send_row.iloc[7] - next_act_event_send_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_send_row.iloc[6] + (current_gap - target_length) ,
                                       next_act_event_send_row.iloc[7] ]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                            gap_list.append(gappers)
                    elif ref_start_send < ref_start_sstart:
                        distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_sstart_row.iloc[6] + (current_gap - target_length),
                                       next_act_event_sstart_row.iloc[7] ]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                            gap_list.append(gappers)
                elif current_ref_pos[1] < current_ref_pos[0]:
                    next_act_event_send = act_compo[act_compo['send'] < current_ref_pos[1]]


                    next_act_event_send = next_act_event_send.sort_values(by=['send'], ascending=False)
                    next_act_event_send_row = next_act_event_send.iloc[0]

                    ref_start_send = next_act_event_send_row.iloc[9]



                    next_act_event_sstart = act_compo[act_compo['sstart'] < current_ref_pos[1]]


                    next_act_event_sstart = next_act_event_sstart.sort_values(by=['sstart'], ascending=False)
                    next_act_event_sstart_row = next_act_event_sstart.iloc[0]



                    ref_start_sstart = next_act_event_sstart_row.iloc[8]

                    if ref_start_sstart == ref_start_send:

                        alignos = [next_act_event_sstart_row.iloc[3], next_act_event_send_row.iloc[3]]
                        max_indy = alignos.index(max(alignos))
                        if max_indy == 0:
                            ref_start_sstart += 1
                        elif max_indy == 1:
                            ref_start_send += 1
                        else:
                            ref_start_send += 1



                    if ref_start_send > ref_start_sstart:
                        distance_hits = abs(next_act_event_send_row.iloc[7] - next_act_event_send_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_send_row.iloc[6] + (current_gap - target_length),
                                       next_act_event_send_row.iloc[7] ]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                            gap_list.append(gappers)
                    elif ref_start_send < ref_start_sstart:
                        distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                        current_gap = current_gap + distance_hits
                        if current_gap >= target_length:
                            gappers = [next_act_event_sstart_row.iloc[6] + (current_gap - target_length),
                                       next_act_event_sstart_row.iloc[7]]
                            gap_list.append(gappers)
                        else:
                            gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                            gap_list.append(gappers)

    return gap_list

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

def isolate_narrow(reccy_hits, pyt_csv, tree, reccy_csv_gubbins, mut_bases_csv, reference_id, flanking_length, contig_loc, nice_ids_tot):
    ## Function to take in the reccy hits csv for a paticular cluster. Then work through the
    ## reccy hits and find the isolate with the fewest snps around the insertion site. This will
    ## then be output in the out_df along with the end and start points for the flanks to be taken
    ## from. Also finds a suitable reference isolate close to the clade, if the actual reference has
    ## the element present within.
    ## Input: reccy_hits: The subset of the reccy hits csv for a particular cluster
    ##        pyt_csv: Subset of the hit allocation csv for this particular cluster
    ##        tree: The node labelled tree for this cluster
    ##        reccy_csv_gubbins: The embl recombination predictions to csv for this cluster
    ##        mut_bases_csv: The mutations along branches csv for this cluster
    ##        reference_id: The reference for this cluster

    isolate_id = []
    mge_bef = []
    mge_aft = []
    before_start_region = []
    before_end_region = []
    after_start_region = []
    after_end_region = []
    snp_count = []
    control_id = []
    control_start = []
    control_end = []
    insert_start = []
    insert_end = []

    reccy_starts_bef_control = []
    reccy_ends_bef_control = []
    reccy_starts_aft_control = []
    reccy_ends_aft_control = []

    regions_bef = []
    regions_aft = []

    out_df = pandas.DataFrame()



    for k in range(len(reccy_hits.index)):
        current_row = reccy_hits.iloc[k]
        isolate_row = current_row['isolate_example']
        insertion_node = current_row['insertion_node']
        start_node_to_reccy = current_row['start_nodes']
        end_node_to_reccy = current_row['end_nodes']
        ref_insert_start = current_row['start_insert']
        ref_insert_end = current_row['end_insert']
        current_insert_name = current_row['insert_name']

        ref_insert_list = [ref_insert_start, ref_insert_end]

        starting_id = isolate_row
        nodes_from_insertion = leaf_tips(tree, example_id=isolate_row, current_node=insertion_node)

        expansio = []
        mge_isolates = []


        ## Find the biggest act span isolate from the reccy hits csv. This finds all the leaf tips from
        ## the reconstructed insertion node. If an isolate is present which doesn't have the mge then
        ## we ignore it.

        for l in range(len(nodes_from_insertion)):
            august_isolate = nodes_from_insertion[l]
            if pyt_csv['id'].isin([august_isolate]).any():
                pyt_row = pyt_csv[pyt_csv['id'] == august_isolate]

                if pyt_row['insert_name'].isin([current_insert_name]).any():
                    starts = pyt_row[['before_sstart', 'before_send']]
                    ends = pyt_row[['after_sstart', 'after_send']]
                    starts_span = abs(starts.iloc[0, 1] - starts.iloc[0, 0])
                    ends_span = abs(ends.iloc[0, 1] - ends.iloc[0, 0])
                    total_span = starts_span + ends_span
                    expansio.append(total_span)
                    mge_isolates.append(august_isolate)
                else:

                    expansio.append(1)

            else:
                expansio.append(1)

        max_span_index = expansio.index(max(expansio))
        max_span_isolate = nodes_from_insertion[max_span_index]

        act_edges = pyt_csv[pyt_csv['id'] == max_span_isolate]
        act_bounders = act_bounds(act_edges)

        total_life_forever = []

        ## Go from the current isolate to the insertion loci, getting the mutations along the branches
        ## Finds this among all the isolates that have the mge presen from the insertion node.

        for m in range(len(mge_isolates)):
            current_isolate = mge_isolates[m]
            current_chain = chains_of_isolates(tree, current_isolate, start_node_to_reccy)
            branch_muts = 0
            reccy_muts = 0
            node_rec = current_isolate in end_node_to_reccy

            if node_rec == False:
                for n in range(len(current_chain)):

                    current_end = current_chain[n]
                    reccy_detected = reccy_csv_gubbins[reccy_csv_gubbins['end_node'] == current_end]
                    reccy_hits_loc = reccy_detected[((reccy_detected['start_rec'] >= act_bounders[0]) & (
                                reccy_detected['start_rec'] <= act_bounders[1])) | (
                                                                (reccy_detected['end_rec'] >= act_bounders[0]) & (
                                                                    reccy_detected['end_rec'] <= act_bounders[1]))]
                    if not reccy_hits_loc.empty:
                        reccy_muts += reccy_hits_loc['snp_number'].sum()

                    mut_changes_sub = mut_bases_csv[mut_bases_csv['end_node'] == current_end]
                    mut_changes = mut_changes_sub.iloc[0, 3:26].sum()
                    branch_muts += mut_changes
                    total_changes = reccy_muts + branch_muts

                total_life_forever.append(total_changes)
            elif node_rec == True:
                current_end = current_isolate
                reccy_detected = reccy_csv_gubbins[reccy_csv_gubbins['end_node'] == current_end]
                reccy_hits_loc = reccy_detected[((reccy_detected['start_rec'] >= act_bounders[0]) & (
                        reccy_detected['start_rec'] <= act_bounders[1])) | (
                                                        (reccy_detected['end_rec'] >= act_bounders[0]) & (
                                                        reccy_detected['end_rec'] <= act_bounders[1]))]
                if not reccy_hits_loc.empty:
                    reccy_muts += reccy_hits_loc['snp_number'].sum()

                mut_changes_sub = mut_bases_csv[mut_bases_csv['end_node'] == current_end]
                mut_changes = mut_changes_sub.iloc[0, 3:26].sum()
                branch_muts += mut_changes
                total_life_forever.append(reccy_muts + branch_muts)



        mge_seq_to_look_at = total_life_forever.index(min(total_life_forever))
        snp_count_indiv = min(total_life_forever)
        mge_id = mge_isolates[mge_seq_to_look_at]
        orig_mge_id = mge_id
        orig_snp_count = snp_count_indiv
        while mge_id in isolate_id:
            total_life_forever.remove(snp_count_indiv)
            mge_isolates.remove(mge_id)
            if len(mge_isolates) == 0:
                mge_id = orig_mge_id
                snp_count_indiv = orig_snp_count
                break
            mge_seq_to_look_at = total_life_forever.index(min(total_life_forever))
            snp_count_indiv = min(total_life_forever)
            mge_id = mge_isolates[mge_seq_to_look_at]



        mge_deets = pyt_csv[pyt_csv['id'] == mge_id]


        ## Check if reference among the tips for this node insertion.

        reference_presence = reference_present(insertion_node=insertion_node,
                                               tree=tree,
                                               isolate_examp=isolate_row,
                                               reference=reference_id)

        if reference_presence == "No" and reference_id not in nice_ids_tot:
            csv_ref_name = "ref" + "!" + mge_id

            control_id.append(csv_ref_name)
            control_start.append(ref_insert_start)
            control_end.append(ref_insert_end)

            if mge_deets.iloc[0, 1] < mge_deets.iloc[0, 2]:
                rec_end_bef = ref_insert_start
                rec_start_aft = ref_insert_end

                rec_start_bef = ref_insert_start - flanking_length
                rec_end_aft = ref_insert_end + flanking_length

                reccy_starts_bef_control.append(rec_start_bef)
                reccy_ends_bef_control.append(rec_end_bef)
                reccy_starts_aft_control.append(rec_start_aft)
                reccy_ends_aft_control.append(rec_end_aft)

            elif mge_deets.iloc[0, 1] > mge_deets.iloc[0, 2]:
                rec_start_bef = ref_insert_end
                rec_end_aft = ref_insert_start

                rec_end_bef = rec_start_bef + flanking_length
                rec_start_aft = rec_end_aft - flanking_length

                reccy_starts_bef_control.append(rec_start_bef)
                reccy_ends_bef_control.append(rec_end_bef)
                reccy_starts_aft_control.append(rec_start_aft)
                reccy_ends_aft_control.append(rec_end_aft)




        elif reference_presence == "Yes" or reference_id in nice_ids_tot:
            cont_id, cont_start, cont_end = outside_control(insertion_node, tree, isolate_row, act_compos,
                                                            ref_insert_list, nice_ids_tot)
            if isinstance(cont_id, int):
                print("Can't extract control for this isolate, skipping: %s" % mge_id)
                continue

            csv_ref_name = cont_id + "!" + mge_id

            control_id.append(csv_ref_name)

            if mge_deets.iloc[0, 1] < mge_deets.iloc[0, 2]:
                rec_end_bef = cont_start
                rec_start_aft = cont_end

                rec_start_bef = rec_end_bef - flanking_length
                rec_end_aft = rec_start_aft + flanking_length

                reccy_starts_bef_control.append(rec_start_bef)
                reccy_ends_bef_control.append(rec_end_bef)
                reccy_starts_aft_control.append(rec_start_aft)
                reccy_ends_aft_control.append(rec_end_aft)

            elif mge_deets.iloc[0, 1] > mge_deets.iloc[0, 2]:
                rec_start_bef = cont_end
                rec_end_aft = cont_start

                rec_end_bef = rec_start_bef + flanking_length



                rec_start_aft = rec_end_aft - flanking_length

                reccy_starts_bef_control.append(rec_start_bef)
                reccy_ends_bef_control.append(rec_end_bef)
                reccy_starts_aft_control.append(rec_start_aft)
                reccy_ends_aft_control.append(rec_end_aft)

        isolate_id.append(mge_id)
        mge_bef.append(mge_deets.iloc[0, 1])
        mge_aft.append(mge_deets.iloc[0, 2])
        before_start_region.append(mge_deets.iloc[0, 5])
        before_end_region.append(mge_deets.iloc[0, 6])
        after_start_region.append(mge_deets.iloc[0, 7])
        after_end_region.append(mge_deets.iloc[0, 8])
        snp_count.append(snp_count_indiv)
        insert_start.append(mge_deets.iloc[0, 3])
        insert_end.append(mge_deets.iloc[0, 4])

    if len(isolate_id) == 0:
        return  "no", "no", "no"


    out_df['isolate_id'] = pandas.Series(data=isolate_id)
    out_df['mge_start'] = pandas.Series(data=mge_bef, index=out_df.index)
    out_df['mge_end'] = pandas.Series(data=mge_aft, index=out_df.index)
    out_df['before_start'] = pandas.Series(data=before_start_region, index=out_df.index)
    out_df['before_end'] = pandas.Series(data=before_end_region, index=out_df.index)
    out_df['after_start'] = pandas.Series(data=after_start_region, index=out_df.index)
    out_df['after_end'] = pandas.Series(data=after_end_region, index=out_df.index)
    out_df['snp_count'] = pandas.Series(data=snp_count, index=out_df.index)
    out_df['insert_start'] = pandas.Series(insert_start, index=out_df.index)
    out_df['insert_end'] = pandas.Series(insert_end, index=out_df.index)

    reccy_starts_bef = []
    reccy_ends_bef = []
    reccy_starts_aft = []
    reccy_ends_aft = []

    ## Now lets work through this out df to get the

    for o in range(len(out_df.index)):
        current_isolate = out_df.iloc[o, 0]
        mge_hits = out_df.iloc[o, [1,2]]
        bef_hits = out_df.iloc[o, [3,4]]
        aft_hits = out_df.iloc[o, [5,6]]
        insert_hits = out_df.iloc[o, [8,9]]
        act_file = act_compos + current_isolate + ".crunch.gz"
        compo_names = ['query', 'subject', 'pid', 'align', 'gap', 'mismatch', 'qstart',
                       'qend', 'sstart', 'send', 'eval', 'bitscore']

        contig_suffix = "#contig_bounds.csv"
        contig_isolate = re.sub("#", "_", current_isolate)
        contig_file_path = contig_loc + contig_isolate + contig_suffix

        contig_file = pandas.read_csv(contig_file_path)


        compo_table = pandas.read_csv(act_file, sep='\t', names=compo_names)

        before_subset = compo_table[(compo_table['sstart'] == bef_hits[0]) & (compo_table['send'] == bef_hits[1])]
        after_subset = compo_table[(compo_table['sstart'] == aft_hits[0]) & (compo_table['send'] == aft_hits[1])]

        ###########################################################################
        ## This is for when the mge hit is in the reference and the whole region ##
        ## will be in the act compo region ########################################
        ###########################################################################

        if before_subset.empty:
            if mge_hits[0] < mge_hits[1]:
                before_subset = compo_table[(compo_table['sstart'] == bef_hits[0]) & (compo_table['send'] == aft_hits[1])]
                after_subset = compo_table[(compo_table['sstart'] == bef_hits[0]) & (compo_table['send'] == aft_hits[1])]
            elif mge_hits[0] > mge_hits[1]:
                before_subset = compo_table[(compo_table['sstart'] == aft_hits[0]) & (compo_table['send'] == bef_hits[1])]
                after_subset = compo_table[(compo_table['sstart'] == aft_hits[0]) & (compo_table['send'] == bef_hits[1])]
        ## It might be the hit info for this isolate comes from combined act hits, so we need to find the closest hits
        ## in together then take the ACT compos from these closest hits in.
        if before_subset.empty or after_subset.empty:

            if mge_hits[0] < mge_hits[1]:
                before_subset = compo_table[compo_table['qend'] == insert_hits[0]]
                after_subset = compo_table[compo_table['qstart'] == insert_hits[1]]
                before_subset = before_subset.sort_values(by='align', ascending=False)
                after_subset = after_subset.sort_values(by='align', ascending=False)
            elif mge_hits[0] > mge_hits[1]:
                before_subset = compo_table[compo_table['qstart'] == insert_hits[1]]
                after_subset = compo_table[compo_table['qend'] == insert_hits[0]]
                before_subset = before_subset.sort_values(by='align', ascending=False)
                after_subset = after_subset.sort_values(by='align', ascending=False)

        ## Need to check if overlap within hit but not throughout hit.

        if before_subset.empty:
             if mge_hits[0] < mge_hits[1]:
                 before_subset = compo_table[compo_table['qend'] <= insert_hits[1]]
                 before_subset = before_subset.sort_values(by='qend', ascending=False)
                 if before_subset.iloc[0,7] > mge_hits[0]:
                     before_subset.iloc[0,7] = mge_hits[0]
             elif mge_hits[0] > mge_hits[1]:
                 before_subset = compo_table[compo_table['qstart'] >= insert_hits[0]]
                 before_subset = before_subset.sort_values(by='qstart', ascending=True)
                 if before_subset.iloc[0, 6] < mge_hits[0]:
                     before_subset.iloc[0,6] = mge_hits[0]

        if after_subset.empty:
             if mge_hits[0] < mge_hits[1]:
                 after_subset = compo_table[compo_table['qstart'] >= insert_hits[0]]
                 after_subset = after_subset.sort_values(by='qstart', ascending=True)
                 if after_subset.iloc[0,6] < mge_hits[1]:
                     after_subset.iloc[0,6] = mge_hits[1]
             elif mge_hits[0] > mge_hits[1]:
                after_subset = compo_table[compo_table['qend'] <= insert_hits[1]]
                after_subset = after_subset.sort_values(by='qend', ascending=False)
                if after_subset.iloc[0,7] > mge_hits[1]:
                    after_subset.iloc[0,7] = mge_hits[1]


        out_df.at[o, 'before_start'] = before_subset.iloc[0, 6]
        out_df.at[o, 'before_end'] = before_subset.iloc[0, 7]
        out_df.at[o, 'after_start'] = after_subset.iloc[0, 6]
        out_df.at[o, 'after_end'] = after_subset.iloc[0, 7]

        ###########################################################################
        ## Lets take 20k around the hit rather than 20k from start ################
        ###########################################################################




        if mge_hits[0] < mge_hits[1]:
            rec_start_bef = before_subset.iloc[0,6]
            rec_end_bef = before_subset.iloc[0,7]
            bef_hits = [rec_end_bef - flanking_length, rec_end_bef]
            rec_start_aft = after_subset.iloc[0,6]
            rec_end_aft = after_subset.iloc[0,7]
            aft_hits = [rec_start_aft ,rec_start_aft + flanking_length]
            bef_contig = contig_checker(contig_file, bef_hits)
            aft_contig = contig_checker(contig_file, aft_hits)



            if bef_contig == 0:
                current_contig_bounds = bounds_of_contig(contig_file,contig_checker(contig_file, mge_hits))
                current_length = rec_end_bef - current_contig_bounds.values[0]
                print("Need to expand before, currently have %s, need %s" % (current_length, flanking_length))
                bef_regions = compo_enlarger(before_subset, "forward","bef",compo_table, flanking_length, current_isolate, current_length)

            else:
                rec_start_bef = rec_end_bef - flanking_length
                bef_regions = [[rec_start_bef, rec_end_bef]]

            if aft_contig == 0:
                current_contig_bounds = bounds_of_contig(contig_file, contig_checker(contig_file, mge_hits))
                current_length = current_contig_bounds.values[1] - rec_start_aft
                print("Need to expand before, currently have %s, need %s" % (current_length, flanking_length))
                aft_regions = compo_enlarger(after_subset, "forward","aft", compo_table,
                                             flanking_length, current_isolate, current_length)
            else:
                rec_end_aft = rec_start_aft + flanking_length
                aft_regions = [[rec_start_aft, rec_end_aft]]


            reccy_starts_bef.append(rec_start_bef)
            reccy_ends_bef.append(rec_end_bef)
            reccy_starts_aft.append(rec_start_aft)
            reccy_ends_aft.append(rec_end_aft)

        elif mge_hits[0] > mge_hits[1]:
            rec_start_bef = before_subset.iloc[0, 6]
            rec_end_bef = before_subset.iloc[0, 7]
            bef_hits = [rec_start_bef, rec_start_bef + flanking_length]
            rec_start_aft = after_subset.iloc[0, 6]
            rec_end_aft = after_subset.iloc[0, 7]
            aft_hits = [rec_end_aft, rec_end_aft - flanking_length]
            bef_contig = contig_checker(contig_file, bef_hits)
            aft_contig = contig_checker(contig_file, aft_hits)


            if bef_contig == 0:
                current_contig_bounds = bounds_of_contig(contig_file, contig_checker(contig_file, mge_hits))
                current_length = current_contig_bounds.values[1] - rec_start_bef
                print("Need to expand reverse before, currently have %s, need %s" % (current_length, flanking_length))
                bef_regions = compo_enlarger(before_subset, "reverse", "bef", compo_table,
                                             flanking_length, current_isolate, current_length)


            else:
                rec_end_bef = rec_start_bef + flanking_length
                bef_regions = [[rec_start_bef, rec_end_bef]]

            if aft_contig == 0:
                current_contig_bounds = bounds_of_contig(contig_file, contig_checker(contig_file, mge_hits))
                current_length = rec_end_aft - current_contig_bounds.values[0]
                print("Need to expand reverse after, currently have %s, need %s" % (current_length, flanking_length))
                aft_regions = compo_enlarger(after_subset, "reverse", "aft", compo_table,
                                             flanking_length, current_isolate, current_length)
            else:
                rec_start_aft = rec_end_aft - flanking_length
                aft_regions = [[rec_start_aft, rec_end_aft]]

            rec_end_bef = rec_start_bef + flanking_length
            rec_start_aft = rec_end_aft - flanking_length

            reccy_starts_bef.append(rec_start_bef)
            reccy_ends_bef.append(rec_end_bef)
            reccy_starts_aft.append(rec_start_aft)
            reccy_ends_aft.append(rec_end_aft)
        regions_bef.append(bef_regions)
        regions_aft.append(aft_regions)







    out_df['bef_rec_start'] = pandas.Series(data=reccy_starts_bef, index=out_df.index)
    out_df['bef_rec_end'] = pandas.Series(data=reccy_ends_bef, index=out_df.index)
    out_df['aft_rec_start'] = pandas.Series(data=reccy_starts_aft, index=out_df.index)
    out_df['aft_rec_end'] = pandas.Series(data=reccy_ends_aft, index=out_df.index)

    out_df['control_id'] = pandas.Series(data=control_id, index=out_df.index)
    out_df['bef_rec_start_cont'] = pandas.Series(data=reccy_starts_bef_control, index=out_df.index)
    out_df['bef_rec_end_cont'] = pandas.Series(data=reccy_ends_bef_control, index=out_df.index)
    out_df['aft_rec_start_cont'] = pandas.Series(data=reccy_starts_aft_control, index=out_df.index)
    out_df['aft_rec_end_cont'] = pandas.Series(data=reccy_ends_aft_control, index=out_df.index)



    return out_df, regions_bef, regions_aft

def extracting_flanks(out_df, out_dir, ref_name, fasta_directory, regions_bef, regions_aft):


    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)




    for posdnuos in range(len(out_df.index)):
        current_row = out_df.iloc[posdnuos]
        current_id = current_row.iloc[0]
        current_start_bef = (current_row['bef_rec_start'] - 1)
        current_end_bef = current_row['bef_rec_end']
        current_start_aft = (current_row['aft_rec_start'] - 1)
        current_end_aft = current_row['aft_rec_end']

        current_bef_gappers = regions_bef[posdnuos]
        current_aft_gappers = regions_aft[posdnuos]


        current_control_id = current_row['control_id']
        current_start_bef_cont = int((current_row['bef_rec_start_cont'] - 1))
        current_end_bef_cont = int(current_row['bef_rec_end_cont'])
        current_start_aft_cont = int((current_row['aft_rec_start_cont'] - 1))
        current_end_aft_cont = int(current_row['aft_rec_end_cont'])


        if "ref!" in current_control_id:
            fasta_file_control = fasta_directory + re.sub("#", "_", ref_name) + ".dna"
            act_control_id = current_control_id



        else:
            act_control_id = re.split("!", current_control_id)[0]
            fasta_file_control = fasta_directory + act_control_id + ".dna"


        fasta_file = fasta_directory + current_id + ".dna"


        if os.path.isfile(fasta_file) == False:
            fasta_file = fasta_directory + current_id + ".contigs_velvet.dna"
            if os.path.isfile(fasta_file) == False:

                dna_id = re.sub("#", "_", current_id)
                fasta_file = fasta_directory + dna_id + ".contigs_velvet.dna"
                if os.path.isfile(fasta_file) == False:
                    fasta_file = fasta_directory + dna_id + ".dna"

                    if os.path.isfile(fasta_file) == False:
                        print(fasta_file)
                        sys.exit("No DNA file for this isolate:")

        if os.path.isfile(fasta_file_control) == False:
            fasta_file_control = fasta_directory + act_control_id + ".contigs_velvet.dna"
            if os.path.isfile(fasta_file_control) == False:

                dna_id = re.sub("#", "_", act_control_id)
                fasta_file_control = fasta_directory + dna_id + ".contigs_velvet.dna"
                if os.path.isfile(fasta_file_control) == False:
                    fasta_file_control = fasta_directory + dna_id + ".dna"

                    if os.path.isfile(fasta_file_control) == False:
                        print(fasta_file_control)
                        sys.exit("No DNA file for this control isolate:")


        if sum(out_df['isolate_id'].str.count(current_id)) > 1:
            current_id = current_id + "_NUM_" + str(posdnuos)
            current_control_id = current_control_id + "_NUM_" + str(posdnuos)


        new_blast_file = out_dir + "/" + current_id + "_whole_blast_seq.fasta"
        new_blast_bef = out_dir + "/" + current_id + "_before_flank.fasta"
        new_blast_aft = out_dir + "/" + current_id + "_after_flank.fasta"
        new_blast_file_control = out_dir + "/" + current_control_id + "_control_whole_blast_seq.fasta"
        new_blast_bef_control = out_dir + "/" + current_control_id + "_control_before_flank.fasta_control"
        new_blast_aft_control = out_dir + "/" + current_control_id + "_control_after_flank.fasta_control"



        with open(new_blast_file, "w") as blasty:
            for seq_record in SeqIO.parse(fasta_file, "fasta"):

                #blasty.write(str(seq_record.id) + "\n")
                #blasty.write(str(seq_record.seq[current_start:current_end]))

                new_seqqys_bef_list = []
                new_seqqys_aft_list = []
                for praha in  range(len(current_bef_gappers)):
                    current_elements = current_bef_gappers[praha]
                    current_seqqys_bef = str(seq_record.seq[(current_elements[0] - 1):current_elements[1]])
                    new_seqqys_bef_list.append(current_seqqys_bef)
                for prague in range(len(current_aft_gappers)):

                    current_elements = current_aft_gappers[prague]

                    current_seqqys_aft = str(seq_record.seq[(current_elements[0] - 1):current_elements[1]])

                    new_seqqys_aft_list.append(current_seqqys_aft)


                new_seqqys_bef = ''.join(new_seqqys_bef_list)
                new_seqqys_aft = ''.join(new_seqqys_aft_list)



                new_seqqys = new_seqqys_bef + new_seqqys_aft
                seqqy_id = str(seq_record.id)
                output_sequence = SeqIO.SeqRecord(Seq(new_seqqys), id=seqqy_id)
                output_before = SeqIO.SeqRecord(Seq(new_seqqys_bef), id=seqqy_id)
                output_after = SeqIO.SeqRecord(Seq(new_seqqys_aft), id=seqqy_id)

                SeqIO.write(output_sequence, new_blast_file, "fasta")
                SeqIO.write(output_before, new_blast_bef, "fasta")
                SeqIO.write(output_after, new_blast_aft, "fasta")

        with open(new_blast_file_control, "w") as blasty:
            for seq_record in SeqIO.parse(fasta_file_control, "fasta"):

                #blasty.write(str(seq_record.id) + "\n")
                #blasty.write(str(seq_record.seq[current_start:current_end]))




                new_seqqys_bef = str(seq_record.seq[current_start_bef_cont:current_end_bef_cont])
                new_seqqys_aft = str(seq_record.seq[current_start_aft_cont:current_end_aft_cont])
                new_seqqys = new_seqqys_bef + new_seqqys_aft
                seqqy_id = str(seq_record.id)
                output_sequence = SeqIO.SeqRecord(Seq(new_seqqys), id=seqqy_id)
                output_before = SeqIO.SeqRecord(Seq(new_seqqys_bef), id=seqqy_id)
                output_after = SeqIO.SeqRecord(Seq(new_seqqys_aft), id=seqqy_id)

                SeqIO.write(output_sequence, new_blast_file_control, "fasta")
                SeqIO.write(output_before, new_blast_bef_control, "fasta")
                SeqIO.write(output_after, new_blast_aft_control, "fasta")

###############################################################################
## So now we'll run through the detected gubbins transformation events to try #
## and get the isolate with the least amount of snp changes from the imported #
## sequence, this will allow us to search BLAST using this isolates sequence ##
## to ensure close similarity to the actual import ############################
###############################################################################

if __name__ == '__main__':
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Flanks Extractor ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    overall_start = time.perf_counter()

    input_args = get_options()
    reccy_hits = input_args.reccy_hits
    pyt_csv = input_args.hit_csv
    act_compos = input_args.act_compos
    flanking_length = input_args.flank_length
    fasta_directory = input_args.dna_dir
    out_dir = input_args.out_dir
    out_name = input_args.out_name
    contig_bounds = input_args.contig_bounds
    prop_hits = input_args.proper_hits
    tot_flanks_csv = pandas.DataFrame()

    pyt_csv = pandas.read_csv(pyt_csv)
    reccy_hits = pandas.read_csv(reccy_hits)
    proper_hits = pandas.read_csv(prop_hits)

    nice_ids_tot = nice_proper_hits_ids(proper_hits.iloc[:,0].tolist())

    base_loc = input_args.gubbins_res
    unique_clusters = reccy_hits['cluster_name'].unique()

    seq_clus = 1
    for cluster in unique_clusters:
        print("On cluster: %s, %s of %s" % (cluster, seq_clus, len(unique_clusters)))
        tic_cluster = time.perf_counter()
        current_dat = reccy_hits[reccy_hits['cluster_name'] == cluster]
        current_pyt = pyt_csv[pyt_csv['cluster_name'] == cluster]
        current_ref_name = current_pyt['ref_name'].iloc[0]
        current_dir = base_loc + cluster
        try:
            cluster_files = os.listdir(current_dir)
        except:
            current_dir = current_dir + "_run_data"
            cluster_files = os.listdir(current_dir)
        tree_indexio = [k for k, s in enumerate(cluster_files) if "node_labelled.final_tree.tre" in s]
        branch_mutations = [k for k, s in enumerate(cluster_files) if "_per_branch_mutations.csv" in s]
        embl_reccy = [k for k, s in enumerate(cluster_files) if "_recombinations.csv" in s]

        tree_loc = current_dir + "/" + cluster_files[tree_indexio[0]]
        embl_branch_loc = current_dir + "/" + cluster_files[branch_mutations[0]]
        embl_rec_loc = current_dir + "/" + cluster_files[embl_reccy[0]]

        tree = dendropy.Tree.get(path=tree_loc, schema="newick", preserve_underscores=True)
        branch_mutations = pandas.read_csv(embl_branch_loc)
        embl_reccy_csv = pandas.read_csv(embl_rec_loc)

        if cluster != "gpsc.136":
            continue

        print(current_dat)
        flanks_csv, regions_bef, regions_aft = isolate_narrow(current_dat, current_pyt, tree, embl_reccy_csv,
                                                              branch_mutations, current_ref_name, flanking_length,
                                                              contig_bounds, nice_ids_tot)
        if isinstance(flanks_csv, str):
            continue

        extract_flanks = extracting_flanks(flanks_csv, out_dir,current_ref_name, fasta_directory, regions_bef, regions_aft)

        tot_flanks_csv = tot_flanks_csv.append(flanks_csv, ignore_index=True, sort=False)
        seq_clus += 1
        print("Extracted %s for cluster: %s" % (len(flanks_csv.index), cluster))

    toc_tot = time.perf_counter()
    tot_flanks_csv.to_csv(path_or_buf=out_name, index=False)

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Flanks extractor took: %s (seconds)" % (toc_tot - overall_start))
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Flanks extracted ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

