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
    parser.add_argument('--gff_csv',required=True, help="reference isolate gff csv", type=str)
    parser.add_argument('--contig_bounds', required=True, help='bounds of contigs used in reconstruction', type = str)
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
    if len(starting_seq_indy) < 1:
        print("No FASTA in GFF for this isolate: %s" % isolate_id)
        return "No hit"

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


                narrowed_gff.iloc[contig_rows_indy,3] = narrowed_gff.iloc[contig_rows_indy, 3] + num_to_add - 1
                narrowed_gff.iloc[contig_rows_indy, 4] =  narrowed_gff.iloc[contig_rows_indy, 4] + num_to_add - 1

                last_indy = max(contig_rows_indy)

        gff_index_to_drop = range((last_indy + 1), (len(narrowed_gff.index) ))
        out_gff_csv = narrowed_gff.drop(gff_index_to_drop)

    ## test out finder ##


    out_gff_csv = out_gff_csv.reset_index(drop=True)

    return out_gff_csv

def outside_control(insertion_node, tree, example_id, act_comp_dir, ref_insertion, pyt_csv,
                    current_res, prev_res, current_gene, gff_csv, contig_loc):
    ## This function identifies control isolates outside of insertion node and the regions
    ## to look at
    res_pyt = pyt_csv[pyt_csv['profile'] == current_res]
    nice_ids_tot = nice_proper_hits_ids(res_pyt.iloc[:, 0].tolist())
    prev_res_pyt = pyt_csv[pyt_csv['profile'] == "S"]
    nice_ids_prev = nice_proper_hits_ids(prev_res_pyt.iloc[:, 0].tolist())

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

    outside_iso = 0
    isolate_before_end = False
    isolate_after_start = False
    for ref in length_dict_sorted:
        if ref in nice_ids_tot and ref not in nice_ids_prev:
            continue

        outside_iso = ref
        contig_suffix = "#contig_bounds.csv"
        contig_isolate = re.sub("#", "_", outside_iso)
        contig_file_path = contig_loc + contig_isolate + contig_suffix
        iso_contig = pandas.read_csv(contig_file_path)
        current_gff_loc, ref_loc, cluster_name = gff_finder(gff_csv, outside_iso, True)

        current_gff = pandas.read_csv(current_gff_loc.iloc[0], sep='\t',
                                      names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                             'attributes'],
                                      header=None)
        k = 1
        if len(iso_contig.index) > 1:
            current_gff = gff_to_dna(current_gff, iso_contig, outside_iso, input_k=k)
            if isinstance(current_gff, str):
                continue

        if current_gene == "pbp1A":
            correct_length, pbp_row = search_for_gene(current_gff, ['pbp1A', 'ponA'], 2159,
                                                      100, False, None)

        elif current_gene == "pbp2B":  ## pbp 2b
            correct_length, pbp_row = search_for_gene(current_gff, ['pbp2B', 'penA'], 2042,
                                                      100, False, None)
        ## pbp 2x
        elif current_gene == "pbp2X":
            correct_length, pbp_row = search_for_gene(current_gff, ['pbp2X', 'pbpX'], 2252,
                                                      100, False, None)


        if pbp_row.iloc[0,6] == "+":
            isolate_before_end = pbp_row.iloc[0,3]
            isolate_after_start = pbp_row.iloc[0,4]
            strand = "forward"
            break
        elif pbp_row.iloc[0,6] == "-":
            isolate_before_end = pbp_row.iloc[0, 3]
            isolate_after_start = pbp_row.iloc[0, 4]
            strand = "reverse"
            break

    return outside_iso, isolate_before_end, isolate_after_start, strand

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
    ## Function to extract regions of the genome where the pbp genes are resident.
    ## Need to find the ACT copmarison to the reference that matches the pbp gene, then
    ## find the appropriate length regions to add on.
    ## Input: start_act: The pbp gene location.
    ##              ori: The orientation of the hit
    ##              pos: The position of the hit, either before or after the gene
    ##        act_compo: The act_comparison csv
    ##    target_length: The flanking region length

    ## First lets find the ACT compo for the pbp gene in question.

    if ori == "forward":
        pbp_act = act_compo[(act_compo['qstart'] <= start_act[0]) &
                            (act_compo['qend'] >= start_act[1])]
        if pbp_act.empty:
            if pos == "bef":
                pbp_act = act_compo[(act_compo['qstart'] <= start_act[0]) &
                                    (act_compo['qend'] >= start_act[0])]
                if pbp_act.empty:
                    print("Warning not using immediate flanks to gene")
                    pbp_act = act_compo[(act_compo['qstart'] <= start_act[0]) &
                                        (act_compo['qend'] <= start_act[1])]
                    pbp_act = pbp_act.sort_values(by='qend', ascending=False)
            elif pos == "aft":
                pbp_act = act_compo[(act_compo['qstart'] <= start_act[1]) &
                                    (act_compo['qend'] >= start_act[1])]
                if pbp_act.empty:
                    print("Warning not using immediate flanks to gene")
                    pbp_act = act_compo[(act_compo['qstart'] >= start_act[0]) &
                                        (act_compo['qend'] >= start_act[1])]
                    pbp_act = pbp_act.sort_values(by='qstart', ascending=True)
        if pos == "bef":
            current_length = start_act[0] - pbp_act['qstart'].iloc[0]
            starting_act = pbp_act.iloc[0]
            starting_gap = [starting_act['qstart'], start_act[0]]
        elif pos == "aft":
            current_length = pbp_act['qend'].iloc[0] - start_act[1]
            starting_act = pbp_act.iloc[0]
            starting_gap = [start_act[1], starting_act['qend']]

    elif ori == "reverse":
        pbp_act = act_compo[(act_compo['qstart'] <= start_act[1]) &
                            (act_compo['qend'] >= start_act[0])]
        if pbp_act.empty:
            if pos == "bef":
                pbp_act = act_compo[(act_compo['qstart'] <= start_act[0]) &
                                    (act_compo['qend'] >= start_act[0])]
                if pbp_act.empty:
                    print("Warning not using immediate flanks to gene")
                    pbp_act = act_compo[(act_compo['qstart'] >= start_act[1]) &
                                        (act_compo['qend'] >= start_act[0])]
                    pbp_act = pbp_act.sort_values(by='qstart', ascending=True)
            elif pos == "aft":
                pbp_act = act_compo[(act_compo['qstart'] <= start_act[1]) &
                                    (act_compo['qend'] >= start_act[1])]
                if pbp_act.empty:
                    print("Warning not using immediate flanks to gene")
                    pbp_act = act_compo[(act_compo['qstart'] <= start_act[0]) &
                                        (act_compo['qend'] <= start_act[1])]
                    pbp_act = pbp_act.sort_values(by='qend', ascending=False)
        if pos == "bef":
            current_length = pbp_act['qend'].iloc[0] - start_act[0]
            starting_act = pbp_act.iloc[0]
            starting_gap = [start_act[0], starting_act['qend']]
        elif pos == "aft":
            current_length = start_act[1] - pbp_act['qstart'].iloc[0]
            starting_act = pbp_act.iloc[0]
            starting_gap = [starting_act['qstart'], start_act[1]]



    ## Define the search direction for the next hits
    if (pos == "aft" and ori == "forward") | (pos == "bef" and ori == "reverse"):
        direction = "downstream"
    else:
        direction = "upstream"

    if starting_act['sstart'] < starting_act['send']:
        act_ori = "forward"
    else:
        act_ori = "reverse"
    #if direction == "downstream":

    print(starting_act)



    gap_list = [starting_gap]

    while current_length < target_length:
        if direction == "downstream":
            if act_ori == "forward":
                ## This catches any reverse hits
                next_act_event_send = act_compo[act_compo['send'] > starting_act['send']]
                next_act_event_send = next_act_event_send.sort_values(by = ['send'], ascending=True)
                next_act_event_send_row = next_act_event_send.iloc[0]

                ref_start_send = next_act_event_send_row.iloc[9]


                ## This for forward hits
                next_act_event_sstart = act_compo[act_compo['sstart'] > starting_act['send']]
                next_act_event_sstart = next_act_event_sstart.sort_values(by = ['sstart'], ascending=True)
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


                if ref_start_send <= ref_start_sstart:
                    distance_hits = abs(next_act_event_send_row.iloc[7] - next_act_event_send_row.iloc[6])
                    current_length = current_length + distance_hits
                    if current_length >= target_length:
                        gappers = [next_act_event_send_row.iloc[6] , next_act_event_send_row.iloc[7] - (current_length - target_length)]
                        gap_list.append(gappers)
                    else:
                        gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                        gap_list.append(gappers)
                        starting_act = next_act_event_send_row
                elif ref_start_send > ref_start_sstart:
                    distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                    current_length = current_length + distance_hits
                    if current_length >= target_length:
                        gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7] - (current_length - target_length)]
                        gap_list.append(gappers)
                    else:
                        gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                        gap_list.append(gappers)
                        starting_act = next_act_event_sstart_row
            elif act_ori == "reverse":
                ## This catches any reverse hits
                next_act_event_send = act_compo[act_compo['send'] < starting_act['send']]
                next_act_event_send = next_act_event_send.sort_values(by=['send'], ascending=False)
                next_act_event_send_row = next_act_event_send.iloc[0]

                ref_start_send = next_act_event_send_row.iloc[9]

                ## This for forward hits
                next_act_event_sstart = act_compo[act_compo['sstart'] < starting_act['send']]
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
                    current_length = current_length + distance_hits
                    if current_length >= target_length:
                        gappers = [next_act_event_send_row.iloc[6] , next_act_event_send_row.iloc[7] - (current_length - target_length)]
                        gap_list.append(gappers)
                    else:
                        gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                        gap_list.append(gappers)
                        starting_act = next_act_event_send_row
                elif ref_start_send < ref_start_sstart:
                    distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                    current_length = current_length + distance_hits
                    if current_length >= target_length:
                        gappers = [next_act_event_send_row.iloc[6] , next_act_event_send_row.iloc[7] - (current_length - target_length)]
                        gap_list.append(gappers)
                    else:
                        gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                        gap_list.append(gappers)
                        starting_act = next_act_event_sstart_row


        elif direction == "upstream":

            if act_ori == "forward":
                print(current_length)
                ## This catches any forward hits
                next_act_event_send = act_compo[act_compo['send'] < starting_act['sstart']]
                next_act_event_send = next_act_event_send.sort_values(by=['send'], ascending=False)
                next_act_event_send_row = next_act_event_send.iloc[0]

                ref_start_send = next_act_event_send_row.iloc[9]

                ## This for reverse hits
                next_act_event_sstart = act_compo[act_compo['sstart'] < starting_act['sstart']]
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

                print(next_act_event_sstart_row )
                print(next_act_event_send_row)
                if ref_start_send >= ref_start_sstart:
                    distance_hits = abs(next_act_event_send_row.iloc[7] - next_act_event_send_row.iloc[6])
                    current_length = current_length + distance_hits
                    if current_length >= target_length:
                        gappers = [next_act_event_send_row.iloc[6] +  (current_length - target_length) , next_act_event_send_row.iloc[7]]
                        gap_list.append(gappers)
                    else:
                        gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                        gap_list.append(gappers)
                        starting_act = next_act_event_send_row
                elif ref_start_send < ref_start_sstart:
                    distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                    current_length = current_length + distance_hits
                    if current_length >= target_length:
                        gappers = [next_act_event_send_row.iloc[6] +  (current_length - target_length) , next_act_event_send_row.iloc[7]]
                        gap_list.append(gappers)
                    else:
                        gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                        gap_list.append(gappers)
                        starting_act = next_act_event_sstart_row
            elif act_ori == "reverse":
                ## This catches any reverse hits
                next_act_event_send = act_compo[act_compo['send'] > starting_act['sstart']]
                next_act_event_send = next_act_event_send.sort_values(by=['send'], ascending=True)
                next_act_event_send_row = next_act_event_send.iloc[0]

                ref_start_send = next_act_event_send_row.iloc[9]

                ## This for forward hits
                next_act_event_sstart = act_compo[act_compo['sstart'] > starting_act['sstart']]
                next_act_event_sstart = next_act_event_sstart.sort_values(by=['sstart'], ascending=True)
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

                if ref_start_send <= ref_start_sstart:
                    distance_hits = abs(next_act_event_send_row.iloc[7] - next_act_event_send_row.iloc[6])
                    current_length = current_length + distance_hits
                    if current_length >= target_length:
                        gappers = [next_act_event_send_row.iloc[6] +  (current_length - target_length) , next_act_event_send_row.iloc[7]]
                        gap_list.append(gappers)
                    else:
                        gappers = [next_act_event_send_row.iloc[6], next_act_event_send_row.iloc[7]]
                        gap_list.append(gappers)
                        starting_act = next_act_event_send_row
                elif ref_start_send > ref_start_sstart:
                    distance_hits = abs(next_act_event_sstart_row.iloc[7] - next_act_event_sstart_row.iloc[6])
                    current_length = current_length + distance_hits
                    if current_length >= target_length:
                        gappers = [next_act_event_send_row.iloc[6] +  (current_length - target_length) , next_act_event_send_row.iloc[7]]
                        gap_list.append(gappers)
                    else:
                        gappers = [next_act_event_sstart_row.iloc[6], next_act_event_sstart_row.iloc[7]]
                        gap_list.append(gappers)
                        starting_act = next_act_event_sstart_row


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

def gff_finder(gff_csv, isolate_id, clus_name):
    ## Function to get the location of an isolates gff file
    isolate_check = isolate_id + "\."
    isolate_rows = gff_csv['isolate'].str.contains(isolate_check)
    isolate_row_indy = isolate_rows.where(isolate_rows == True)
    isolate_row_indy = isolate_row_indy.index[isolate_row_indy == True].tolist()
    if len(isolate_row_indy) != 1:
        print(isolate_id)
        print(isolate_row_indy)
    #print(gff_csv.head(), isolate_row_indy)
    isolate_loc = gff_csv.iloc[isolate_row_indy,0]
    isolate_ref = gff_csv.iloc[isolate_row_indy,1]
    if clus_name:
        cluster_name = gff_csv.iloc[isolate_row_indy,2]
    else:
        cluster_name = "BOO!"

    return isolate_loc, isolate_ref, cluster_name

def search_for_gene(ref_in,name,gene_length,tol,correct_length,gene_rower):

    for gene in name:
        if not correct_length:
            gene_row = ref_in['attributes'].str.contains(gene)
            gene_row_indy = gene_row.where(gene_row == True)
            gene_row_indy = gene_row_indy.index[gene_row_indy == True].tolist()
            gene_rower = ref_in.iloc[gene_row_indy]
            if gene_rower.empty == False:
                gene_len = [abs(int(gene_rower.iloc[0,4]) - int(gene_rower.iloc[0,3]))]

                overhang = [gene_length - tol, gene_length + tol]

                if overhang[0] <= gene_len[0] <= overhang[1]:
                    correct_length = True
                else:
                    sys.stderr.write('Found gene' + gene + ' but wrong length: ' + str(gene_len[0]) + ', expected: ' + str(gene_length) + '\n')

    return correct_length,gene_rower

def isolate_narrow(reccy_hits, pyt_csv, tree, reccy_csv_gubbins, mut_bases_csv, reference_id,
                   flanking_length, contig_loc, gff_csv):
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
    insert_start = []
    insert_end = []
    previous_resistance = []
    resistance = []
    cluster = []
    node_tips = []

    reccy_starts_bef_control = []
    reccy_ends_bef_control = []
    reccy_starts_aft_control = []
    reccy_ends_aft_control = []
    pbp_insert_start_cont = []
    pbp_insert_end_cont = []

    regions_bef = []
    regions_aft = []
    gene_name = []

    out_df = pandas.DataFrame()

    ## This loops through the reccy hits, finding the isolate in the clade with the fewest SNPs present around the
    ## pbp insertion point, appending them and their insert points in the reference to the out df. Also finds the
    ## control isolates position.

    for k in range(len(reccy_hits.index)):
        current_row = reccy_hits.iloc[k]
        isolate_row = current_row['isolate_example']
        insertion_node = current_row['insertion_node']
        start_node_to_reccy = current_row['start_nodes']
        end_node_to_reccy = current_row['end_nodes']
        ref_insert_start = current_row['start_insert']
        ref_insert_end = current_row['end_insert']
        current_resistance = current_row['resistance']
        current_prev_res = current_row['previous_resistance']
        lower_bounds = ref_insert_start - flanking_length
        upper_bounds = ref_insert_end + flanking_length
        current_gene = current_row['gene']
        num_tips = current_row['insertion_node_tip_nums']
        current_cluster = current_row['cluster_name']


        ref_insert_list = [ref_insert_start, ref_insert_end]
        nodes_from_insertion = leaf_tips(tree, example_id=isolate_row, current_node=insertion_node)

        res_isolates = []


        ## So lets run through all the leaf tips from this inseretion node, and narrow down to only
        ## the ones with the gain of resistance for this insertion.

        for l in range(len(nodes_from_insertion)):
            august_isolate = nodes_from_insertion[l]
            if pyt_csv['isolate'].isin([august_isolate]).any():
                pyt_row = pyt_csv[pyt_csv['isolate'] == august_isolate]
                if pyt_row['profile'].isin([current_resistance]).any():
                    res_isolates.append(august_isolate)



        act_bounders = [lower_bounds, upper_bounds]

        total_life_forever = []

        ## Goes through all the leaf tips for this inseriton and then selects the isolate with the
        ## fewest snps around the insertion loci
        ##

        for m in range(len(res_isolates)):
            current_isolate = res_isolates[m]
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
                                                                    reccy_detected['end_rec'] <= act_bounders[1])) |
                                                    ((reccy_detected['start_rec'] <= act_bounders[0]) &
                                                     (reccy_detected['end_rec'] >= act_bounders[1]))]
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
                                                        reccy_detected['end_rec'] <= act_bounders[1]))|
                                                    ((reccy_detected['start_rec'] <= act_bounders[0]) &
                                                     (reccy_detected['end_rec'] >= act_bounders[1]))]
                if not reccy_hits_loc.empty:
                    reccy_muts += reccy_hits_loc['snp_number'].sum()

                mut_changes_sub = mut_bases_csv[mut_bases_csv['end_node'] == current_end]
                mut_changes = mut_changes_sub.iloc[0, 3:26].sum()
                branch_muts += mut_changes
                total_life_forever.append(reccy_muts + branch_muts)



        res_seq_to_look_at = total_life_forever.index(min(total_life_forever))
        snp_count_indiv = min(total_life_forever)
        res_id = res_isolates[res_seq_to_look_at]
        orig_mge_id = res_id
        orig_snp_count = snp_count_indiv
        while res_id in isolate_id:
            total_life_forever.remove(snp_count_indiv)
            res_isolates.remove(res_id)
            if len(res_isolates) == 0:
                res_id = orig_mge_id
                snp_count_indiv = orig_snp_count
                break
            mge_seq_to_look_at = total_life_forever.index(min(total_life_forever))
            snp_count_indiv = min(total_life_forever)
            res_id = res_isolates[mge_seq_to_look_at]



        res_deets = pyt_csv[pyt_csv['isolate'] == res_id]
        contig_suffix = "#contig_bounds.csv"
        contig_isolate = re.sub("#", "_", res_id)
        contig_file_path = contig_loc + contig_isolate + contig_suffix
        iso_contig = pandas.read_csv(contig_file_path)

        current_gff_loc, ref_loc, cluster_name = gff_finder(gff_csv, res_id, True)

        current_gff = pandas.read_csv(current_gff_loc.iloc[0], sep='\t',
                                      names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                             'attributes'],
                                      header=None)
        input_k = 1
        if len(iso_contig.index) > 1:
            current_gff = gff_to_dna(current_gff, iso_contig, res_id, input_k=input_k)
            if isinstance(current_gff, str):
                continue

        if current_gene == "pbp1A":
            correct_length, pbp_row = search_for_gene(current_gff, ['pbp1A', 'ponA'], 2159,
                                                           100, False, None)

        elif current_gene == "pbp2B":    ## pbp 2b
            correct_length, pbp_row = search_for_gene(current_gff, ['pbp2B', 'penA'], 2042,
                                                           100, False, None)
        ## pbp 2x
        elif current_gene == "pbp2X":
            correct_length, pbp_row = search_for_gene(current_gff, ['pbp2X', 'pbpX'], 2252,
                                                           100, False, None)

        ## Get the locs of the pbps
        if pbp_row.iloc[0,6] == "+":
            pbp_start = pbp_row.iloc[0, 3]
            pbp_end = pbp_row.iloc[0,4]
            pbp_bef_start = pbp_start - flanking_length
            pbp_bef_end = pbp_start
            pbp_aft_start = pbp_end
            pbp_aft_end = pbp_end + flanking_length
        elif pbp_row.iloc[0,6] == "-":
            pbp_start = pbp_row.iloc[0, 4]
            pbp_end = pbp_row.iloc[0, 3]
            pbp_bef_start = pbp_start
            pbp_bef_end = pbp_start + flanking_length
            pbp_aft_start = pbp_end - flanking_length
            pbp_aft_end = pbp_end
        ## Check if reference among the tips for this node insertion.

        cont_id, cont_start, cont_end, cont_strand = outside_control(insertion_node, tree, isolate_row, act_compos,
                                                        ref_insert_list, pyt_csv, current_resistance, current_prev_res,
                                                        current_gene, gff_csv, contig_loc)
        if isinstance(cont_id, int):
            print("Can't extract control for this isolate, skipping: %s" % mge_id)
            continue

        csv_ref_name = cont_id + "!" + res_id

        control_id.append(csv_ref_name)



        if cont_strand == "forward":
            rec_end_bef = cont_start
            rec_start_aft = cont_end
            rec_start_bef = rec_end_bef - flanking_length
            rec_end_aft = rec_start_aft + flanking_length

        elif cont_strand == "reverse":
            rec_end_bef = cont_end + flanking_length
            rec_start_aft = cont_start - flanking_length
            rec_start_bef = cont_end
            rec_end_aft = cont_start



        reccy_starts_bef_control.append(rec_start_bef)
        reccy_ends_bef_control.append(rec_end_bef)
        reccy_starts_aft_control.append(rec_start_aft)
        reccy_ends_aft_control.append(rec_end_aft)
        pbp_insert_start_cont.append(cont_start)
        pbp_insert_end_cont.append(cont_end)


        isolate_id.append(res_id)
        mge_bef.append(pbp_start)
        mge_aft.append(pbp_end)
        before_start_region.append(pbp_bef_start)
        before_end_region.append(pbp_bef_end)
        after_start_region.append(pbp_aft_start)
        after_end_region.append(pbp_aft_end)
        snp_count.append(snp_count_indiv)
        insert_start.append(ref_insert_start)
        insert_end.append(ref_insert_end)
        gene_name.append(current_gene)
        previous_resistance.append(current_prev_res)
        resistance.append(current_resistance)
        cluster.append(current_cluster)
        node_tips.append(num_tips)

        if res_id == "11511_7#11":
            print("@~@~@~@~~@~@~@~@~@~@~@~@~@~@~@~@~@~~@~@~@~@~@~@~~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@")
            print("PBP before: ", pbp_start)
            print("PBP after: ", pbp_end)
            print("Before start:", pbp_bef_start)
            print("Before end:", pbp_bef_end)
            print("After start:", pbp_aft_start)
            print("After end:", pbp_aft_end)
            print("@~@~@~@~~@~@~@~@~@~@~@~@~@~@~@~@~@~~@~@~@~@~@~@~~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@")

    if len(isolate_id) == 0:
        return "no", "no", "no"




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
    out_df['gene'] = pandas.Series(gene_name, index=out_df.index)
    out_df['previous_resistance'] = pandas.Series(previous_resistance, index=out_df.index)
    out_df['resistance'] = pandas.Series(resistance, index = out_df.index)
    out_df['Num_tips'] = pandas.Series(node_tips, index=out_df.index)
    out_df['cluster'] = pandas.Series(cluster, index=out_df.index)
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



        bef_contig = contig_checker(contig_file, bef_hits)
        aft_contig = contig_checker(contig_file, aft_hits)
        if mge_hits[1] > mge_hits[0]:
            pbp_ori = "forward"
        elif mge_hits[0] > mge_hits[1]:
            pbp_ori = "reverse"



        if bef_contig == 0:
            current_contig_bounds = bounds_of_contig(contig_file,contig_checker(contig_file, mge_hits))
            print(rec_end_bef)
            print(current_contig_bounds)
            current_length = rec_end_bef - current_contig_bounds.values[0]
            print("Need to expand before, currently have %s, need %s" % (current_length, flanking_length))
            bef_regions = [bef_hits.tolist()]
            #bef_regions = compo_enlarger(mge_hits, pbp_ori,"bef",compo_table, flanking_length, current_isolate, current_length)

        else:

            bef_regions = [bef_hits.tolist()]

        if aft_contig == 0:
            current_contig_bounds = bounds_of_contig(contig_file, contig_checker(contig_file, mge_hits))
            current_length = min((current_contig_bounds[1] - mge_hits[1]), (mge_hits[1] - current_contig_bounds[0]))
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(current_isolate)
            print(current_contig_bounds)
            print(mge_hits)
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("Need to expand after, currently have %s, need %s" % (current_length, flanking_length))
            aft_regions = [aft_hits.tolist()]
            #aft_regions = compo_enlarger(mge_hits, pbp_ori,"aft", compo_table, flanking_length, current_isolate, current_length)
            print(aft_regions)
        else:
            aft_regions = [aft_hits.tolist()]

        if current_isolate == "11511_7#11":
            print("@~@~@~@~~@~@~@~@~@~@~@~@~@~@~@~@~@~~@~@~@~@~@~@~~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@")
            print("PBP gappers: ", aft_regions)
            print("@~@~@~@~~@~@~@~@~@~@~@~@~@~@~@~@~@~~@~@~@~@~@~@~~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@")
            #sys.exit(1)

        reccy_starts_bef.append(rec_start_bef)
        reccy_ends_bef.append(rec_end_bef)
        reccy_starts_aft.append(rec_start_aft)
        reccy_ends_aft.append(rec_end_aft)


        regions_bef.append(bef_regions)
        regions_aft.append(aft_regions)








    out_df['bef_rec_start'] = pandas.Series(data=bef_hits.iloc[0], index=out_df.index)
    out_df['bef_rec_end'] = pandas.Series(data=bef_hits.iloc[1], index=out_df.index)
    out_df['aft_rec_start'] = pandas.Series(data=aft_hits.iloc[0], index=out_df.index)
    out_df['aft_rec_end'] = pandas.Series(data=aft_hits.iloc[1], index=out_df.index)

    out_df['control_id'] = pandas.Series(data=control_id, index=out_df.index)
    out_df['bef_rec_start_cont'] = pandas.Series(data=reccy_starts_bef_control, index=out_df.index)
    out_df['bef_rec_end_cont'] = pandas.Series(data=reccy_ends_bef_control, index=out_df.index)
    out_df['aft_rec_start_cont'] = pandas.Series(data=reccy_starts_aft_control, index=out_df.index)
    out_df['aft_rec_end_cont'] = pandas.Series(data=reccy_ends_aft_control, index=out_df.index)
    out_df['cont_pbp_start'] = pandas.Series(data=pbp_insert_start_cont, index=out_df.index)
    out_df['cont_pbp_end'] = pandas.Series(data=pbp_insert_end_cont, index=out_df.index)



    return out_df, regions_bef, regions_aft

def extracting_flanks(out_df, out_dir, ref_name, fasta_directory, regions_bef, regions_aft):


    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)




    for posdnuos in range(len(out_df.index)):
        current_row = out_df.iloc[posdnuos]
        current_id = current_row.iloc[0]
        current_gene_name = current_row['gene']
        pbp_start = current_row['mge_start']
        pbp_end = current_row['mge_end']


        if pbp_start < pbp_end:
            pbp_regions = [pbp_start, pbp_end]
        else:
            pbp_regions = [pbp_end, pbp_start]

        current_bef_gappers = regions_bef[posdnuos]
        current_aft_gappers = regions_aft[posdnuos]


        current_control_id = current_row['control_id']
        current_start_bef_cont = int((current_row['bef_rec_start_cont'] - 1))
        current_end_bef_cont = int(current_row['bef_rec_end_cont'])
        current_start_aft_cont = int((current_row['aft_rec_start_cont'] - 1))
        current_end_aft_cont = int(current_row['aft_rec_end_cont'])
        cont_pbp_start = int(current_row['cont_pbp_start'])
        cont_pbp_end = int(current_row['cont_pbp_end'])

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


        new_blast_file = out_dir + "/" + current_id + "_" + current_gene_name + "_whole_blast_seq.fasta"
        new_blast_bef = out_dir + "/" + current_id + "_" + current_gene_name + "_before_flank.fasta"
        new_blast_aft = out_dir + "/" + current_id + "_" + current_gene_name + "_after_flank.fasta"
        new_blast_pbp = out_dir + "/" + current_id + "_" + current_gene_name + "_pbp_sequence.fasta"
        new_blast_file_control = out_dir + "/" + current_control_id + "_" + current_gene_name + "_control_whole_blast_seq.fasta"
        new_blast_bef_control = out_dir + "/" + current_control_id + "_" + current_gene_name + "_control_before_flank.fasta_control"
        new_blast_aft_control = out_dir + "/" + current_control_id + "_" + current_gene_name + "_control_after_flank.fasta_control"
        new_blast_pbp_cont = out_dir + "/" + current_control_id + "_" + current_gene_name + "_control_pbp_sequence.fasta"

        if current_control_id == "10050_2#15!6569_4#17":
            print("~~~~~~~~~~~~~~~~~~~~~ after control locs ~~~~~~~~~~~~~~~~~~~~~~")
            print(current_start_aft_cont, current_end_aft_cont)

        with open(new_blast_file, "w") as blasty:
            for seq_record in SeqIO.parse(fasta_file, "fasta"):

                #blasty.write(str(seq_record.id) + "\n")
                #blasty.write(str(seq_record.seq[current_start:current_end]))

                new_seqqys_bef_list = []
                new_seqqys_aft_list = []
                for praha in  range(len(current_bef_gappers)):
                    current_elements = current_bef_gappers[praha]
                    current_seqqys_bef = str(seq_record.seq[(int(current_elements[0] - 1)):int(current_elements[1])])
                    new_seqqys_bef_list.append(current_seqqys_bef)
                for prague in range(len(current_aft_gappers)):

                    current_elements = current_aft_gappers[prague]

                    current_seqqys_aft = str(seq_record.seq[(int(current_elements[0] - 1)):int(current_elements[1])])

                    new_seqqys_aft_list.append(current_seqqys_aft)

                pbp_seqqys = str(seq_record.seq[int(pbp_regions[0] - 1): int(pbp_regions[1])])


                new_seqqys_bef = ''.join(new_seqqys_bef_list)
                new_seqqys_aft = ''.join(new_seqqys_aft_list)



                new_seqqys = new_seqqys_bef + new_seqqys_aft
                seqqy_id = str(seq_record.id)
                output_sequence = SeqIO.SeqRecord(Seq(new_seqqys), id=seqqy_id)
                output_before = SeqIO.SeqRecord(Seq(new_seqqys_bef), id=seqqy_id)
                output_after = SeqIO.SeqRecord(Seq(new_seqqys_aft), id=seqqy_id)
                pbp_seq = SeqIO.SeqRecord(Seq(pbp_seqqys), id = seqqy_id)

                SeqIO.write(output_sequence, new_blast_file, "fasta")
                SeqIO.write(output_before, new_blast_bef, "fasta")
                SeqIO.write(output_after, new_blast_aft, "fasta")
                SeqIO.write(pbp_seq, new_blast_pbp, "fasta")

        with open(new_blast_file_control, "w") as blasty:
            for seq_record in SeqIO.parse(fasta_file_control, "fasta"):

                #blasty.write(str(seq_record.id) + "\n")
                #blasty.write(str(seq_record.seq[current_start:current_end]))




                new_seqqys_bef = str(seq_record.seq[int(current_start_bef_cont):int(current_end_bef_cont)])
                new_seqqys_aft = str(seq_record.seq[int(current_start_aft_cont):int(current_end_aft_cont)])
                new_pbp_seqqys = str(seq_record.seq[int(cont_pbp_start - 1 ): int(cont_pbp_end)])
                new_seqqys = new_seqqys_bef + new_seqqys_aft
                seqqy_id = str(seq_record.id)
                output_sequence = SeqIO.SeqRecord(Seq(new_seqqys), id=seqqy_id)
                output_before = SeqIO.SeqRecord(Seq(new_seqqys_bef), id=seqqy_id)
                output_after = SeqIO.SeqRecord(Seq(new_seqqys_aft), id=seqqy_id)
                cont_pbp_output = SeqIO.SeqRecord(Seq(new_pbp_seqqys), id = seqqy_id)

                SeqIO.write(output_sequence, new_blast_file_control, "fasta")
                SeqIO.write(output_before, new_blast_bef_control, "fasta")
                SeqIO.write(output_after, new_blast_aft_control, "fasta")
                SeqIO.write(cont_pbp_output, new_blast_pbp_cont, "fasta")

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
    gff_csv = input_args.gff_csv

    tot_flanks_csv = pandas.DataFrame()

    pyt_csv = pandas.read_csv(pyt_csv)
    reccy_hits = pandas.read_csv(reccy_hits)
    gff_csv = pandas.read_csv(gff_csv)
    gff_csv.columns = ['isolate', 'reference', 'cluster']
    base_loc = input_args.gubbins_res
    unique_clusters = reccy_hits['cluster_name'].unique()

    seq_clus = 1
    for cluster in unique_clusters:
        print("On cluster: %s, %s of %s" % (cluster, seq_clus, len(unique_clusters)))
        tic_cluster = time.perf_counter()
        current_dat = reccy_hits[reccy_hits['cluster_name'] == cluster]
        current_pyt = pyt_csv[pyt_csv['cluster'] == cluster]
        current_gff_rows = gff_csv
        current_ref_loc = current_gff_rows['reference'].iloc[0]
        ## get the current ref name for contig bounds
        ref_name = os.path.basename(current_ref_loc)
        ref_name = re.sub("\..*[a-zA-Z]*$", "", ref_name)
        current_ref_name = ref_name

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
        print("Narrowing isolates:")
        narrow_start = time.perf_counter()

        flanks_csv, regions_bef, regions_aft = isolate_narrow(current_dat, current_pyt, tree, embl_reccy_csv,
                                                              branch_mutations, current_ref_name, flanking_length,
                                                              contig_bounds, current_gff_rows)
        if isinstance(flanks_csv, str):
            continue

        narrow_end = time.perf_counter()
        print("Took: %s" % (narrow_end - narrow_start))

        extract_flanks = extracting_flanks(flanks_csv, out_dir,current_ref_name, fasta_directory, regions_bef, regions_aft)

        tot_flanks_csv = tot_flanks_csv.append(flanks_csv, ignore_index=True, sort=False)
        seq_clus += 1
        print("Extracted %s for cluster: %s" % (len(flanks_csv.index), cluster))

    toc_tot = time.perf_counter()
    tot_flanks_csv.to_csv(path_or_buf=out_name, index=False)


    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Flanks extractor took: %s (seconds)" % (toc_tot - overall_start))
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Flanks extracted ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

