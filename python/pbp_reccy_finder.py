import pandas
from itertools import product
import numpy
import sys
import re
import difflib
import dendropy
import argparse
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
import time
from pastml.acr import pastml_pipeline

###############################################################################
## Load up the files ##########################################################
###############################################################################
def get_options():

    purpose = '''This is a python script to output details about specific recombinations that bring
    in pbps. Usage: reccy_finder.py  <gubbins_res_locs> <pbp_profiles> <gff_csv> <contig_bounds> <out_name> '''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='reccy_detector.py')

    parser.add_argument('--gubbins_res', required=True, help='Directory where all cluster dirs of gubbins res are stored"', type=str)
    parser.add_argument('--pbp_profiles', required=True, help='Hit csv file from hit allocator', type=str)
    parser.add_argument('--gff_csv', required=True, help="Location of the gff csv containing isolate and reference locations", type=str)
    parser.add_argument('--contig_bounds', required=True, help="contig bound location", type=str)
    parser.add_argument('--out_name', required=True, help= 'Prefix to append to out out_files', type = str)

    args = parser.parse_args()

    return args

def search_function_2(value, df, col_name):
    current_col = df[col_name]
    indexio = current_col.loc[current_col == value].index
    if len(indexio) > 1:
        indexio = indexio[0]
    if indexio.size > 0:
        row_indy = indexio

        return row_indy

    else:
        return "Not_here"

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start + len(needle))
        n -= 1
    return start

def previous_node_finder(tree, example_id, current_node):
    node_seq = []


    tip_node = tree.find_node_with_taxon_label(label=example_id)
    parent_node = tip_node.parent_node
    node_seq.append(example_id)

    while parent_node != None:
        new_parent = parent_node.parent_node
        old_parent = parent_node
        parent_node = new_parent

        node_seq.append(old_parent.label)

    current_node_index = node_seq.index(current_node)
    new_node = node_seq[current_node_index + 1]



    return new_node

def leaf_tip_nums(tree, example_id, current_node):
    tip_node = tree.find_node_with_taxon_label(label=example_id)
    parent_node = tip_node.parent_node

    if parent_node.label == current_node:
        tip_numbers = len(parent_node.leaf_nodes())
    else:
        while parent_node.label != current_node:
            parent_node = parent_node.parent_node

            if parent_node.label == current_node:
                tip_numbers = len(parent_node.leaf_nodes())

    return tip_numbers

def length_checker(rows):
    lengos = []
    for k in range(len(rows.index)):
        length = rows.iloc[k, 4] - rows.iloc[k, 3]
        lengos.append(length)

    row_to_keep = lengos.index(max(lengos))

    df_to_keep = rows.iloc[[row_to_keep]]

    return df_to_keep

def state_finder(node_id, character_file):
    tree_states_node = character_file[character_file['node'] == node_id]
    return tree_states_node['pbp'].values[0]

def node_reconstruct(tree_loc, hit_csv):
    ## Function to reconstruct the insertion node of the individual cluster for a specific gps cluster


    tree = dendropy.Tree.get(path=tree_loc,
                             schema="newick", preserve_underscores=True)


    ## hit df ##
    cluster_csv = hit_csv.copy()
    cluster_csv = cluster_csv.reset_index(drop=True)

    # print(tree.taxon_namespace)
    ## Now we'll create the mapping of the tree isolates to the cluster nums, with 0 being no hit.
    cluster_num_col = ["ND"]
    tree_names = ["start"]
    for taxon in tree.taxon_namespace:
        current_name = taxon.label
        current_name_string = str(current_name)
        indy = search_function_2(current_name_string, cluster_csv, 'isolate')
        if isinstance(indy, str):
            multi_hit_name = current_name_string + "_1"
            indy_2 = search_function_2(multi_hit_name, cluster_csv, 'isolate')
            if isinstance(indy_2, str):
                cluster_num_col.append("")
            else:
                cluster_val = cluster_csv.iloc[indy_2, cluster_csv.columns.get_loc("profile")]
                cluster_num_col.append(cluster_val)
        else:

            cluster_val = cluster_csv.iloc[indy.values[0], cluster_csv.columns.get_loc("profile")]
            cluster_num_col.append(cluster_val)
        tree_names.append(current_name_string)

        node_changer = tree.find_node_with_taxon_label(label=current_name)
        node_changer.taxon.label = dendropy.Taxon(label=current_name_string)

    cluster_num_col = cluster_num_col[1:]
    tree_names = tree_names[1:]

    ## Heres the dictionary with our insertion numbers and the tree tip ids
    mega_characters = dict(zip(tree_names, cluster_num_col))

    with open("./pbp_fasta_data.tsv", 'w') as outfile:
        csv_writer = csv.writer(outfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(["id"] + ["pbp"])
        for k, v in mega_characters.items():
            csv_writer.writerow([k] + [v])


    data_file = "./pbp_fasta_data.tsv"
    columns = ['pbp']

    tic_ml = time.perf_counter()
    pastml_pipeline(data=data_file, columns=columns, name_column='pbp',
                    tree=tree_loc, verbose=False,
                    out_data="./pbp_states_res.txt",
                    work_dir="./")
    toc_ml = time.perf_counter()

    print("ML recon took this long: %s (seconds)" % (round(toc_ml - tic_ml)))

    tree_states = pandas.read_csv("./pbp_states_res.txt", sep="\t")

    # SeqIO.convert("./fasta_data.tsv",
    #               "tab", "./fasta_data.fasta",
    #               "fasta")
    #
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(path=tree_loc,
                             schema="newick", preserve_underscores=True,
                             taxon_namespace = taxa)
    #
    # data_mega = dendropy.StandardCharacterMatrix.get_from_path("./fasta_data.fasta",
    #                                                            "fasta", taxon_namespace=taxa)
    # taxon_state_sets_map = data_mega.taxon_state_sets_map(gaps_as_missing=True)
    #
    # tree = dendropy.Tree.get_from_path(tree_loc,
    #                                    schema="newick", preserve_underscores=True, taxon_namespace=taxa)

    #os.remove("./fasta_data.tsv")
    #os.remove("./fasta_data.fasta")

    # score = fitch_down_pass(postorder_nodes=tree.postorder_node_iter(),
    #                         taxon_state_sets_map=taxon_state_sets_map)
    # fitch_up_pass(tree.postorder_node_iter())

    ###########################################################################
    ## Ok so now the tree has been reconstructed with the fitch up pass down ##
    ## pass algorithm, we'll go through the cluster_csv (our hit csv) and    ##
    ## append the node reconstructed for insertion ############################
    ###########################################################################


    node_nums = []
    init_isolate = []
    cluster_name = []
    node_balance = []
    MGE = []
    non_MGE = []
    objective_balance = []
    tree_ids = []
    node_home = []
    reccy_on_node = []
    previous_state = []

    for row in range(len(cluster_csv.index)):
        current_row = cluster_csv.iloc[row]
        current_id = str(current_row['isolate'])
        current_clust_name = str(current_row['cluster'])
        if current_id.count("_") >= 2:
            second_occy = find_nth(current_id, "_", 2)
            current_id = current_id[0:second_occy]

        tree_id = current_id




        tree_node = tree.find_node_with_taxon_label(label=tree_id)
        if tree_node == None:
            print("Tree labels don't match isolate ids from fastas")
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            exit()

        taxon_state = state_finder(tree_id, tree_states)

        parent_node = tree_node.parent_node
        init_run = 0
        parent_state = state_finder(parent_node.label, tree_states)


        # if current_id == "22841_3#12":
        #     print(parent_state)

        if parent_state != taxon_state:
            children = 1
            tot_children = len(parent_node.leaf_nodes())
            newer_parent = parent_node
            non_children = tot_children - children
            node_home.append(newer_parent.label)
            reccy_on_node.append("No")
            previous_state.append(parent_state)
        else:
            while parent_state == taxon_state:
                if init_run == 0:
                    new_parent = parent_node.parent_node
                    old_parent = parent_node
                    init_run += 1
                elif init_run > 0:
                    new_parent = newer_parent.parent_node
                    old_parent = newer_parent

                newer_parent = new_parent

                if newer_parent == None:
                    parent_state = 30000
                    newer_parent = old_parent
                    children = len(old_parent.leaf_nodes())
                    tot_children = children
                    non_children = tot_children - children
                    node_home.append(old_parent.label)
                    reccy_on_node.append("Yes")
                    previous_state.append(taxon_state)
                else:
                    parent_state = state_finder(newer_parent.label, tree_states)
                    if parent_state != taxon_state:
                        children = len(old_parent.leaf_nodes())
                        tot_children = len(newer_parent.leaf_nodes())
                        non_children = tot_children - children
                        node_home.append(old_parent.label)
                        reccy_on_node.append("Yes")
                        previous_state.append(parent_state)

        if newer_parent.label not in node_nums or current_clust_name not in cluster_name:
            node_nums.append(newer_parent.label)
            init_isolate.append(current_id)
            cluster_name.append(current_clust_name)
            node_balance.append(children / tot_children)
            MGE.append(children)
            non_MGE.append(tot_children)
            stat_val = (abs(children - non_children) / tot_children)
            objective_balance.append(stat_val)

    out_balance_csv = pandas.DataFrame({'isolate_id': init_isolate,
                                        'hit_loc': cluster_name,
                                        'node_balance': node_balance,
                                        'MGE_isolates': MGE,
                                        'total_isolates': non_MGE,
                                        'test_stat': objective_balance,
                                        'node_id': node_nums})

    cluster_csv['insertion_node'] = pandas.Series(node_home, index=cluster_csv.index)
    cluster_csv['node_rec'] = pandas.Series(reccy_on_node, index=cluster_csv.index)
    cluster_csv['previous_state'] = pandas.Series(previous_state, index=cluster_csv.index)

    return cluster_csv

def reference_gff_clean(ref_gff_loc):
    ## Function to clean up the recombination prediction gff file to allow for parsing through

    ## Load up the file
    ref_gff_tsv = pandas.read_csv(ref_gff_loc, sep='\t',
                                  names=["type", "prog", "class", "start", "end", "trent", "alexander", "arnold",
                                         "attribute"],
                                  header=None)
    ## Set up the output
    start_nodes = []
    end_nodes = []
    example_taxa = []
    snipper = []

    ## Loop through the rows

    for k in range(len(ref_gff_tsv.index)):
        ## Tidy up the node names
        attribute_row = ref_gff_tsv.iloc[k]['attribute']
        attribute_row = str(attribute_row)
        nodes = attribute_row.split(";neg")[0]
        nodes_now = nodes.split("=\"")[-1]
        nodes_both = nodes_now.strip("\"")
        start_node = nodes_both.split("->")[0]
        end_node = nodes_both.split("->")[-1]

        ## Now we'll get the first taxa object

        taxa = attribute_row.split("taxa=")[-1]
        taxa_no_leading_whitespace = re.sub("^\\D*", "", taxa)
        taxa_narrow = taxa_no_leading_whitespace.split("\";snp_count")[0]
        taxa_narrower = taxa_narrow.split("  ")[0]
        if taxa_narrower.count("_") > 1:
            taxa_even_narrower = taxa_narrower.split("_.")[0]
            taxa_narrower = taxa_even_narrower
        new_name = taxa_narrower
        example_taxa.append(new_name)

        ## Now we'll get the snp count

        snippers = attribute_row.split(";snp_count=")[-1]
        snippers_2 = re.sub(";", "", snippers)
        snippers_3 = re.sub("\"", "", snippers_2)

        snipper.append(snippers_3)

        # if

        start_nodes.append(start_node)
        end_nodes.append(end_node)

    ref_gff_tsv['start_node'] = pandas.Series(data=start_nodes, index=ref_gff_tsv.index)
    ref_gff_tsv['end_node'] = pandas.Series(data=end_nodes, index=ref_gff_tsv.index)
    ref_gff_tsv['example_id'] = pandas.Series(data=example_taxa, index=ref_gff_tsv.index)
    ref_gff_tsv['snp_count'] = pandas.Series(data=snipper, index=ref_gff_tsv.index)

    return ref_gff_tsv

def reccy_finder(ref_gff, closest_vals, node_label, example_id, node_rec, tree):
    new_row = []
    start_insert = []
    end_insert = []
    start_node = []
    end_node = []
    start_gubb = []
    end_gubb = []
    gubb_length = []
    isolate_idz = []
    snippy_count = []
    insertion_node = []
    missing_gubbins = []

    midpoint = (closest_vals[1] + closest_vals[0]) / 2


    if node_rec == "Yes":
        end_node_subset_gff = ref_gff[(ref_gff['end_node'] == node_label)]


        if end_node_subset_gff.empty:

            ## We'll check if the previous node has any reccys now

            example_isolate = example_id

            if "ROOT" in node_label:
                out_df = pandas.DataFrame()
                out_df['start_insert'] = pandas.Series(data=closest_vals[0])
                out_df['end_insert'] = pandas.Series(data=closest_vals[1], index=out_df.index)
                out_df['isolate_id'] = pandas.Series(data=example_id, index=out_df.index)
                out_df['insertion_node'] = pandas.Series(data=node_label, index=out_df.index)

                return out_df, False


            new_node = previous_node_finder(tree, example_isolate, node_label)

            new_end_node_subset = ref_gff[ref_gff['end_node'] == new_node]
            new_end_node_subset_loccy = new_end_node_subset[(new_end_node_subset['start'] <= closest_vals[0]) & (new_end_node_subset['end'] >= closest_vals[1])]

            if new_end_node_subset_loccy.empty:
                ## Lets try with a slightly larger margin around the closer vals
                new_end_node_subset_loccy = new_end_node_subset[(new_end_node_subset['start'] <= midpoint) & (new_end_node_subset['end'] >= midpoint)]

                if new_end_node_subset_loccy.empty:
                    new_end_node_subset_loccy = new_end_node_subset[
                        (new_end_node_subset['start'] <= (closest_vals[0] + 100)) & (
                                new_end_node_subset['end'] >= (closest_vals[1] - 100))]

                if new_end_node_subset_loccy.empty:
                    #print("Definitely no reccy here 1")
                    I_am_lorde = "lalalaaa"
                else:
                    if len(new_end_node_subset_loccy.index) > 1:
                        single_row = length_checker(new_end_node_subset_loccy)
                    else:
                        single_row = new_end_node_subset_loccy
                    start_insert.append(closest_vals[0])
                    end_insert.append(closest_vals[1])
                    start_node.append(single_row.iloc[0, 9])
                    end_node.append(single_row.iloc[0, 10])
                    start_gubb.append(single_row.iloc[0, 3])
                    end_gubb.append(single_row.iloc[0, 4])
                    gubb_length.append(end_gubb[-1] - start_gubb[0])
                    isolate_idz.append(example_id)
                    snippy_count.append(single_row.iloc[0, 12])
                    insertion_node.append(node_label)




            else:
                if len(new_end_node_subset_loccy.index) > 1:
                    single_row = length_checker(new_end_node_subset_loccy)
                else:
                    single_row = new_end_node_subset_loccy
                start_insert.append(closest_vals[0])
                end_insert.append(closest_vals[1])
                start_node.append(single_row.iloc[0, 9])
                end_node.append(single_row.iloc[0, 10])
                start_gubb.append(single_row.iloc[0, 3])
                end_gubb.append(single_row.iloc[0, 4])
                gubb_length.append(end_gubb[-1] - start_gubb[0])
                isolate_idz.append(example_id)
                snippy_count.append(single_row.iloc[0, 12])
                insertion_node.append(node_label)


        else:
            end_node_subset_loc = end_node_subset_gff[(end_node_subset_gff['start'] <= closest_vals[0]) & (end_node_subset_gff['end'] >= closest_vals[1])]


            if end_node_subset_loc.empty:
                ## Lets use a more lenient val range now then if this doesn't work we'll use the previous node

                new_end_node_subset_loccy = end_node_subset_gff[(end_node_subset_gff['start'] <= midpoint) & (end_node_subset_gff['end'] >= midpoint)]

                if new_end_node_subset_loccy.empty:
                    new_end_node_subset_loccy = end_node_subset_loc[
                        (end_node_subset_loc['start'] <= (closest_vals[0] + 100)) & (
                                end_node_subset_loc['end'] >= (closest_vals[1] - 100))]

                if new_end_node_subset_loccy.empty:

                    # print("No detected reccy in this region", node_label, closest_vals)

                    ## Now we'll look at one node furtherup as sometimes we see these recombination events
                    example_isolate = example_id

                    if "ROOT" in node_label:
                        out_df = pandas.DataFrame()
                        out_df['start_insert'] = pandas.Series(data=closest_vals[0])
                        out_df['end_insert'] = pandas.Series(data=closest_vals[1], index=out_df.index)
                        out_df['isolate_id'] = pandas.Series(data=example_id, index=out_df.index)
                        out_df['insertion_node'] = pandas.Series(data=node_label, index=out_df.index)

                        return out_df, False

                    new_node = previous_node_finder(tree, example_isolate, node_label)

                    new_end_node_subset = ref_gff[ref_gff['end_node'] == new_node]
                    new_end_node_subset_loccy = new_end_node_subset[
                        (new_end_node_subset['start'] <= closest_vals[0]) & (
                                    new_end_node_subset['end'] >= closest_vals[1])]

                    if new_end_node_subset_loccy.empty:
                        ## Lets try a more lenient window
                        new_end_node_subset_loccy = new_end_node_subset[
                            (new_end_node_subset['start'] <= midpoint) & (
                                    new_end_node_subset['end'] >= midpoint)]

                        if new_end_node_subset_loccy.empty:
                            new_end_node_subset_loccy = new_end_node_subset[
                                (new_end_node_subset['start'] <= (closest_vals[0] + 100)) & (
                                        new_end_node_subset['end'] >= (closest_vals[1] - 100))]

                        if new_end_node_subset_loccy.empty:
                            #print("Definitely no reccy here 2", new_node, node_label)
                            feeling_good = "on a wednesday"
                        else:
                            if len(new_end_node_subset_loccy.index) > 1:
                                single_row = length_checker(new_end_node_subset_loccy)
                            else:
                                single_row = new_end_node_subset_loccy
                            start_insert.append(closest_vals[0])
                            end_insert.append(closest_vals[1])
                            start_node.append(single_row.iloc[0, 9])
                            end_node.append(single_row.iloc[0, 10])
                            start_gubb.append(single_row.iloc[0, 3])
                            end_gubb.append(single_row.iloc[0, 4])
                            gubb_length.append(end_gubb[-1] - start_gubb[0])
                            isolate_idz.append(example_id)
                            snippy_count.append(single_row.iloc[0, 12])
                            insertion_node.append(node_label)


                else:
                    if len(new_end_node_subset_loccy.index) > 1:
                        single_row = length_checker(new_end_node_subset_loccy)
                    else:
                        single_row = new_end_node_subset_loccy
                    start_insert.append(closest_vals[0])
                    end_insert.append(closest_vals[1])
                    start_node.append(single_row.iloc[0, 9])
                    end_node.append(single_row.iloc[0, 10])
                    start_gubb.append(single_row.iloc[0, 3])
                    end_gubb.append(single_row.iloc[0, 4])
                    gubb_length.append(end_gubb[-1] - start_gubb[0])
                    isolate_idz.append(example_id)
                    snippy_count.append(single_row.iloc[0, 12])
                    insertion_node.append(node_label)

            else:
                if len(end_node_subset_loc.index) > 1:
                    single_row = length_checker(end_node_subset_loc)
                else:
                    single_row = end_node_subset_loc
                start_insert.append(closest_vals[0])
                end_insert.append(closest_vals[1])
                start_node.append(single_row.iloc[0, 9])
                end_node.append(single_row.iloc[0, 10])
                start_gubb.append(single_row.iloc[0, 3])
                end_gubb.append(single_row.iloc[0, 4])
                gubb_length.append(end_gubb[-1] - start_gubb[0])
                isolate_idz.append(example_id)
                snippy_count.append(single_row.iloc[0, 12])
                insertion_node.append(node_label)
    elif node_rec == "No":
            end_node_subset_gff = ref_gff[(ref_gff['start_node'] == node_label)]

            if end_node_subset_gff.empty:
                new_end_node_subset = ref_gff[(ref_gff['end_node'] == node_label)]

                new_end_node_subset_loccy = new_end_node_subset[(new_end_node_subset['start'] <= closest_vals[0]) & (new_end_node_subset['end'] >= closest_vals[1])]

                if new_end_node_subset_loccy.empty:
                    new_end_node_subset_loccy = new_end_node_subset[
                        (new_end_node_subset['start'] <= midpoint) & (
                                    new_end_node_subset['end'] >= midpoint)]

                    if new_end_node_subset_loccy.empty:
                        new_end_node_subset_loccy = new_end_node_subset[
                            (new_end_node_subset['start'] <= (closest_vals[0] + 100)) & (
                                    new_end_node_subset['end'] >= (closest_vals[1] - 100))]

                    if new_end_node_subset_loccy.empty:
                        #print("No defo no reccy here")
                        lorde = "lorde lorde"
                    else:
                        if len(new_end_node_subset_loccy.index) > 1:
                            single_row = length_checker(new_end_node_subset_loccy)
                        else:
                            single_row = new_end_node_subset_loccy
                        start_insert.append(closest_vals[0])
                        end_insert.append(closest_vals[1])
                        start_node.append(single_row.iloc[0, 9])
                        end_node.append(single_row.iloc[0, 10])
                        start_gubb.append(single_row.iloc[0, 3])
                        end_gubb.append(single_row.iloc[0, 4])
                        gubb_length.append(end_gubb[-1] - start_gubb[0])
                        isolate_idz.append(example_id)
                        snippy_count.append(single_row.iloc[0, 12])
                        insertion_node.append(node_label)


                else:
                    if len(new_end_node_subset_loccy.index) > 1:
                        single_row = length_checker(new_end_node_subset_loccy)
                    else:
                        single_row = new_end_node_subset_loccy
                    start_insert.append(closest_vals[0])
                    end_insert.append(closest_vals[1])
                    start_node.append(single_row.iloc[0, 9])
                    end_node.append(single_row.iloc[0, 10])
                    start_gubb.append(single_row.iloc[0, 3])
                    end_gubb.append(single_row.iloc[0, 4])
                    gubb_length.append(end_gubb[-1] - start_gubb[0])
                    isolate_idz.append(example_id)
                    snippy_count.append(single_row.iloc[0, 12])
                    insertion_node.append(node_label)


            else:

                end_node_subset_loc = end_node_subset_gff[
                    (end_node_subset_gff['start'] <= closest_vals[0]) & (end_node_subset_gff['end'] >= closest_vals[1])]
                if end_node_subset_loc.empty:

                    end_node_subset_loc = end_node_subset_gff[(end_node_subset_gff['start'] <= midpoint) & (
                                end_node_subset_gff['end'] >= midpoint)]

                    if end_node_subset_loc.empty:
                        end_node_subset_loc = end_node_subset_gff[(end_node_subset_gff['start'] <= (closest_vals[0] + 100)) & (
                                end_node_subset_gff['end'] >= (closest_vals[1] - 100))]

                    if end_node_subset_loc.empty:


                        new_end_node_subset = ref_gff[(ref_gff['end_node'] == node_label)]

                        new_end_node_subset_loccy = new_end_node_subset[
                            (new_end_node_subset['start'] <= closest_vals[0]) & (
                                        new_end_node_subset['end'] >= closest_vals[1])]

                        if new_end_node_subset_loccy.empty:
                            new_end_node_subset_loccy = new_end_node_subset[
                                (new_end_node_subset['start'] <= midpoint) & (
                                            new_end_node_subset['end'] >= midpoint)]

                            if new_end_node_subset_loccy.empty:
                                new_end_node_subset_loccy = new_end_node_subset[
                                    (new_end_node_subset['start'] <= (closest_vals[0] + 100)) & (
                                            new_end_node_subset['end'] >= (closest_vals[1] - 100))]

                            if new_end_node_subset_loccy.empty:
                                #print("No defo no reccy here", node_label)
                                yeah_yeah = "yeaah I am Lorde"
                            else:
                                if len(new_end_node_subset_loccy.index) > 1:
                                    single_row = length_checker(new_end_node_subset_loccy)
                                else:
                                    single_row = new_end_node_subset_loccy
                                start_insert.append(closest_vals[0])
                                end_insert.append(closest_vals[1])
                                start_node.append(single_row.iloc[0, 9])
                                end_node.append(single_row.iloc[0, 10])
                                start_gubb.append(single_row.iloc[0, 3])
                                end_gubb.append(single_row.iloc[0, 4])
                                gubb_length.append(end_gubb[-1] - start_gubb[0])
                                isolate_idz.append(example_id)
                                snippy_count.append(single_row.iloc[0, 12])
                                insertion_node.append(node_label)




                        else:
                            if len(new_end_node_subset_loccy.index) > 1:
                                single_row = length_checker(new_end_node_subset_loccy)
                            else:
                                single_row = new_end_node_subset_loccy
                            start_insert.append(closest_vals[0])
                            end_insert.append(closest_vals[1])
                            start_node.append(single_row.iloc[0, 9])
                            end_node.append(single_row.iloc[0, 10])
                            start_gubb.append(single_row.iloc[0, 3])
                            end_gubb.append(single_row.iloc[0, 4])
                            gubb_length.append(end_gubb[-1] - start_gubb[0])
                            isolate_idz.append(example_id)
                            snippy_count.append(single_row.iloc[0, 12])
                            insertion_node.append(node_label)
                    else:

                        if len(end_node_subset_loc.index) > 1:
                            single_row = length_checker(end_node_subset_loc)
                        else:
                            single_row = end_node_subset_loc

                        start_insert.append(closest_vals[0])
                        end_insert.append(closest_vals[1])
                        start_node.append(single_row.iloc[0, 9])
                        end_node.append(single_row.iloc[0, 10])
                        start_gubb.append(single_row.iloc[0, 3])
                        end_gubb.append(single_row.iloc[0, 4])
                        gubb_length.append(end_gubb[-1] - start_gubb[0])
                        isolate_idz.append(example_id)
                        snippy_count.append(single_row.iloc[0, 12])
                        insertion_node.append(node_label)

                else:

                    if len(end_node_subset_loc.index) > 1:
                        single_row = length_checker(end_node_subset_loc)
                    else:
                        single_row = end_node_subset_loc

                    start_insert.append(closest_vals[0])
                    end_insert.append(closest_vals[1])
                    start_node.append(single_row.iloc[0, 9])
                    end_node.append(single_row.iloc[0, 10])
                    start_gubb.append(single_row.iloc[0, 3])
                    end_gubb.append(single_row.iloc[0, 4])
                    gubb_length.append(end_gubb[-1] - start_gubb[0])
                    isolate_idz.append(example_id)
                    snippy_count.append(single_row.iloc[0, 12])
                    insertion_node.append(node_label)
    if len(end_gubb) == 1:

        out_df = pandas.DataFrame()
        out_df['start_insert'] = pandas.Series(data=start_insert)
        out_df['end_insert'] = pandas.Series(data=end_insert, index=out_df.index)
        out_df['start_nodes'] = pandas.Series(data=start_node, index=out_df.index)
        out_df['end_nodes'] = pandas.Series(data=end_node, index=out_df.index)
        out_df['start_gubb'] = pandas.Series(data=start_gubb, index=out_df.index)
        out_df['end_gubb'] = pandas.Series(data=end_gubb, index=out_df.index)
        out_df['gub_length'] = pandas.Series(data=gubb_length, index=out_df.index)
        out_df['isolate_example'] = pandas.Series(data=isolate_idz, index=out_df.index)
        out_df['snp_count'] = pandas.Series(data=snippy_count, index=out_df.index)
        out_df['insertion_node'] = pandas.Series(data=insertion_node, index=out_df.index)

        return out_df, True

    else:
        out_df = pandas.DataFrame()
        out_df['start_insert'] = pandas.Series(data=closest_vals[0])
        out_df['end_insert'] = pandas.Series(data=closest_vals[1], index=out_df.index)
        out_df['isolate_id'] = pandas.Series(data=example_id, index=out_df.index)
        out_df['insertion_node'] = pandas.Series(data=node_label, index=out_df.index)

        return  out_df, False

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

def reccy_main(tree_loc, ref_gff_tsv, hit_csv, gene_row, gene_name, ):
    ## Function to search through the hit_csv, which now has the insertion nodes of its cluster attached,
    ## and find out if the hits are within a reconstructed recombination event (hence likely to be transformation)


    tree = dendropy.Tree.get(path=tree_loc,
                             schema="newick", preserve_underscores=True)

    tree_tip_names = []

    for taxon in tree.taxon_namespace:
        current_name = taxon.label
        tree_tip_names.append(current_name)

    ###############################################################################
    ## Now we'll look through the python hits csv to find the nodes we need to ####
    ## look into ##################################################################
    ###############################################################################

    ## This is the hit csv with attached insertion node data
    python_hits_csv = hit_csv.copy()

    node_count = python_hits_csv['insertion_node'].value_counts()


    # print(len(node_count))

    # sorted(product(arr1, arr2), key=lambda t: abs(t[0]-t[1]))[0]
    reccy_df_columns = ['start_insert', 'end_insert', 'start_nodes', 'end_nodes', 'start_gubb',
                        'end_gubb', 'gub_length', 'isolate_example', 'snp_count', 'insertion_node',
                        'insertion_node_tip_nums', 'end_node_tip_nums', 'resistance', 'previous_resistance']

    reccy_df = pandas.DataFrame(columns=reccy_df_columns)


    start_inserts = []
    end_inserts = []
    start_nodes = []
    end_nodes = []
    start_gubbs = []
    end_gubbs = []
    gubb_lengths = []
    isolate_idzs = []
    snippy_counts = []
    insertion_nodes = []
    inny_node_tips = []
    end_node_tips = []
    resistance = []
    prev_resistance = []

    non_reccy_df_columns = ['start_insert', 'end_insert', 'isolate_id', 'insertion_node', 'insertion_node_tip_nums',
                            'resistance', 'previous_resistance']
    non_reccy_df = pandas.DataFrame(columns=non_reccy_df_columns)

    non_start_insert = []
    non_end_insert = []
    non_example_isolate = []
    non_insertion_node = []
    non_inny_node_tips = []
    non_inny_insert_name = []
    non_inny_prev_res = []

    for k in range(len(node_count)):
        # print(k)
        node_label = node_count.index[k]
        noders = python_hits_csv[python_hits_csv['insertion_node'] == node_label]


        check_if_same_insertion_loci = noders['cluster'].value_counts()


        if len(check_if_same_insertion_loci) == 1:
            noders_row = noders.iloc[0]
            before_locs = gene_row.iloc[0,3]
            after_locs = gene_row.iloc[0, 4]
            example_id = noders.iloc[0, 0]
            previous_state = noders['previous_state'].iloc[0]

            closest_vals = [before_locs, after_locs]
            closest_vals = sorted(closest_vals)

            ## First we check if we should use the node as the end node or the start node
            if noders_row['node_rec'] == "Yes":
                new_row, add_row = reccy_finder(ref_gff_tsv, closest_vals, node_label, example_id,
                                                "Yes", tree)

                if add_row == True:
                    start_inserts.append(new_row.iloc[0, 0])
                    end_inserts.append(new_row.iloc[0, 1])
                    start_nodes.append(new_row.iloc[0, 2])
                    end_nodes.append(new_row.iloc[0, 3])
                    start_gubbs.append(new_row.iloc[0, 4])
                    end_gubbs.append(new_row.iloc[0, 5])
                    gubb_lengths.append(new_row.iloc[0, 6])
                    isolate_idzs.append(new_row.iloc[0, 7])
                    snippy_counts.append(new_row.iloc[0, 8])
                    insertion_nodes.append(new_row.iloc[0, 9])
                    inny_node_tips.append(leaf_tip_nums(tree, example_id, current_node=insertion_nodes[-1]))
                    end_node_tips.append(leaf_tip_nums(tree, example_id, current_node=end_nodes[-1]))
                    resistance.append(noders.iloc[0, noders.columns.get_loc('profile')])
                    prev_resistance.append(previous_state)




                else:
                    non_start_insert.append(new_row.iloc[0, 0])
                    non_end_insert.append(new_row.iloc[0, 1])
                    non_example_isolate.append(new_row.iloc[0, 2])
                    non_insertion_node.append(new_row.iloc[0, 3])
                    non_inny_node_tips.append(leaf_tip_nums(tree, example_id, current_node=non_insertion_node[-1]))
                    non_inny_insert_name.append(noders.iloc[0, noders.columns.get_loc('profile')])
                    non_inny_prev_res.append(previous_state)


            else:
                new_row, add_row = reccy_finder(ref_gff_tsv, closest_vals, node_label, example_id,
                                                "No", tree)
                if add_row == True:
                    start_inserts.append(new_row.iloc[0, 0])
                    end_inserts.append(new_row.iloc[0, 1])
                    start_nodes.append(new_row.iloc[0, 2])
                    end_nodes.append(new_row.iloc[0, 3])
                    start_gubbs.append(new_row.iloc[0, 4])
                    end_gubbs.append(new_row.iloc[0, 5])
                    gubb_lengths.append(new_row.iloc[0, 6])
                    isolate_idzs.append(new_row.iloc[0, 7])
                    snippy_counts.append(new_row.iloc[0, 8])
                    insertion_nodes.append(new_row.iloc[0, 9])
                    inny_node_tips.append(leaf_tip_nums(tree, example_id, current_node=insertion_nodes[-1]))
                    end_node_tips.append(leaf_tip_nums(tree, example_id, current_node=start_nodes[-1]))
                    resistance.append(noders.iloc[0, noders.columns.get_loc('profile')])
                    prev_resistance.append(previous_state)





                else:
                    non_start_insert.append(new_row.iloc[0, 0])
                    non_end_insert.append(new_row.iloc[0, 1])
                    non_example_isolate.append(new_row.iloc[0, 2])
                    non_insertion_node.append(new_row.iloc[0, 3])
                    non_inny_node_tips.append(leaf_tip_nums(tree, example_id, current_node=non_insertion_node[-1]))
                    non_inny_insert_name.append(noders.iloc[0, noders.columns.get_loc('profile')])
                    non_inny_prev_res.append(previous_state)

            before_loc_list = before_locs.tolist()
            after_loc_list = after_locs.tolist()

            closest_vals = sorted(closest_vals)
            non_start_insert.append(closest_vals[0])
            non_end_insert.append(closest_vals[1])
            non_example_isolate.append(example_id)
            non_insertion_node.append(node_label)
            non_inny_node_tips.append(leaf_tip_nums(tree, example_id, current_node=non_insertion_node[-1]))
            non_inny_insert_name.append(noders.iloc[0, noders.columns.get_loc('profile')])
            non_inny_prev_res.append(previous_state)
            ## new row here

        else:
            for k in range(len(check_if_same_insertion_loci)):
                current_name = check_if_same_insertion_loci.index[k]
                subset_noders = noders[noders['cluster_names'] == current_name]
                noders_row = subset_noders.iloc[0]
                before_locs = gene_row.iloc[0, 3]
                after_locs = gene_row.iloc[0, 4]
                example_id = subset_noders.iloc[0, 0]
                previous_state = subset_noders['previous_state'].iloc[0]

                closest_vals = [before_locs, after_locs]
                closest_vals = sorted(closest_vals)


                ## First we check if we should use the node as the end node or the start node
                if noders_row['node_rec'] == "Yes":
                    new_row, add_row = reccy_finder(ref_gff_tsv, closest_vals, node_label, example_id,
                                                    "Yes", tree)
                    if add_row == True:
                        start_inserts.append(new_row.iloc[0, 0])
                        end_inserts.append(new_row.iloc[0, 1])
                        start_nodes.append(new_row.iloc[0, 2])
                        end_nodes.append(new_row.iloc[0, 3])
                        start_gubbs.append(new_row.iloc[0, 4])
                        end_gubbs.append(new_row.iloc[0, 5])
                        gubb_lengths.append(new_row.iloc[0, 6])
                        isolate_idzs.append(new_row.iloc[0, 7])
                        snippy_counts.append(new_row.iloc[0, 8])
                        insertion_nodes.append(new_row.iloc[0, 9])
                        inny_node_tips.append(leaf_tip_nums(tree, example_id, current_node=insertion_nodes[-1]))
                        end_node_tips.append(leaf_tip_nums(tree, example_id, current_node=end_nodes[-1]))
                        resistance.append(noders.iloc[0, noders.columns.get_loc('profile')])
                        prev_resistance.append(previous_state)





                    else:
                        non_start_insert.append(new_row.iloc[0, 0])
                        non_end_insert.append(new_row.iloc[0, 1])
                        non_example_isolate.append(new_row.iloc[0, 2])
                        non_insertion_node.append(new_row.iloc[0, 3])
                        non_inny_node_tips.append(leaf_tip_nums(tree, example_id, current_node=non_insertion_node[-1]))
                        non_inny_insert_name.append(noders.iloc[0, noders.columns.get_loc('profile')])
                        non_inny_prev_res.append(previous_state)


                else:
                    new_row, add_row = reccy_finder(ref_gff_tsv, closest_vals, node_label, example_id,
                                                    "No", tree)
                    if add_row == True:
                        start_inserts.append(new_row.iloc[0, 0])
                        end_inserts.append(new_row.iloc[0, 1])
                        start_nodes.append(new_row.iloc[0, 2])
                        end_nodes.append(new_row.iloc[0, 3])
                        start_gubbs.append(new_row.iloc[0, 4])
                        end_gubbs.append(new_row.iloc[0, 5])
                        gubb_lengths.append(new_row.iloc[0, 6])
                        isolate_idzs.append(new_row.iloc[0, 7])
                        snippy_counts.append(new_row.iloc[0, 8])
                        insertion_nodes.append(new_row.iloc[0, 9])
                        inny_node_tips.append(leaf_tip_nums(tree, example_id, current_node=insertion_nodes[-1]))
                        end_node_tips.append(leaf_tip_nums(tree, example_id, current_node=start_nodes[-1]))
                        resistance.append(noders.iloc[0, noders.columns.get_loc('profile')])
                        prev_resistance.append(previous_state)


                    else:
                        non_start_insert.append(new_row.iloc[0, 0])
                        non_end_insert.append(new_row.iloc[0, 1])
                        non_example_isolate.append(new_row.iloc[0, 2])
                        non_insertion_node.append(new_row.iloc[0, 3])
                        non_inny_node_tips.append(leaf_tip_nums(tree, example_id, current_node=non_insertion_node[-1]))
                        non_inny_insert_name.append(noders.iloc[0, noders.columns.get_loc('profile')])
                        non_inny_prev_res.append(previous_state)


    gubb_lengths = [x + 1 for x in gubb_lengths]



    if len(start_inserts) > 0:

        reccy_df['start_insert'] = pandas.Series(data=start_inserts)
        reccy_df['end_insert'] = pandas.Series(data=end_inserts, index=reccy_df.index)
        reccy_df['start_nodes'] = pandas.Series(data=start_nodes, index=reccy_df.index)
        reccy_df['end_nodes'] = pandas.Series(data=end_nodes, index=reccy_df.index)
        reccy_df['start_gubb'] = pandas.Series(data=start_gubbs, index=reccy_df.index)
        reccy_df['end_gubb'] = pandas.Series(data=end_gubbs, index=reccy_df.index)
        reccy_df['gub_length'] = pandas.Series(data=gubb_lengths, index=reccy_df.index)
        reccy_df['isolate_example'] = pandas.Series(data=isolate_idzs, index=reccy_df.index)
        reccy_df['snp_count'] = pandas.Series(data=snippy_counts, index=reccy_df.index)
        reccy_df['insertion_node'] = pandas.Series(data=insertion_nodes, index=reccy_df.index)
        reccy_df['insertion_node_tip_nums'] = pandas.Series(data=inny_node_tips, index=reccy_df.index)
        reccy_df['end_node_tip_nums'] = pandas.Series(data=end_node_tips, index=reccy_df.index)
        reccy_df['resistance'] = pandas.Series(data=resistance, index=reccy_df.index)
        reccy_df['gene'] = gene_name
        reccy_df['previous_resistance'] = pandas.Series(data=prev_resistance, index=reccy_df.index)
        reccy_df = reccy_df.drop_duplicates(subset=['start_nodes', 'end_nodes'], keep='first')

    if len(non_start_insert) > 0:
        non_reccy_df['start_insert'] = pandas.Series(data=non_start_insert)
        non_reccy_df['end_insert'] = pandas.Series(data=non_end_insert, index=non_reccy_df.index)
        non_reccy_df['isolate_id'] = pandas.Series(data=non_example_isolate, index=non_reccy_df.index)
        non_reccy_df['insertion_node'] = pandas.Series(data=non_insertion_node, index=non_reccy_df.index)
        non_reccy_df['insertion_node_tip_nums'] = pandas.Series(data=non_inny_node_tips, index=non_reccy_df.index)
        non_reccy_df['resistance'] = pandas.Series(data=non_inny_insert_name, index=non_reccy_df.index)
        non_reccy_df['gene'] = gene_name
        non_reccy_df['previous_resistance'] = pandas.Series(data=non_inny_prev_res, index=non_reccy_df.index)

    return reccy_df, non_reccy_df

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

if __name__ == '__main__':
    tic_tot = time.perf_counter()
    input_args = get_options()

    ## Loop through all the clusters to get the reccy csv for each one then append together
    ## pbp profiles should have three columns: id, profile, cluster
    pbp_profiles = pandas.read_csv(input_args.pbp_profiles)
    pbp_profiles.columns = ['isolate','profile','cluster']

    gff_csv = pandas.read_csv(input_args.gff_csv)
    gff_csv.columns = ['isolate', 'reference', 'cluster']
    base_loc = input_args.gubbins_res
    unique_clusters = pbp_profiles['cluster'].unique()
    contig_loc = input_args.contig_bounds

    tot_reccy_csv = pandas.DataFrame()
    tot_non_reccy = pandas.DataFrame()
    hit_df_new = pandas.DataFrame()

    seq_clus = 1
    for cluster in unique_clusters:
        print("On cluster: %s, %s of %s" % (cluster, seq_clus, len(unique_clusters)))
        tic_cluster = time.perf_counter()
        ## Narrow down the pbp and gff csvs to just the current cluster
        current_dat = pbp_profiles[pbp_profiles['cluster'] == cluster]
        current_gff_rows = gff_csv[gff_csv['cluster'] == cluster]
        current_gff_rows = current_gff_rows.reset_index(drop=True)

        current_ref_loc = current_gff_rows['reference'].iloc[0]
        ## get the current ref name for contig bounds
        ref_name = os.path.basename(current_ref_loc )
        ref_name = re.sub("\..*[a-zA-Z]*$", "", ref_name)
        contig_suffix = "#contig_bounds.csv"
        contig_isolate = re.sub("#", "_", ref_name)
        contig_file_path = contig_loc + contig_isolate + contig_suffix
        ref_contig = pandas.read_csv(contig_file_path)
        current_dir = base_loc + cluster
        print(current_dir)
        try:
            cluster_files = os.listdir(current_dir)
        except:
            current_dir = current_dir + "_run_data"
            cluster_files = os.listdir(current_dir)

        ## Get the tree and the recombination predictions for the current cluster
        tree_indexio = [k for k, s in enumerate(cluster_files) if "node_labelled.final_tree.tre" in s]
        reccy_index =  [k for k, s in enumerate(cluster_files) if "recombination_predictions.gff" in s]

        tree_loc = current_dir + "/" + cluster_files[tree_indexio[0]]
        gff_loc = current_dir + "/" + cluster_files[reccy_index[0]]

        current_dat = node_reconstruct(tree_loc, current_dat)

        rec_gff = reference_gff_clean(gff_loc)
        ## Now we've appended the likely gain of resistance nodes in the collection, we need
        ## to work through the three pbp genes and see if there is any recombination events
        ## around these at these nodes



        ref_gff = pandas.read_csv(current_ref_loc, sep='\t',
                                      names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                             'attributes'],
                                      header=None)
        if len(ref_contig.index) > 1:
            current_gff = gff_to_dna(ref_gff, ref_contig, ref_name, input_k=1)
            if isinstance(current_gff, str):
                continue
        ## pbp 1a
        correct_length_1a, pbp1a_row = search_for_gene(ref_gff, ['pbp1A','ponA'], 2159,
                                                     100, False, None)
        ## pbp 2b
        correct_length_2b, pbp2b_row = search_for_gene(ref_gff, ['pbp2B', 'penA'], 2042,
                                                    100, False, None)
        ## pbp 2x
        correct_length_2x, pbp2x_row = search_for_gene(ref_gff, ['pbp2X', 'pbpX'], 2252,
                                                    100, False, None)

        all_present = correct_length_1a and correct_length_2b and correct_length_2x

        if not all_present:
            missing_genes = ""
            if not correct_length_1a:
                missing_genes = missing_genes + " pbp1A "
            if not correct_length_2b:
                missing_genes = missing_genes + " pbp2B "
            if not correct_length_2x:
                missing_genes = missing_genes + " pbp2X "
            print("Missing: %s for cluster: %s" % (missing_genes, cluster))
            continue


        ## pbp 1a reccy search

        pbp1a_reccy_hits, pbp1a_reccy_misses = reccy_main(tree_loc, rec_gff, current_dat, pbp1a_row, "pbp1A")

        ## pbp 2b reccy search

        pbp2b_reccy_hits, pbp2b_reccy_misses = reccy_main(tree_loc, rec_gff, current_dat, pbp2b_row, "pbp2B")

        ## pbp2x reccy search

        pbp2x_reccy_hits, pbp2x_reccy_misses = reccy_main(tree_loc, rec_gff, current_dat, pbp2x_row, "pbp2X")

        insertion_nodes = len(current_dat['insertion_node'].unique())

        current_reccy_hits = pandas.concat([pbp1a_reccy_hits, pbp2b_reccy_hits, pbp2x_reccy_hits],
                                           axis=0, sort=False, ignore_index=True)
        current_reccy_misses = pandas.concat([pbp1a_reccy_misses, pbp2b_reccy_misses, pbp2x_reccy_misses],
                                             axis=0, sort=False, ignore_index=True)

        current_reccy_hits['cluster_name'] = cluster
        current_reccy_misses['cluster_name'] = cluster

        tot_reccy_csv = tot_reccy_csv.append(current_reccy_hits, sort= False, ignore_index=True)
        tot_non_reccy = tot_non_reccy.append(current_reccy_misses, sort = False, ignore_index=True)
        hit_df_new = hit_df_new.append(current_dat, sort= False, ignore_index=True)
        toc_cluster = time.perf_counter()
        seq_clus += 1
        print("Found %s genes in recombinations, %s outside recombination, total state alterations %s" %
              (len(current_reccy_hits.index), len(current_reccy_misses.index),insertion_nodes))
        print("Took this long for %s, %s (seconds)" % (cluster, (toc_cluster - tic_cluster)))
        #sys.exit(1)



    reccy_name = input_args.out_name + "_reccy_hits.csv"
    non_reccy_name = input_args.out_name + "_non_reccy_hits.csv"
    hits_df_name = input_args.out_name + "_hits_df.csv"

    tot_reccy_csv.to_csv(path_or_buf=reccy_name, index=False)
    tot_non_reccy.to_csv(path_or_buf=non_reccy_name, index=False)
    hit_df_new.to_csv(path_or_buf=hits_df_name, index=False)
    toc_tot = time.perf_counter()

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Reccy finder took: %s (seconds)" % (toc_tot - tic_tot))
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Reccies found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
