import dendropy
from dendropy.model.parsimony import fitch_down_pass
from dendropy.model.parsimony import fitch_up_pass
import pandas
import csv
import numpy
import sys
import re
from Bio import Phylo
from Bio import SeqIO
from statistics import mean
from statistics import stdev
from scipy.stats import wilcoxon
from scipy.stats import mannwhitneyu
import bisect
import operator
import argparse
import time
import os
###############################################################################
## Load up the trees ##########################################################
###############################################################################

def get_options():

    purpose = '''This is a python script to output each branch's mutation rate.
    Usage: node_balance_and_mutation_get.py <gubbins_dir> <hit_locs_csv> <out_csv> '''
    parser = argparse.ArgumentParser(description=purpose, prog='node_mutation.py')

    parser.add_argument('--gubbins_res', required=True, help='Directory where all cluster dirs of gubbins res are stored"', type=str)
    parser.add_argument('--hit_locs_csv', required=True, help='Hit csv file from hit allocator', type=str)
    parser.add_argument('--out_name', required=True, help='Prefix to append to out out_files', type=str)

    args = parser.parse_args()

    return args

def search_function(value, df, col_name):
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

def node_reconstruct(tree_loc, hit_csv):
    ## Function to reconstruct the insertion node of the individual cluster for a specific gps cluster


    tree = dendropy.Tree.get(path=tree_loc,
                             schema="newick", preserve_underscores=True)


    ## hit df ##
    cluster_csv = hit_csv.copy()
    cluster_csv = cluster_csv.reset_index(drop=True)

    # print(tree.taxon_namespace)
    ## Now we'll create the mapping of the tree isolates to the cluster nums, with 0 being no hit.
    cluster_num_col = [1]
    tree_names = ["start"]
    for taxon in tree.taxon_namespace:
        current_name = taxon.label
        current_name_string = str(current_name)
        indy = search_function(current_name_string, cluster_csv, 'id')
        if isinstance(indy, str):
            multi_hit_name = current_name_string + "_1"
            indy_2 = search_function(multi_hit_name, cluster_csv, 'id')
            if isinstance(indy_2, str):
                cluster_num_col.append(0)
            else:
                cluster_val = cluster_csv.iloc[indy_2, cluster_csv.columns.get_loc("insert_name")]
                cluster_num_col.append(int(cluster_val))
        else:

            cluster_val = cluster_csv.iloc[indy.values[0], cluster_csv.columns.get_loc("insert_name")]
            cluster_num_col.append(int(cluster_val))
        tree_names.append(current_name_string)

        node_changer = tree.find_node_with_taxon_label(label=current_name)
        node_changer.taxon.label = dendropy.Taxon(label=current_name_string)

    cluster_num_col = cluster_num_col[1:]
    tree_names = tree_names[1:]

    ## Heres the dictionary with our insertion numbers and the tree tip ids
    mega_characters = dict(zip(tree_names, cluster_num_col))

    with open("./fasta_data.tsv", 'w') as outfile:
        csv_writer = csv.writer(outfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        for k, v in mega_characters.items():
            csv_writer.writerow([k] + [v])

    SeqIO.convert("./fasta_data.tsv",
                  "tab", "./fasta_data.fasta",
                  "fasta")

    taxa = dendropy.TaxonNamespace()

    data_mega = dendropy.StandardCharacterMatrix.get_from_path("./fasta_data.fasta",
                                                               "fasta", taxon_namespace=taxa)
    taxon_state_sets_map = data_mega.taxon_state_sets_map(gaps_as_missing=True)

    tree = dendropy.Tree.get_from_path(tree_loc,
                                       schema="newick", preserve_underscores=True, taxon_namespace=taxa)

    os.remove("./fasta_data.tsv")
    os.remove("./fasta_data.fasta")

    score = fitch_down_pass(postorder_nodes=tree.postorder_node_iter(),
                            taxon_state_sets_map=taxon_state_sets_map)
    fitch_up_pass(tree.postorder_node_iter())

    ###########################################################################
    ## Ok so now the tree has been reconstructed with the fitch up pass down ##
    ## pass algorithm, we'll return the node labelled tree    #################
    ###########################################################################


    return tree

def branch_mutations(tree, hit_csv, embl_csv, embl_reccy):
    total_isolates = []

    for taxon in tree.leaf_node_iter():
        leaf_label = taxon.taxon

        # leaf_label = re.sub("\"","", leaf_label)
        total_isolates.append(leaf_label)

    start_id = []
    isolates_ids = []
    start_node = []
    finish_node = []
    tag = []
    AT = []
    AC = []
    AG = []
    TA = []
    TC = []
    TG = []
    CA = []
    CT = []
    CG = []
    GA = []
    GT = []
    GC = []
    AT_outside = []
    AC_outside = []
    AG_outside = []
    TA_outside = []
    TC_outside = []
    TG_outside = []
    CA_outside = []
    CT_outside = []
    CG_outside = []
    GA_outside = []
    GT_outside = []
    GC_outside = []
    cluster_num = []
    node_chain_len = []
    end_nodes = []
    particular_end_node = []
    node_already_tested = []
    edge_lengths = []


    for isolate in range(len(total_isolates)):
        ## This loops through all the isolates to get all the branches
        ##

        current_row = total_isolates[isolate]
        isolate_test = str(current_row)
        isolate_test = re.sub("\'", "", isolate_test)

        nodes_in_insertion = []
        edge_length = []

        ## First part of the loop establishes the chain of nodes leading to the insertion
        ## node for the element.

        isolate_id = str(current_row)

        current_row = re.sub("\'", "", isolate_id)
        nodes_in_insertion.append(isolate_test)
        tree_node = tree.find_node_with_taxon_label(label=current_row)
        # print(tree_node)
        edge_length.append(tree_node.edge.length)
        parent_node = tree_node.parent_node


        while parent_node != None:
            new_parent = parent_node.parent_node
            old_parent = parent_node
            parent_node = new_parent

            nodes_in_insertion.append(old_parent.label)
            edge_length.append(old_parent.edge.length)

            if parent_node == None:
                end_node = old_parent.label
                if end_node not in end_nodes:
                    end_nodes.append(end_node)

        for nodeys in range(len(nodes_in_insertion) - 1):
            current_end_id = nodes_in_insertion[nodeys]
            if current_end_id not in node_already_tested:
                source_node = nodes_in_insertion[nodeys + 1]
                target_node = current_end_id
                start_node.append(source_node)
                finish_node.append(target_node)
                tag.append("No")

                subset_embl_csv = embl_csv[embl_csv['end_node'] == current_end_id]

                A_T = 0
                A_C = 0
                A_G = 0
                T_A = 0
                T_C = 0
                T_G = 0
                C_A = 0
                C_T = 0
                C_G = 0
                G_A = 0
                G_T = 0
                G_C = 0
                A_T_outside = 0
                A_C_outside = 0
                A_G_outside = 0
                T_A_outside = 0
                T_C_outside = 0
                T_G_outside = 0
                C_A_outside = 0
                C_T_outside = 0
                C_G_outside = 0
                G_A_outside = 0
                G_T_outside = 0
                G_C_outside = 0

                if subset_embl_csv.empty:
                    A_T += 0
                    A_C += 0
                    A_G += 0
                    T_A += 0
                    T_C += 0
                    T_G += 0
                    C_A += 0
                    C_T += 0
                    C_G += 0
                    G_A += 0
                    G_T += 0
                    G_C += 0
                    A_T_outside += 0
                    A_C_outside += 0
                    A_G_outside += 0
                    T_A_outside += 0
                    T_C_outside += 0
                    T_G_outside += 0
                    C_A_outside += 0
                    C_T_outside += 0
                    C_G_outside += 0
                    G_A_outside += 0
                    G_T_outside += 0
                    G_C_outside += 0



                else:

                    A_to_T = len(subset_embl_csv[(subset_embl_csv['start_base'] == "A") & (
                                subset_embl_csv['end_base'] == "T")].index)
                    A_to_G = len(subset_embl_csv[(subset_embl_csv['start_base'] == "A") & (
                                subset_embl_csv['end_base'] == "G")].index)
                    A_to_C = len(subset_embl_csv[(subset_embl_csv['start_base'] == "A") & (
                                subset_embl_csv['end_base'] == "C")].index)
                    T_to_A = len(subset_embl_csv[(subset_embl_csv['start_base'] == "T") & (
                                subset_embl_csv['end_base'] == "A")].index)
                    T_to_G = len(subset_embl_csv[(subset_embl_csv['start_base'] == "T") & (
                                subset_embl_csv['end_base'] == "G")].index)
                    T_to_C = len(subset_embl_csv[(subset_embl_csv['start_base'] == "T") & (
                                subset_embl_csv['end_base'] == "C")].index)
                    G_to_A = len(subset_embl_csv[(subset_embl_csv['start_base'] == "G") & (
                                subset_embl_csv['end_base'] == "A")].index)
                    G_to_T = len(subset_embl_csv[(subset_embl_csv['start_base'] == "G") & (
                                subset_embl_csv['end_base'] == "T")].index)
                    G_to_C = len(subset_embl_csv[(subset_embl_csv['start_base'] == "G") & (
                                subset_embl_csv['end_base'] == "C")].index)
                    C_to_A = len(subset_embl_csv[(subset_embl_csv['start_base'] == "C") & (
                                subset_embl_csv['end_base'] == "A")].index)
                    C_to_T = len(subset_embl_csv[(subset_embl_csv['start_base'] == "C") & (
                                subset_embl_csv['end_base'] == "T")].index)
                    C_to_G = len(subset_embl_csv[(subset_embl_csv['start_base'] == "C") & (
                                subset_embl_csv['end_base'] == "G")].index)

                    A_T += (A_to_T)
                    A_C += (A_to_C)
                    A_G += (A_to_G)
                    T_A += (T_to_A)
                    T_C += (T_to_C)
                    T_G += (T_to_G)
                    C_A += (C_to_A)
                    C_T += (C_to_T)
                    C_G += (C_to_G)
                    G_A += (G_to_A)
                    G_T += (G_to_T)
                    G_C += (G_to_C)

                    subset_reccy_csv = embl_reccy[embl_reccy['end_node'] == current_end_id]

                    if subset_reccy_csv.empty:
                        A_T_outside += A_to_T
                        A_C_outside += A_to_C
                        A_G_outside += A_to_G
                        T_A_outside += T_to_A
                        T_C_outside += T_to_C
                        T_G_outside += T_to_G
                        C_A_outside += C_to_A
                        C_T_outside += C_to_T
                        C_G_outside += C_to_G
                        G_A_outside += G_to_A
                        G_T_outside += G_to_T
                        G_C_outside += G_to_C
                    else:
                        for reccy_row in range(len(subset_reccy_csv.index)):
                            current_start_reccy = subset_reccy_csv.iloc[reccy_row, 2]
                            current_end_reccy = subset_reccy_csv.iloc[reccy_row, 3]
                            for mut_row in range(len(subset_embl_csv.index)):
                                base_pos = subset_embl_csv.iloc[mut_row, 4]
                                start_base = subset_embl_csv.iloc[mut_row, 2]
                                end_base = subset_embl_csv.iloc[mut_row, 3]
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "A" and end_base == "T":
                                    A_to_T -= 1
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "A" and end_base == "C":
                                    A_to_C -= 1
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "A" and end_base == "G":
                                    A_to_G -= 1
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "T" and end_base == "A":
                                    T_to_A -= 1
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "T" and end_base == "C":
                                    T_to_C -= 1
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "T" and end_base == "G":
                                    T_to_G -= 1
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "C" and end_base == "A":
                                    C_to_A -= 1
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "C" and end_base == "T":
                                    C_to_T -= 1
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "C" and end_base == "G":
                                    C_to_G -= 1
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "G" and end_base == "A":
                                    G_to_A -= 1
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "G" and end_base == "T":
                                    G_to_T -= 1
                                if base_pos >= current_start_reccy and base_pos <= current_end_reccy and start_base == "G" and end_base == "C":
                                    G_to_C -= 1

                        A_T_outside += A_to_T
                        A_C_outside += A_to_C
                        A_G_outside += A_to_G
                        T_A_outside += T_to_A
                        T_C_outside += T_to_C
                        T_G_outside += T_to_G
                        C_A_outside += C_to_A
                        C_T_outside += C_to_T
                        C_G_outside += C_to_G
                        G_A_outside += G_to_A
                        G_T_outside += G_to_T
                        G_C_outside += G_to_C

                AT.append(A_T)
                AC.append(A_C)
                AG.append(A_G)
                TA.append(T_A)
                TC.append(T_C)
                TG.append(T_G)
                CA.append(C_A)
                CT.append(C_T)
                CG.append(C_G)
                GA.append(G_A)
                GT.append(G_T)
                GC.append(G_C)
                AT_outside.append(A_T_outside)
                AC_outside.append(A_C_outside)
                AG_outside.append(A_G_outside)
                TA_outside.append(T_A_outside)
                TC_outside.append(T_C_outside)
                TG_outside.append(T_G_outside)
                CA_outside.append(C_A_outside)

                CT_outside.append(C_T_outside)
                CG_outside.append(C_G_outside)
                GA_outside.append(G_A_outside)
                GT_outside.append(G_T_outside)
                GC_outside.append(G_C_outside)
                edge_lengths.append(edge_length[nodeys])
                node_already_tested.append(current_end_id)
                start_id.append(current_row)

        if isolate % 10 == 0:
            print(isolate / len(total_isolates) * 100)

    non_mge_mutations_out = pandas.DataFrame({'start_node': start_node,
                                              'end_node': finish_node,
                                              'tag': tag,
                                              'A-T': AT,
                                              'A-C': AC,
                                              'A-G': AG,
                                              'T-A': TA,
                                              'T-C': TC,
                                              'T-G': TG,
                                              'C-A': CA,
                                              'C-T': CT,
                                              'C-G': CG,
                                              'G-A': GA,
                                              'G-T': GT,
                                              'G-C': GC,
                                              'A-T_clonal': AT_outside,
                                              'A-C_clonal': AC_outside,
                                              'A-G_clonal': AG_outside,
                                              'T-A_clonal': TA_outside,
                                              'T-C_clonal': TC_outside,
                                              'T-G_clonal': TG_outside,
                                              'C-A_clonal': CA_outside,
                                              'C-T_clonal': CT_outside,
                                              'C-G_clonal': CG_outside,
                                              'G-A_clonal': GA_outside,
                                              'G-T_clonal': GT_outside,
                                              'G-C_clonal': GC_outside,
                                              'branch_lengths': edge_lengths,
                                              'starting_isolate': start_id})

    return non_mge_mutations_out


if __name__ == '__main__':

    start_overall = time.perf_counter()

    input_args = get_options()
    ###############################################################################
    ## Load up the csv and then run through each of the clusters to check through #
    ## the res ####################################################################
    ###############################################################################

    cluster_csv = pandas.read_csv(input_args.hit_locs_csv)

    base_loc = input_args.gubbins_res
    unique_clusters = cluster_csv['cluster_name'].unique()

    tot_reccy_csv = pandas.DataFrame()
    tot_non_reccy = pandas.DataFrame()

    seq_clus = 1
    for cluster in unique_clusters:
        print("On cluster: %s, %s of %s" % (cluster, seq_clus, len(unique_clusters)))
        tic_cluster = time.perf_counter()
        current_dat = cluster_csv[cluster_csv['cluster_name'] == cluster]
        current_dir = base_loc + cluster
        cluster_files = os.listdir(current_dir)
        tree_indexio = [k for k, s in enumerate(cluster_files) if "node_labelled.final_tree.tre" in s]
        embl_branch = [k for k, s in enumerate(cluster_files) if "_branch_base.csv" in s]
        embl_reccy = [k for k, s in enumerate(cluster_files) if "_recombinations.csv" in s]

        tree_loc = current_dir + "/" + cluster_files[tree_indexio[0]]
        embl_branch_loc = current_dir + "/" + cluster_files[embl_branch[0]]
        embl_rec_loc = current_dir + "/" + cluster_files[embl_reccy[0]]

        embl_csv = pandas.read_csv(embl_branch_loc)
        embl_reccy_csv = pandas.read_csv(embl_rec_loc)

        ## So now we've got all the files we need for this particular cluster, we'll run through
        ## the tree to get the nodes labelled with the inserts. Then we'll run through the
        ## branches and get the summary of the mutations present

        tree = node_reconstruct(tree_loc, current_dat)
        branches_csv = branch_mutations(tree,current_dat,embl_csv, embl_reccy_csv)

        branches_csv['cluster_name'] = cluster


        branches_out_loc = current_dir + "/" + cluster + "_per_branch_mutations.csv"
        branches_csv.to_csv(path_or_buf=branches_out_loc, index=False)
        seq_clus += 1

    toc_tot = time.perf_counter()

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Mutation finder took: %s (seconds)" % (toc_tot - start_overall))
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Branch mutations found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")




