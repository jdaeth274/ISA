import re
import pandas
import numpy
import sys
import argparse


def get_options():
    purpose = '''This is a python script to intake a csv of an already created library and 
     combine it with a newly created library for wider searching  
    Usage: python lib_combiner.py <new_lib> <old_lib>  <out_csv>'''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='pen_checker_cdc.py')

    parser.add_argument('--new_lib', required=True, help='Out merged BLAST csv (required)', type=str)
    parser.add_argument('--old_lib', required=True, help='Isolates gff and references gff csv (required)', type=str)
    parser.add_argument('--out_csv', required=True, help='Out Library csv name (required)', type=str)

    args = parser.parse_args()

    return args


def library_integrator(library_csv, prospective_csv, isolate_id):
    ## Function to decide whether to merge a hits data into the library csv.
    ## Input: library_csv: The current library csv to search for matches against
    ##        prospective_csv: The prospective set of results to be potentially merged
    ## This will first work on the basis of hits sharing almost identical blast matches
    ## Then we'll look into flank region composition and finally the total insert composition
    ## to see if this is a novel hit or not.

    lib_new = library_csv.copy()
    ids_to_drop = []

    ## lets narrow down first by number of genes in the element if this is not exact I think its
    ## fair to assume this would be novel.
    #'before_flank_gene', 'after_flank_gene', 'before_flank_avg',
    #'after_flank_avg'




    ## So there is a hit with the same number of mge genes, let now check the element length +- 2 bp
    ## Needs to be a hit with no length +- 2bp and no genes +- 1

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
                        novel_hit = True
                    else:
                        remain_48 = pandas.concat([gene_hits, length_hits], ignore_index=True, sort = False)
                        remain_48 = remain_48.drop_duplicates()
                        before_gene_name = prospective_csv['before_gene_name'][0]
                        after_gene_name = prospective_csv['after_gene_name'][0]

                        gene_name_matches = remain_48[(remain_48['before_gene_name'] == before_gene_name) &\
                                                      (remain_48['after_gene_name'] == after_gene_name)]
                        if gene_name_matches.empty:
                            novel_hit = True
                        else:

                            flanks_length = remain_48[remain_48['flanks_length'] > prospective_csv['flanks_length'][0]]


                            if flanks_length.empty:
                                ids_to_drop = remain_48['id'].tolist()

                                novel_hit = True
                            else:
                                novel_hit = False


                else:
                    novel_hit = True
            else:
                novel_hit = True
    else:
        novel_hit = True

    if len(ids_to_drop) > 0:
        indies_to_drop = lib_new[lib_new['id'].isin(ids_to_drop)].index
        lib_new = lib_new.drop(indies_to_drop)
        lib_new = lib_new.reset_index(drop=True)
    if novel_hit:
        lib_new = lib_new.append(prospective_csv, ignore_index = True, sort = False)
        lib_new = lib_new.reset_index(drop=True)




    return(lib_new)

class library_df():
    def __init__(self, df_loc):
        self.library = pandas.read_csv(df_loc)

    @property
    def nrow(self):
        return (len(self.library.index))


if __name__ == '__main__':

    ## get options
    input_args = get_options()

    old_lib = library_df(input_args.old_lib)
    new_lib = library_df(input_args.new_lib)
    print("New library of length %s" % new_lib.nrow)
    print("Old library of length %s" % old_lib.nrow)
    start_length = new_lib.nrow

    for index, row in old_lib.library.iterrows():
        df_row = pandas.DataFrame(row).transpose().reset_index(drop=True)
        new_lib.library = library_integrator(new_lib.library, df_row, row['id'])

    added_rows = new_lib.nrow - start_length

    print("Added in %s row(s), library is now has %s row(s)" % (added_rows, new_lib.nrow))

    new_lib.library.to_csv(path_or_buf=input_args.out_csv, index=False)




