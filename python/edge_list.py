import pandas
import numpy
import argparse

def get_options():
    purpose = '''This script intakes the hit allocator csv of hit locations for an MGE and then returns an edge list 
        of the hits in order to form a cytoscape graph 
        Usage: python edge_list.py --hit_csv <hit_csv> --out_name <out_csv>'''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='pen_checker_cdc.py')

    parser.add_argument('--hit_csv', required=True, help='Hit csv from hit_allocator script"', type=str)
    parser.add_argument('--out_name', required=True, help='Out edge listname (required)', type=str)

    args = parser.parse_args()

    return args

def col_creator(hit_csv):
    ## Function to create a df with the columns representing the insertion points of an element.
    ## Input: hit_csv: The hit csv output from hit_allocator script
    colnames_hits = hit_csv['insert_name'].unique().tolist()

    cols_df = numpy.zeros(shape=(len(hit_csv.index), len(colnames_hits)))
    cols_df = pandas.DataFrame(cols_df)
    hit_col_num = 0

    for hit in colnames_hits:
        hit_isos = hit_csv[hit_csv['insert_name'] == hit]
        hit_ids = hit_isos['id']
        hit_col = hit_ids.append(pandas.Series(numpy.repeat(numpy.nan, (len(hit_csv.index) - len(hit_isos.index)))), ignore_index=True)
        cols_df.iloc[:,hit_col_num] = pandas.Series(hit_col, index= cols_df.index)
        hit_col_num += 1

    return cols_df

if __name__ == '__main__':
    input_args = get_options()

    results_numpy = numpy.zeros(shape=(1, 2))
    resultslist = pandas.DataFrame(data=results_numpy, columns=['source', 'target'])

    hit_csv = pandas.read_csv(input_args.hit_csv)

    cluster_cols = col_creator(hit_csv)


    for busta in range(len(cluster_cols.columns)):
        current_col = cluster_cols.iloc[:, busta]
        current_col = current_col.dropna()
        current_col = list(current_col)
        if len(current_col) > 1:
            for person in current_col:
                myindex = current_col.index(person)
                newlist = current_col[:myindex] + current_col[myindex + 1:]  # make a new temp list without the person in it
                for item in newlist:
                    mytuple = pandas.Series(data=[person, item], index=['source', 'target'])
                    backtuple = pandas.Series(data=[item, person])
                    back_test = (resultslist['source'].isin([item]).any()) & (
                        resultslist['target'].isin([person]).any())
                    if back_test == False:  # remove any reversed duplicates
                        resultslist = resultslist.append(mytuple, ignore_index=True)
        else:
            mytuple = pandas.Series(data=[current_col[0], numpy.nan], index=['source', 'target'])
            resultslist = resultslist.append(mytuple, ignore_index=True)

    resultslist = resultslist.iloc[1:]

    resultslist.to_csv(input_args.out_name, index=False, sep="\t")