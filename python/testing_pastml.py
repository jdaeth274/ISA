from pastml.acr import pastml_pipeline
import ete3
import time
import dendropy
import pandas

data = "/home/jd2117/Dropbox/phd/insertion_site_analysis/pmen_mega_reccy_test/pmen9_fasta_.csv"

tree = "/home/jd2117/Dropbox/phd/insertion_site_analysis/data/pmen_run/pmen9/PMEN9_gubbins.node_labelled.final_tree.tre"

columns = ['mega']

html_compressed = "/home/jd2117/Dropbox/phd/insertion_site_analysis/pmen9_html.html"


tic_ml = time.perf_counter()
pastml_pipeline(data = data, columns = columns, name_column = 'mega',
                 tree = tree, verbose = True,
                 out_data = "./states_res.txt",
                work_dir = "./")
toc_ml = time.perf_counter()

print("ML recon took this long: %s (seconds)" % (round(toc_ml - tic_ml)))

pmen9_tree = dendropy.Tree.get(path=tree,
                             schema="newick", preserve_underscores=True)

tree_states = pandas.read_csv("./states_res.txt", sep="\t")

for node in pmen9_tree:

    tree_states_node = tree_states[tree_states['node'] == node.label]
    if not tree_states_node.empty:
        print(tree_states_node['mega'].values[0])


