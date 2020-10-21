import pandas
import re
import subprocess


reference_loc = pandas.read_csv("/rds/general/project/bacterial_evo_genomics/live/gps_annotations_4_2_2020/gps_gubbins_runs/tot_missing_isolates.tsv", sep= "\t")
references_tot = pandas.read_csv("./gps_run_data/gps_reference_isolate_fasta.csv")


act_fasta_list = pandas.DataFrame()
dna_list = []

for k in range(len(reference_loc.index)):
	current_isolate = reference_loc.iloc[k,0]
	current_cluster = reference_loc.iloc[k,1]
	current_pos_fa = "../missing_refs/" + current_isolate + ".contigs_velvet.fa"
	new_pos_fa = current_isolate + ".velvet.fasta"
	new_pos_gff = current_isolate + ".velvet.gff"
	current_pos_gff = "../missing_refs/" + current_isolate + ".velvet.gff"
	cluster_num = re.sub("gpsc\.","",current_cluster)
	mv_fa = "cp " + current_pos_fa + " ../cluster_" + cluster_num + "_list_gffs/" + new_pos_fa
	mv_gff = "cp " + current_pos_fa + " ../cluster_" + cluster_num + "_list_gffs/" + new_pos_gff
	subprocess.run(mv_fa, shell=True)
	subprocess.run(mv_gff, shell=True)
	ref_fasta = references_tot[references_tot['cluster_name'] == current_cluster].iloc[0,1]
	new_fasta_line = pandas.DataFrame()
	new_fasta_loc = "/rds/general/project/bacterial_evo_genomics/live/gps_annotations_4_2_2020/cluster_" + cluster_num + "_list_gffs/" + new_pos_fa
	new_fasta_line['fasta_loc'] = pandas.Series(new_fasta_loc)
	new_fasta_line['refernce'] = pandas.Series(ref_fasta, index=new_fasta_line.index)
	dna_list.append(new_fasta_loc)
	act_fasta_list = act_fasta_list.append(new_fasta_line, ignore_index=True, sort=False)

act_fasta_list.to_csv(path_or_buf="./missing_act_csv.csv", index=False)
with open("./missing_fasta_list.txt", mode='wt', encoding='utf-8') as myfile:
	myfile.write('\n'.join(dna_list))

print("Finished")







