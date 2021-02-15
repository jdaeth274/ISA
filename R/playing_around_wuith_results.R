###############################################################################
## Reading in the hit df for gps for mega #####################################

require(ggtree)

gps_mega_hits <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/gps_mega_hits_df.csv",
                          stringsAsFactors = FALSE)


## look at the gpsc 10 res

gps_mega_10_hits <- gps_mega_hits

gps_tn916_missing <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/tn916_res/gps_tn916_run_free/gps_tn916_missing_df.csv",
                              stringsAsFactors = FALSE)
gps_tn916_assigned <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/tn916_res/gps_tn916_run_free/gps_tn916_hits_df.csv",
                               stringsAsFactors = FALSE)

gps_tn916_reccy <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/tn916_res/gps_tn916_run_free/gps_tn916_reccy_hits.csv",
                            stringsAsFactors = FALSE) %>% mutate(reccy = "Yes") %>% rename(isolate_id = isolate_example)
gps_tn916_reccy_miss <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/tn916_res/gps_tn916_run_free/gps_tn916_non_reccy_hits.csv",
                            stringsAsFactors = FALSE) %>% mutate(reccy = "No")

gps_tn916_assigned <- gps_tn916_assigned %>% mutate(insert_node = paste(cluster_name, insertion_node, sep = "-")) 

gps_tn916_assigned %>% count(insert_node) %>% nrow()  
gps_tn916_assigned %>% count(cluster_name) %>% nrow()


merged_916_insertions <- bind_rows(gps_tn916_reccy, gps_tn916_reccy_miss)
merged_916_insertions %>% count(cluster_name) %>% nrow()
merged_916_insertions %>% count(insert_name) %>% arrange(desc(n)) %>% head()
merged_916_insertions %>% count(reccy)
merged_916_insertions[which.max(merged_916_insertions$insertion_node_tip_nums),]
count(gps_tn916_missing, reason)

count(gps_tn916_assigned, insert_name) %>% arrange(desc(n)) %>% head()
gps_tn916_assigned %>% filter(cluster_name == "gpsc.1") %>% count(insert_name) %>% arrange(desc(n)) %>% head()
gps_tn916_assigned %>% filter(insert_name == 106) %>% count(cluster_name) %>% arrange(desc(n)) %>% head()

gps_tn916_library <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/tn916_res/gps_tn916_run_free/gps_tn916_library.csv",
                              stringsAsFactors = FALSE)

gps_tn916_assigned[gps_tn916_assigned$insert_name == "121",]

ggplot(data = gps_tn916_assigned) +
  geom_histogram(aes(x = insert_length), binwidth = 5000) 
  
  


###############################
## gps mega summaries 

gps_mega_hits <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/mega_res/gps_mega_run4/gps_mega_hits_df.csv", 
                          stringsAsFactors = FALSE)
gps_mega_misses <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/mega_res/gps_mega_run4/gps_mega_missing_df.csv",
                            stringsAsFactors = FALSE)
gps_mega_reccy_hits <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/mega_res/gps_mega_run4/gps_mega_reccy_hits.csv",
                                stringsAsFactors = FALSE) %>% mutate(reccy = "Yes") %>% rename(isolate_id = isolate_example)
gps_mega_reccy_misses <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/mega_res/gps_mega_run4/gps_mega_non_reccy_hits.csv",
                                stringsAsFactors = FALSE) %>% mutate(reccy = "No")


merged_insertions <- bind_rows(gps_mega_reccy_hits, gps_mega_reccy_misses)
count(merged_insertions, insert_name) %>% arrange(desc(n)) %>% head()
count(merged_insertions, reccy)

count(gps_mega_misses, reason)
count(gps_mega_hits, insert_name) %>% arrange(desc(n)) %>% head()

gps_mega_hits <- gps_mega_hits %>% mutate(insert_node = paste(cluster_name, insertion_node, sep = "-"))

colnames(gps_mega_reccy_hits)
colnames(gps_mega_reccy_misses)

count(gps_mega_hits, insert_node) %>% nrow()

gps_tag <- gps_mega_hits %>% filter(insert_name == 28) %>% count(cluster_name) %>% nrow()

gps_mega_hits <- mutate(gps_mega_hits, assigned = "Yes")
gps_mega_misses <- mutate(gps_mega_misses, assigned = "No")

gps_mega_tot <- bind_rows(gps_mega_hits, gps_mega_misses)

gps_tag <- gps_mega_hits %>% filter(insert_name == 25) %>% count(cluster_name) %>% nrow()


################################
## tn916 summaries 

gps_tn916_assigned <- mutate(gps_tn916_assigned, assigned = "Yes")
gps_tn916_missing <- mutate(gps_tn916_missing, assigned = "No")

gps_tn916_tot <- bind_rows(gps_tn916_assigned, gps_tn916_missing)



gps_mega_gpsc1 <- gps_mega_hits[gps_mega_hits$cluster_name == "gpsc.1",]

gpsc_1_tree <- read.tree("~/Dropbox/phd/insertion_site_analysis/data/gpsc.1_run_data/gpsc.1.node_labelled.final_tree.tre")

gpsc_1_csv <- cbind.data.frame(gpsc_1_tree$tip.label, rep("ND", length(gpsc_1_tree$tip.label)), rep("ND",length(gpsc_1_tree$tip.label)))
colnames(gpsc_1_csv) <- c("id","mega","tn916")

gps_mega <- gps_mega_hits[,c("id","insert_name")] %>% rename(mega = insert_name)
gps_916 <- gps_tn916_assigned[,c("id","insert_name")] %>% rename(tn916 = insert_name)

gpsc_1_csv <- left_join(gpsc_1_csv, gps_mega, by = "id") %>% mutate(mega.y = ifelse(is.na(mega.y),mega.x,mega.y)) %>% select(-mega.x) %>%
  rename(MEGA__autocolour = mega.y) %>% left_join(gps_916, by = "id") %>% mutate(tn916.y = ifelse(is.na(tn916.y),tn916.x,tn916.y)) %>%
  rename(Tn916__autocolour = tn916.y) %>% select(-tn916.x) %>% mutate(Tn916 = ifelse(Tn916__autocolour == "179", "Tn2010", Tn916__autocolour)) %>%
  mutate(Tn916 = ifelse(Tn916__autocolour == "205", "Tn2009", Tn916)) %>% mutate(Tn916 = ifelse(Tn916 %in% c("Tn2010","Tn2009","ND"),Tn916,"Other Tn916")) %>%
  mutate(Tn916 = ifelse(Tn916 == "ND", "Absent",Tn916)) %>% mutate(Tn1207.1 = ifelse(MEGA__autocolour == "ND", "Absent","Present")) %>% 
  mutate(bar_val = 2)
  
head(gpsc_1_csv)
gpsc_1_tree <- read.tree("~/Dropbox/phd/insertion_site_analysis/data/gpsc.1_run_data/gpsc.1.node_labelled.final_tree.tre")
clade_groupings <- c(Tn2010 = 922, Tn2009 = 907)
tree <- groupClade(gpsc_1_tree, clade_groupings)
cols <- c(Tn2010 = "red", Tn2009 = "blue")
p <- ggtree(tree, aes(color = group))

p <- ggtree(tree, aes(color = group), ladderize = FALSE) %>% rotate(rootnode(tree)) + 
  geom_treescale(x = 0, y = 1, width = 0.002) + 
  scale_color_manual(values = c(cols, "black"), na.value = "black", name = "Lineage",
                     breaks = c("Tn2010", "Tn2009")) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15)))
  
p
p <- p +scale_color_manual(values = c(cols, "black"), na.value = "black", name = "Lineage",
                      breaks = c("Tn2010", "Tn2009")) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
  theme_tree2(legend.position = c(.1, .88))
p
## visualize the tree with tip labels and tree scale
p <- ggtree(tree, aes(color = group), ladderize = FALSE) %>% scale_color_manual(values = c(cols, "black"), na.value = "black", name = "Lineage",
                     breaks = c("Tn2010", "Tn2009")) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
  theme_tree2(legend.position = c(.1, .88))

bar_data_gpsc1 <- gpsc_1_csv[,c("id","bar_val")]

gpsc_1_tree_plot <- ggtree(gpsc_1_tree, ladderize = FALSE)
gpsc_1_tree_plot <- gpsc_1_tree_plot %<+% gpsc_1_csv + geom_tippoint(aes(colour=Tn916))
gpsc_1_tree_plot <- gpsc_1_tree_plot + geom_facet(panel = "Mega", data = bar_data_gpsc1, geom = ggstance::geom_barh, 
                                                  aes(x = bar_val, color = Tn1207.1, fill = Tn1207.1), 
                                                  stat = "identity", width = .6) 

gpsc_1_tree_plot
write.csv(gpsc_1_csv, file = "~/Dropbox/phd/insertion_site_analysis/data/gpsc.1_run_data/gpsc1_micro_csv.csv",
          row.names = FALSE, quote = FALSE)


pmen9_epi_csv <- left_join(x = pmen9_epi_csv,y = pmen_mega_narrow, by = "id") %>% 
  mutate(insert_name_mega = ifelse(is.na(insert_name_mega), 0, insert_name_mega)) %>% rename(insert_name_mega__autocolour = insert_name_mega) %>% left_join(y = pmen_tn916_narrow, by = "id") %>% 
  mutate(insert_name_tn916 = ifelse(is.na(insert_name_tn916), 0, insert_name_tn916)) %>% rename(insert_name_tn916__autocolour = insert_name_tn916) %>%
  left_join(pmen_sulfa, by = "id") %>% mutate(sulfa__autocolour = ifelse(is.na(sulfa__autocolour), "ND", sulfa__autocolour))



################################################################################
## Merging the hits dfs ########################################################
################################################################################

gps_tn916_tot <- gps_tn916_tot %>% mutate(Tn916 = "Tn916")
gps_mega_tot <- gps_mega_tot %>% mutate(Mega = "Mega")

binded_rows <- bind_rows(gps_tn916_tot, gps_mega_tot)

merged_insertions_tot <- full_join(gps_tn916_tot, gps_mega_tot, by = "id") %>% mutate(Tn916 = ifelse(is.na(Tn916),"",Tn916)) %>%
  mutate(Mega = ifelse(is.na(Mega),"",Mega)) %>% mutate(element = paste(Mega, Tn916, sep = "#")) %>%
  mutate(Cluster_name = ifelse(is.na(cluster_name.x),cluster_name.y, cluster_name.x))

count(merged_insertions_tot, id) %>% nrow()

count(merged_insertions_tot, element)



###############################################################################
## Cluster prevalences ########################################################
###############################################################################

cluster_list <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/gps_reference_isolate_gff.csv",
                         stringsAsFactors = FALSE)

cluster_sizes <- count(cluster_list, cluster_name) %>% rename(size = n)

## Mega prevalences 

mega_cluster_nums <- count(gps_mega_tot, cluster_name) %>% rename(mega = n)

mega_bound <- left_join(mega_cluster_nums, cluster_sizes, by = "cluster_name") %>% mutate(prev = mega/size)

mean(mega_bound$prev)

## tn916 prevalences 

tn916_cluster_nums <- count(gps_tn916_tot, cluster_name) %>% rename(Tn916 = n)

tn916_bound <- left_join(tn916_cluster_nums, cluster_sizes, by = "cluster_name") %>% mutate(prev = Tn916/size)

mean(tn916_bound$prev)


dual_inserts <- merged_insertions_tot %>% group_by(id) %>% mutate(element = paste(Mega, Tn916, sep = "#"))


cluster_fasta_gps <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/gps_reference_isolate_gff.csv",
                              stringsAsFactors = FALSE)

cluster_fasta_gps$one <- 1
cluster_fasta_nums <- aggregate(one ~ cluster_name, cluster_fasta_gps, FUN = sum)

gps_tn916_count <- count(gps_tn916_tot, cluster_name)

tn916_prev <- NULL
for(cluster in gps_tn916_count$cluster_name){
  prevalence <- gps_tn916_count[gps_tn916_count$cluster_name == cluster,2] /
    cluster_fasta_nums[cluster_fasta_nums$cluster_name == cluster,2]
  tn916_prev <- append(tn916_prev, prevalence)
  
}

gps_mega_count <- count(gps_mega_tot, cluster_name)

mega_prev <- NULL
for(cluster in gps_mega_count$cluster_name){
  prevalence <- gps_mega_count[gps_mega_count$cluster_name == cluster,2] /
    cluster_fasta_nums[cluster_fasta_nums$cluster_name == cluster,2]
  mega_prev <- append(mega_prev, prevalence)
  
}

##rlmcd
nrow(gps_mega_hits[gps_mega_hits$insert_name %in% c(7,8,12,18,22,30,32),])





###############################################################################
## Check distribution of pen and MEGA #########################################
###############################################################################

gps_pbp_profs <- read.csv("~/Dropbox/phd/elife_paper/data/gps_pbp_profiles.csv",
                          stringsAsFactors = FALSE)

gps_mega_tot$mega <- "Yes"

gps_mega_binding <- gps_mega_tot[,c("id","mega")]

gps_mega_pen <- left_join(gps_pbp_profs, gps_mega_binding, by = c("id" = "id")) %>%
  mutate(mega = ifelse(is.na(mega), "No","Yes"))

gps_mega_pen_yes <- gps_mega_pen[gps_mega_pen$mega == "Yes",]
















