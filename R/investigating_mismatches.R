###############################################################################
## Checking the new merging protocol versus the old one #######################
###############################################################################

require(dplyr)

old_missing_csv <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_pub_run/pmen_tn916_missing_df.csv",
                            stringsAsFactors = FALSE)
new_missing_csv_pmen_only <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_merging_mk2/pmen_tn916_merging_missing_df.csv",
                            stringsAsFactors = FALSE)
new_missing_csv_gps <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_merging_mk3/pmen_tn916_merging_missing_df.csv",
                                      stringsAsFactors = FALSE)
new_missing_csv_gps_2 <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_merging_mk3/pmen_tn916_merging_missing_df.csv",
                                  stringsAsFactors = FALSE)

library_df_pmen <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_merging_mk2/pmen_tn916_merging_library.csv",
                       stringsAsFactors = FALSE)

library_df_cdc <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_merging_mk3/pmen_tn916_merging_library.csv",
                            stringsAsFactors = FALSE)



## Are these the same ids?

identical(sort(old_missing_csv$id), sort(new_missing_csv_pmen_only$id))

## Now they're not  what's missing?

length(old_missing_csv$id[which(old_missing_csv$id %in% new_missing_csv$id)])
old_missing_csv$id[which(!(old_missing_csv$id %in% new_missing_csv$id))]
new_missing_csv_pmen_only$id[which(!(new_missing_csv_pmen_only$id %in% new_missing_csv_gps$id))]

## So now it looks like adding in the GPS isolates has only added 2 further isolates into the 
## hits department. Lets look at the combination of genes the pipeline is missing. 

dplyr::count(old_missing_csv, reason)
dplyr::count(new_missing_csv_gps, reason)


length_misses <- new_missing_csv_gps[new_missing_csv_gps$reason == "No gene name and length matches",]
dplyr::count(length_misses, merged)

gps_missing <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/tn916_res/gps_tn916_run_free/gps_tn916_missing_df.csv",
                        stringsAsFactors = FALSE)


## Ok, so it looks like all the gains from the merged are now 
## shifted into the no gene name matches. Lets look at these 
## missing gene names 

old_merged <- old_missing_csv[old_missing_csv$reason == "MERGED No good hits before and after",]
gene_no_match <- new_missing_csv_gps[new_missing_csv_gps$reason == "No gene name matches" & 
                                   new_missing_csv_gps$id %in% old_merged$id,]
dplyr::count(gene_no_match, before_gene_name)
dplyr::count(gene_no_match, after_gene_name)

## Combo of genes 
gene_no_match <- gene_no_match %>% mutate(combo = paste(before_gene_name, after_gene_name, sep = "-"))
dplyr::count(gene_no_match, combo)
metpers <- gene_no_match[gene_no_match$combo == "metP-NONE",]
old_missing_csv[old_missing_csv$id %in% metpers$id,c("mge_length","reason")]  








