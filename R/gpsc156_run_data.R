###############################################################################
## gpsc156 investigations #####################################################
###############################################################################

require(ape)
require(dplyr)

###############################################################################
## Load up the tree, create the microreact csv, check where tag inserting in ##
## this grouping, why we could be getting two isolates for the same region? ###
###############################################################################

gpsc156_tree <- read.tree("~/Documents/phd_downloads/gps_downloads/gpsc.156_run_data/gpsc.156.node_labelled.final_tree.tre")

hits_df <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/mega_res/gps_mega_run4/gps_mega_hits_df.csv",
                    stringsAsFactors = FALSE) %>% filter(cluster_name == "gpsc.156")



which(gpsc156_tree$node.label == "internal_21")

gpsc_156_csv <- cbind.data.frame(gpsc156_tree$tip.label, rep("No",length(gpsc156_tree$tip.label)))
colnames(gpsc_156_csv) <- c("id","mega__autocolour")

gpsc_156_csv <- gpsc_156_csv %>% mutate(mega__autocolour = ifelse(id %in% hits_df$id, "Yes","No"))

write.csv(gpsc_156_csv,file = "~/Documents/phd_downloads/gps_downloads/gpsc.156_run_data/gpsc156_micro.csv",
          row.names = FALSE)
