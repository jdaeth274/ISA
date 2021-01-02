###############################################################################
## To be run in /rds/general/project/bacterial_evo_genomics/live/gps_annotations_4_2_2020/gps_pbp_runs

require(dplyr)

missing_files <- list.files("./", "__missing_isolates.txt")

missing_iso <- NULL

for(k in missing_files){
  current_files <- readLines(k)
  
  current_files <- sub(".velvet.gff", "",current_files)
  
  missing_isos <- unique(current_files)
  current_cluster <- sub("__missing_isolates.txt","",sub("cluster_","gpsc.",k))
  
  current_df <- cbind.data.frame(missing_isos, rep(current_cluster))
  colnames(current_df) <- c("id", "cluster_name")
  
  missing_iso <- bind_rows(missing_iso, current_df)
  
}


hit_files <- list.files("./", "__pbp_predictions.csv")

gps_hits <- NULL

for(k in hit_files){
  current_hits <- read.csv(k, stringsAsFactors = FALSE) %>% select(c(id, penicillin_cat))
  
  current_cluster <- sub("__pbp_predictions.csv","",sub("cluster_","gpsc.",k))
  
  current_hits$cluster_name <- current_cluster
  
  gps_hits <- bind_rows(gps_hits, current_hits)
  
}


write.csv(missing_iso, file = "./gps_missing_pbp_isolates.csv", row.names = FALSE)

write.csv(gps_hits, "./gps_pbp_profiles.csv", row.names = FALSE)



