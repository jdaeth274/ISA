###############################################################################
## Script to check the missing isolates based on gubbins trees has to be run ##
##  in directory with tree and fasta list files ###############################
###############################################################################

require(ape, quietly = TRUE)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

out_file <- commandArgs(trailingOnly = TRUE)[1]

run_files <- list.files(path = ".", full.names = TRUE)

tree_file <- run_files[grep(".node_labelled.final_tree.tre",run_files)]
fasta_list_file <- run_files[grep("_fasta_list.txt",run_files)]

cluster_tree <- read.tree(tree_file)
cluster_ids <- sub("\\.[a-z,A-Z].*$","",basename(readLines(fasta_list_file)))

missing_from_tree <- cluster_tree$tip.label[which(!(cluster_tree$tip.label %in% cluster_ids))]

cluster_name <- sub("_.*$","",fasta_list_file)
cluster_name <- sub("\\./","",cluster_name)

if(length(missing_from_tree) > 0){
  
  
  missing_df <- cbind.data.frame(missing_from_tree, rep(cluster_name, length(missing_from_tree)))
  write.table(missing_df, file = out_file, 
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  
}
cat(paste("Finished cluster:", cluster_name))