###############################################################################
## PBP recombination pmen results #############################################
###############################################################################

require(dplyr)
require(ggpubr)

graph_getter_pbp <- function(dir_to_files,graph_name){
  
  #browser()
  last_character <- base::substr(dir_to_files,
                                 nchar(dir_to_files),
                                 nchar(dir_to_files))
  if(last_character != "/"){
    dir_to_files <- paste(dir_to_files, "/", sep = "")
  }
  summo_csv_file <- list.files(dir_to_files, pattern = "*species_compo*")
  summary_csv <- paste(dir_to_files, summo_csv_file,sep = "")
  
  summo_csv <- read.csv(summary_csv, stringsAsFactors = FALSE)
  blast_results <- list.files(dir_to_files, pattern = "*list*.csv")

  summo_csv <- summo_csv %>% mutate(pre_resistance = ifelse(pre_resistance == "reference","ref",pre_resistance)) %>%
    mutate(resistance = ifelse(resistance == "reference", "ref",resistance))
  summo_csv$insertion_point <- paste(summo_csv$pre_resistance, summo_csv$resistance, summo_csv$gene_name, sep = "-")
  
  counted_summo <- plyr::count(summo_csv$insertion_point)  
  
  #############################################################################
  ## Ok, so now we've got the summary_csv, lets go through each of the ########
  ## insertion sites, getting overall pneumo weight by bitscore, rank of top ##
  ## pneumo hit and then top hit bit / bit of top pneumo hit ##################
  #############################################################################
  
  for(k in 1:nrow(counted_summo)){
    #browser()
    current_insertion <- as.character(counted_summo[k, 1])
    subsetted_data <- summo_csv[summo_csv$insertion_point == current_insertion,]
    
    if(k == 1){
      results_df <- bitscore_res(subsetted_data, blast_results,
                                 current_insertion, dir_to_files)
    }else{
      new_df <- bitscore_res(subsetted_data, blast_results,
                             current_insertion, dir_to_files)
      results_df <- rbind.data.frame(results_df, new_df)
    }
    
  }
  
  total_run_through <- bitscore_res(summo_csv, blast_results,
                                    "total", dir_to_files)
  
  #browser()
  total_hist_weighted <- ggplot(data = total_run_through) + geom_violin(aes(x = insertion_loci, 
                                                                             y = weighted_bit)) + 
    labs(title = paste(graph_name, "Weighted total"))
  total_hist_ranked <- ggplot(data = total_run_through) + geom_violin(aes(x = insertion_loci,
                                                                           y = ranked_pneumo)) +
    labs(title = paste(graph_name, "ranked pneumo"))
  total_bit_ranked <- ggplot(data = total_run_through) + geom_violin(aes(x = insertion_loci,
                                                                          y = bit_pneumo)) +
    labs(title = paste(graph_name, "bitscore pneumo"))
  
  arranged_plots <- ggarrange(total_hist_weighted, total_hist_ranked, total_bit_ranked,
                              ncol = 3, nrow = 1, align = "hv")
  
  
  hist_weighted <- ggplot(data = results_df) + geom_boxplot(aes(x = insertion_loci, 
                                                                y = weighted_bit)) + 
    labs(title = paste(graph_name, "Weighted total"))
  hist_ranked <- ggplot(data = results_df) + geom_boxplot(aes(x = insertion_loci,
                                                              y = ranked_pneumo)) +
    labs(title = paste(graph_name, "ranked pneumo"))
  bit_ranked <- ggplot(data = results_df) + geom_violin(aes(x = insertion_loci,
                                                             y = bit_pneumo)) +
    labs(title = paste(graph_name, "bitscore pneumo")) + geom_point(aes(x = insertion_loci,
                                                                        y = bit_pneumo, color = insertion_loci),
                                                                    position = position_jitter(width = 0.1, height = 0))
  
  
  
  
  
  return(list(by_insertion = results_df, total = total_run_through, total_plots = arranged_plots,
              weighted_hist = hist_weighted, ranked_hist = hist_ranked, ranked_bit = bit_ranked))
}

bitscore_res <- function(data_subset, blast_res_list, insertion_loci, dir_to_files){
  
  #browser()
  isolate_file <- paste(data_subset$isolate,"_",data_subset$gene_name, "_species_list.csv",sep = "")
  
  files_to_open <- blast_res_list[which(blast_res_list %in% isolate_file)]
  files_to_open <- paste(dir_to_files, files_to_open, sep = "")
  
  #############################################################################
  ## So now we've got the vector of files to open in files_to_open, we'll loop#
  ## through these and get the weight by bitscores, rank of top_pneumo, and ###
  ## bitscore of top pneumo/top hit for each of the files #####################
  #############################################################################
  #browser()
  output_data <- data.frame(data = matrix(NA, nrow = length(files_to_open),
                                          ncol = 5))
  colnames(output_data) <- c("isolate","insertion_loci","weighted_bit","ranked_pneumo",
                             "bit_pneumo")
  output_data$isolate <- data_subset$isolate
  output_data$insertion_loci <- rep(insertion_loci, nrow(output_data))
  
  for(k in 1:length(files_to_open)){
    current_file <- read.csv(files_to_open[k])
    current_file_ordered <- current_file[order(current_file$bitscore),]
    
    pneumo_rows <- pneumo_finder(current_file)
    
    
    if(length(pneumo_rows) == 0){
      output_data[k,3:5] <- 0
    }else{
      ## weighted bitscores first ##
      bit_score_tot <- sum(current_file$bitscore)
      pneumo_bits <- sum(current_file$bitscore[pneumo_rows])
      output_data[k, 3] <- pneumo_bits / bit_score_tot
      
      ## ranked pneumo second ##
      
      ranking <- (nrow(current_file) - pneumo_rows[1] + 1) / nrow(current_file)
      output_data[k, 4] <- ranking
      
      ## bit_top_pnuemo ##
      
      rank_pneumo_bit <- current_file$bitscore[pneumo_rows[1]]
      top_bit <- current_file$bitscore[1] 
      output_data[k, 5] <- rank_pneumo_bit / top_bit
      
      
      
      
      
      
    }
    
  }
  
  return(output_data)
  
  
  
}

pneumo_finder <- function(data_set){
  
  #############################################################################
  ## This function finds the position of pneumo isolates in the blast results #
  #############################################################################
  
  pneumo_rows <- NULL
  for(k in 1:nrow(data_set)){
    pneumo_grep <- grepl("strep_pneumo", data_set[k,2])
    gpsc_grep <- grepl("GPSC", data_set[k, 2])
    
    if(pneumo_grep == TRUE | gpsc_grep == TRUE){
      pneumo_rows <- c(pneumo_rows, k)
    }
  }
  
  return(pneumo_rows)
}


pbps_only <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/pbp_res/gps_pbp_extraction/500_pbp_seqs_res/gps_pbp_species_compo_pbp.csv",
                      stringsAsFactors = FALSE)
flanks_extracted <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/pbp_res/gps_pbp_extraction/gps_pbp_500_flanks_extracted.csv",
                             stringsAsFactors = FALSE)

pbp_res <- graph_getter_pbp("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/pbp_res/gps_pbp_extraction/500_pbp_seqs_res/",
                        graph_name = "PBP genes")
pbp_res$ranked_bit

results_df <- pbp_res$by_insertion %>% filter(insertion_loci %in% c("S-R-pbp1A","S-R-pbp2X","S-R-pbp2B")) %>% mutate(Gene = insertion_loci)



bit_ranked <- ggplot(data = results_df) + geom_violin(aes(x = Gene,
                                                          y = bit_pneumo)) + geom_point(aes(x = insertion_loci,
                                                                      y = bit_pneumo, color = Gene),
                                                                  position = position_jitter(width = 0.1, height = 0)) +
  labs(y = "Score") 
  
bit_ranked

## Work out the mean bit_pneumo scores for each of the genes 

bit_df <- results_df[,c(2,5)]

bit_df <- aggregate(bit_pneumo ~ insertion_loci, bit_df, FUN = median)
bit_df

## Plot out the R-S changes too 

results_df <- pbp_res$by_insertion %>% filter(insertion_loci %in% c("R-S-pbp1A","R-S-pbp2X","R-S-pbp2B")) %>% mutate(Gene = insertion_loci)

bit_ranked <- ggplot(data = results_df) + geom_violin(aes(x = Gene,
                                                          y = bit_pneumo)) + geom_point(aes(x = insertion_loci,
                                                                                            y = bit_pneumo, color = Gene),
                                                                                        position = position_jitter(width = 0.1, height = 0)) +
  labs(y = "Score") 

bit_ranked




