###############################################################################
## Outputting flank scores for different cluster scores #######################
###############################################################################

require(ggplot2, quietly = TRUE)
require(stringr, quietly = TRUE)
require(ggpubr, quietly = TRUE)
require(tictoc, quietly = TRUE)

graph_getter <- function(dir_to_files,graph_name){
  
  
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
  
  counted_summo <- plyr::count(summo_csv$insertion_point)  
  
  #############################################################################
  ## Ok, so now we've got the summary_csv, lets go through each of the ########
  ## insertion sites, getting overall pneumo weight by bitscore, rank of top ##
  ## pneumo hit and then top hit bit / bit of top pneumo hit ##################
  #############################################################################
  
  for(k in 1:nrow(counted_summo)){
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
  
  
  total_hist_weighted <- ggplot(data = total_run_through) + geom_boxplot(aes(x = insertion_loci, 
                                                                             y = weighted_bit)) + 
    labs(title = paste(graph_name, "Weighted total"))
  total_hist_ranked <- ggplot(data = total_run_through) + geom_boxplot(aes(x = insertion_loci,
                                                                           y = ranked_pneumo)) +
    labs(title = paste(graph_name, "ranked pneumo"))
  total_bit_ranked <- ggplot(data = total_run_through) + geom_boxplot(aes(x = insertion_loci,
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
  bit_ranked <- ggplot(data = results_df) + geom_boxplot(aes(x = insertion_loci,
                                                             y = bit_pneumo)) +
    labs(title = paste(graph_name, "bitscore pneumo"))
  
  
  
  
  
  return(list(by_insertion = results_df, total = total_run_through, total_plots = arranged_plots,
              weighted_hist = hist_weighted, ranked_hist = hist_ranked, ranked_bit = bit_ranked))
}

bitscore_res <- function(data_subset, blast_res_list, insertion_loci, dir_to_files){
  
  
  isolate_file <- paste(data_subset$isolate, "_species_list.csv",sep = "")
  files_to_open <- blast_res_list[which(blast_res_list %in% isolate_file)]
  files_to_open <- paste(dir_to_files, files_to_open, sep = "")
  
  #############################################################################
  ## So now we've got the vector of files to open in files_to_open, we'll loop#
  ## through these and get the weight by bitscores, rank of top_pneumo, and ###
  ## bitscore of top pneumo/top hit for each of the files #####################
  #############################################################################
  
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

series_display <- function(list_of_results, prefix, flanks_vector, cluster, region = "total", 
                           reference = TRUE, sd = TRUE){
  
  graph_df <- data.frame(data = matrix(NA, ncol = 9, nrow = (length(list_of_results) * 2)))
  colnames(graph_df) <- c("weighted_bit","ranked_pneumo","bit_pneumo","cluster","flank",
                          "region","sd_1","sd_2", "sd_3")
  
  ## for references flank has to be repped 6 times before input. 
  graph_df$flank <- flanks_vector
  graph_df$cluster <- rep(cluster, nrow(graph_df))
  new_regions <- paste(region, "Control")
  regionz <- rep("ND", length(region)*2)
  for(i in 1:((length(region) * 2))){
    if((i %% 2) == 0)
      regionz[i] <- new_regions[(i/2)]
    else{
      i_val <- as.integer(((i/2) + 0.5))
      regionz[i] <- region[i_val]
    }
    
  }
  
  region <- c(region, new_regions)
  graph_df$region <- rep(regionz, (nrow(graph_df) / length(region)))
  
  k_vec <- seq(1,(length(list_of_results)*2), by = 2)
  
  for(k in 1:length(list_of_results)){
    current_k <- k_vec[k]
    current_res <- list_of_results[[k]]
    if(cluster == "total"){
      data_set <- current_res$total
      ref_ids <- grep("!", data_set$isolate)
      
      
      graph_df$weighted_bit[current_k] <- mean(data_set$weighted_bit[-ref_ids])
      graph_df$sd_1[current_k] <- stats::sd(data_set$weighted_bit[-ref_ids])
      graph_df$ranked_pneumo[current_k] <- mean(data_set$ranked_pneumo[-ref_ids])
      graph_df$sd_2[current_k] <- stats::sd(data_set$ranked_pneumo[-ref_ids])
      graph_df$bit_pneumo[current_k] <- mean(data_set$bit_pneumo[-ref_ids])
      graph_df$sd_3[current_k] <- stats::sd(data_set$bit_pneumo[-ref_ids])
      
      graph_df$weighted_bit[current_k + 1] <- mean(data_set$weighted_bit[ref_ids])
      graph_df$sd_1[current_k + 1] <- stats::sd(data_set$weighted_bit[ref_ids])
      graph_df$ranked_pneumo[current_k+ 1] <- mean(data_set$ranked_pneumo[ref_ids])
      graph_df$sd_2[current_k+ 1] <- stats::sd(data_set$ranked_pneumo[ref_ids])
      graph_df$bit_pneumo[current_k+ 1] <- mean(data_set$bit_pneumo[ref_ids])
      graph_df$sd_3[current_k+ 1] <- stats::sd(data_set$bit_pneumo[ref_ids])
      
      
    }else{
      reference_db <-  current_res$total
      ref_ids <- grep("!", reference_db$isolate)
      reference_db <-  reference_db[ref_ids,]
      
      
      data_set <- current_res$by_insertion[current_res$by_insertion$insertion_loci == cluster,]
      ref_isolates <- paste("ref!", data_set$isolate, sep = "")
      
      reference_isolates <- which(reference_db$isolate %in% ref_isolates)
      reference_db <- reference_db[reference_isolates,]
      
      
      
      graph_df$weighted_bit[current_k] <- mean(data_set$weighted_bit)
      graph_df$sd_1[current_k] <- stats::sd(data_set$weighted_bit)
      graph_df$ranked_pneumo[current_k] <- mean(data_set$ranked_pneumo)
      graph_df$sd_2[current_k] <- stats::sd(data_set$ranked_pneumo)
      graph_df$bit_pneumo[current_k] <- mean(data_set$bit_pneumo)
      graph_df$sd_3[current_k] <- stats::sd(data_set$bit_pneumo)
      
      graph_df$weighted_bit[current_k + 1] <- mean(reference_db$weighted_bit)
      graph_df$sd_1[current_k + 1] <- stats::sd(reference_db$weighted_bit)
      graph_df$ranked_pneumo[current_k+ 1] <- mean(reference_db$ranked_pneumo)
      graph_df$sd_2[current_k+ 1] <- stats::sd(reference_db$ranked_pneumo)
      graph_df$bit_pneumo[current_k+ 1] <- mean(reference_db$bit_pneumo)
      graph_df$sd_3[current_k+ 1] <- stats::sd(reference_db$bit_pneumo)
      
    }
  }
  
  if(sd){
    
    compo_plot_weighted <- ggplot(data = graph_df, aes(x = flank, y = weighted_bit, colour = region)) + geom_point() +
      labs(x = "Flank length", y = "Weighted bit score", title = paste(prefix,  "weighted bit")) + 
      geom_pointrange(aes(ymin = (weighted_bit - sd_1), ymax = (weighted_bit + sd_1)))
    
    compo_plot_ranked <- ggplot(data = graph_df, aes(x = flank, y = ranked_pneumo, colour = region)) + geom_point() +
      labs(x = "Flank length", y = "ranked pneumo", title = paste(prefix,  "ranked pneumo")) + 
      geom_pointrange(aes(ymin = (ranked_pneumo - sd_2), ymax = (ranked_pneumo + sd_2)))
    
    compo_plot_pneumo <- ggplot(data = graph_df, aes(x = flank, y = bit_pneumo, colour = region)) + geom_point() +
      labs(x = "Flank length", y = "ranked pneumo", title = paste(prefix,  "bit pneumo")) + 
      geom_pointrange(aes(ymin = (bit_pneumo - sd_3), ymax = (bit_pneumo + sd_3)))
  }else{
    compo_plot_weighted <- ggplot(data = graph_df, aes(x = flank, y = weighted_bit, colour = region)) + geom_point(size = 2) +
      labs(x = "Flank length", y = "Weighted bit score", title = paste(prefix,  "weighted bit"))
    
    compo_plot_ranked <- ggplot(data = graph_df, aes(x = flank, y = ranked_pneumo, colour = region)) + geom_point(size = 2) +
      labs(x = "Flank length", y = "ranked pneumo", title = paste(prefix,  "ranked pneumo")) 
    
    compo_plot_pneumo <- ggplot(data = graph_df, aes(x = flank, y = bit_pneumo, colour = region)) + geom_point(size = 4) +
      labs(x = "Flank length", y = "Score", title = paste(prefix,  "blast score")) +
      theme_bw()
  }
  
  return(list(graph_data = graph_df, weighted_plot = compo_plot_weighted,
              ranked_plot = compo_plot_ranked, pneumo_plot = compo_plot_pneumo))
  
  
}

folder_to_res <- function(results_folder, flanks_veccy, graph_name){
  ## Function to take in the output from the flanks_only_search and output the over flanks graphs 
  ## Results folder should be the main run folder, flanks veccy the main 
  
  folders <- list.files(results_folder, full.names = TRUE, include.dirs = TRUE)
  whole_blast_res <- folders[grep("_blast_results", folders)]
  before_blast_res <- folders[grep("_before_flank_blast_res", folders)]
  after_blast_res <- folders[grep("_after_flank_blast_res", folders)]
  
  ## whole blast res 
  tic("Running whole flanks sum up")
  pmen3_list <- list()
  
  for(k in whole_blast_res){
    current_graphed_res <- graph_getter(k, graph_name)
    pmen3_list[[length(pmen3_list) + 1]] <- current_graphed_res
    
  }
  
  pmen3_mega_total <- series_display(pmen3_list, graph_name, flanks_veccy, "total")
  toc()
  ## before blast res 
  tic("Before flanks sum up")
  pmen3_list_bef <- list()
  
  for(k in before_blast_res){
    current_graphed_res <- graph_getter(k, graph_name)
    pmen3_list_bef[[length(pmen3_list_bef) + 1]] <- current_graphed_res
    
  }
  
  before_name <- paste(graph_name, "before flanks")
  pmen3_mega_before <- series_display(pmen3_list_bef, before_name, flanks_veccy, "total")
  toc()
  
  ## after blast res 
  tic("After flanks sum up")
  pmen3_list_aft <- list()
  
  for(k in after_blast_res){
    current_graphed_res <- graph_getter(k, graph_name)
    pmen3_list_aft[[length(pmen3_list_aft) + 1]] <- current_graphed_res
    
  }
  
  after_name <- paste(graph_name, "after flanks")
  pmen3_mega_after <- series_display(pmen3_list_aft, after_name, flanks_veccy, "total")
  toc()
  return(list(whole_res = pmen3_mega_total, before_res = pmen3_mega_before, after_res = pmen3_mega_after))
  
}

###############################################################################
## system arguments: 1: End flank length 
##                   2: Directory of blast results 
##                   3: pdf out 
input_args <- commandArgs(trailingOnly = TRUE)

flanks_vec <- seq(500, input_args[1], 500)

pmen_mega_res <- folder_to_res(results_folder = input_args[2],
                               flanks_veccy = flanks_vec, 
                               graph_name = "PMEN MEGA")

pdf(file = input_args[3], paper = "a4r", width = 12, height = 7)

print(pmen_mega_res$whole_res$pneumo_plot)

print(pmen_mega_res$before_res$pneumo_plot)

print(pmen_mega_res$after_res$pneumo_plot)


dev.off()





