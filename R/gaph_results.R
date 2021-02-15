###############################################################################
## Outputting flank scores for different cluster scores #######################
###############################################################################

require(ggplot2, quietly = TRUE)
require(stringr, quietly = TRUE)
require(ggpubr, quietly = TRUE)
require(tictoc, quietly = TRUE)
require(dplyr, quietly = TRUE)

graph_getter <- function(dir_to_files,graph_name, current_flanks, insert_name){
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

  counted_summo <- plyr::count(summo_csv$insertion_point)  
  
  if(insert_name != "total"){
    counted_summo <- counted_summo[counted_summo$x %in% c(insert_name, "reference"),]
    narrowed_refs <- summo_csv[summo_csv$insertion_point == insert_name, 1]
    refs_to_use <- paste(narrowed_refs, "$",sep = "")
    summo_csv <- summo_csv[grepl(paste(refs_to_use, collapse = "|"), summo_csv[,1]),]
  }
  
  
  #############################################################################
  ## Ok, so now we've got the summary_csv, lets go through each of the ########
  ## insertion sites, getting overall pneumo weight by bitscore, rank of top ##
  ## pneumo hit and then top hit bit / bit of top pneumo hit ##################
  #############################################################################
  
  for(k in 1:nrow(counted_summo)){
    #browser()
    tic(paste(as.character(k),"row done"))
    current_insertion <- as.character(counted_summo[k, 1])
    subsetted_data <- summo_csv[summo_csv$insertion_point == current_insertion,]
    
    if(k == 1){
      results_df <- bitscore_res(subsetted_data, blast_results,
                                 current_insertion, dir_to_files, current_flanks)
      toc()
    }else{
      new_df <- bitscore_res(subsetted_data, blast_results,
                             current_insertion, dir_to_files, current_flanks)
      results_df <- rbind.data.frame(results_df, new_df)
      toc()
    }
    
  }
  
  total_run_through <- bitscore_res(summo_csv, blast_results,
                                    "total", dir_to_files, current_flanks)
  
  #browser()
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

bitscore_res <- function(data_subset, blast_res_list, insertion_loci, dir_to_files, current_flanks){
  
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
                           reference = TRUE, sd = FALSE){
  
  graph_df <- data.frame(data = matrix(NA, ncol = 9, nrow = (length(list_of_results) * 2)))
  colnames(graph_df) <- c("weighted_bit","ranked_pneumo","bit_pneumo","cluster","flank",
                          "region","lqr","uqr", "sd_3")
  
  graph_df2 <- NULL
  
  ## for references flank has to be repped 6 times before input. 
  graph_df$flank <- rep(flanks_vector, each = 2)
  graph_df$cluster <- rep(cluster, nrow(graph_df))
  new_regions <- paste(region, "Control")
  region <- paste(cluster, region)
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
      data_set2 <- current_res$by_insertion
      ref_ids <- grep("!", data_set$isolate)
      ref_ids2 <- grep("!", data_set2$isolate)
      
      
      graph_df$weighted_bit[current_k] <- mean(data_set$weighted_bit[-ref_ids])
      graph_df$lqr[current_k] <- quantile(data_set$bit_pneumo[-ref_ids], probs = 0.05)
      graph_df$ranked_pneumo[current_k] <- mean(data_set$ranked_pneumo[-ref_ids])
      graph_df$uqr[current_k] <- quantile(data_set$bit_pneumo[-ref_ids], probs = 0.95)
      graph_df$bit_pneumo[current_k] <- median(data_set$bit_pneumo[-ref_ids])
      graph_df$sd_3[current_k] <- stats::sd(data_set$bit_pneumo[-ref_ids])
      
      graph_df$weighted_bit[current_k + 1] <- mean(data_set$weighted_bit[ref_ids])
      graph_df$lqr[current_k + 1] <- quantile(data_set$bit_pneumo[ref_ids], probs = 0.05)
      graph_df$ranked_pneumo[current_k+ 1] <- mean(data_set$ranked_pneumo[ref_ids])
      graph_df$uqr[current_k+ 1] <- quantile(data_set$bit_pneumo[ref_ids], probs = 0.95)
      graph_df$bit_pneumo[current_k+ 1] <- median(data_set$bit_pneumo[ref_ids])
      graph_df$sd_3[current_k+ 1] <- stats::sd(data_set$bit_pneumo[ref_ids])
      
      graphing_rows <- data_set2[,c("isolate","bit_pneumo", "insertion_loci")]
      graphing_rows$control <- "Actual"
      graphing_rows$control[ref_ids2] <- "Control"
      graph_df2 <- bind_rows(graph_df2, graphing_rows)
      
      
    }else{

      reference_db <-  current_res$total
      ref_ids <- grep("!", reference_db$isolate)
      reference_db <-  reference_db[ref_ids,]
      
      
      data_set <- current_res$by_insertion[current_res$by_insertion$insertion_loci == cluster,]
      ref_isolates <- paste("!", data_set$isolate, sep = "")
      ref_matches <- unique(grep(paste(ref_isolates, collapse = "|"), reference_db$isolate))
      
      reference_db <- reference_db[ref_matches,]
      reference_db$insertion_loci <- cluster
      data_set2 <- bind_rows(data_set, reference_db)
      ref_ids2 <- grep("!",data_set2$isolate)
#      browser()
      graph_df$weighted_bit[current_k] <- mean(data_set$weighted_bit)
      graph_df$lqr[current_k] <- quantile(data_set$bit_pneumo, probs = 0.25)[1]
      graph_df$ranked_pneumo[current_k] <- mean(data_set$ranked_pneumo)
      graph_df$uqr[current_k] <- quantile(data_set$bit_pneumo, probs = 0.75)[1]
      graph_df$bit_pneumo[current_k] <- median(data_set$bit_pneumo)
      graph_df$sd_3[current_k] <- stats::sd(data_set$bit_pneumo)
      
      graph_df$weighted_bit[current_k + 1] <- mean(reference_db$weighted_bit)
      graph_df$lqr[current_k + 1] <- quantile(reference_db$bit_pneumo, probs = 0.25)
      graph_df$ranked_pneumo[current_k+ 1] <- mean(reference_db$ranked_pneumo)
      graph_df$uqr[current_k+ 1] <- quantile(reference_db$bit_pneumo, probs = 0.75)
      graph_df$bit_pneumo[current_k+ 1] <- median(reference_db$bit_pneumo)
      graph_df$sd_3[current_k+ 1] <- stats::sd(reference_db$bit_pneumo)
      
      
      graphing_rows <- data_set2[,c("isolate","bit_pneumo", "insertion_loci")]
      graphing_rows$control <- "Actual"
      graphing_rows$control[ref_ids2] <- "Control"
      graph_df2 <- bind_rows(graph_df2, graphing_rows) 
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
      labs(x = "Flank length", y = "Weighted bit score") + theme(legend.position = "none")#, title = paste(prefix,  "weighted bit"))
    
    compo_plot_ranked <- ggplot(data = graph_df, aes(x = flank, y = ranked_pneumo, colour = region)) + geom_point(size = 2) +
      labs(x = "Flank length", y = "ranked pneumo") + theme(legend.position = "none")#, title = paste(prefix,  "ranked pneumo")) 
    
    compo_plot_pneumo <- ggplot(data = graph_df, aes(x = flank, y = bit_pneumo, colour = region)) + geom_point(size = 2) +
      geom_line() + geom_ribbon(aes(ymin = lqr, ymax = uqr, fill = region), colour = NA, alpha = 0.25) +
      labs(x = "Flank length", y = "Score") + theme(legend.position = "none") 
      
  }
  
  return(list(graph_data = graph_df, weighted_plot = compo_plot_weighted,
              ranked_plot = compo_plot_ranked, pneumo_plot = compo_plot_pneumo, total_df = graph_df2))
  
  
}

folder_to_res <- function(results_folder, graph_name, insert_name){
  ## Function to take in the output from the flanks_only_search and output the over flanks graphs 
  ## Results folder should be the main run folder, flanks veccy the main 
  
  folders <- list.files(results_folder, full.names = TRUE, include.dirs = TRUE)
  folders <- sub("//","/",folders)
  whole_blast_res <- folders[grep("_blast_results", folders)]
  before_blast_res <- folders[grep("_before_flank_blast_res", folders)]
  after_blast_res <- folders[grep("_after_flank_blast_res", folders)]
  
  ## whole blast res 
  tic("Running whole flanks sum up")
  pmen3_list <- list()
  flanks_veccy <- as.integer(sub("_blast_results","",basename(whole_blast_res)))
  current_iter <- 1
  for(k in whole_blast_res){
    #next
    current_graphed_res <- graph_getter(k, graph_name, flanks_veccy[current_iter], insert_name)
    pmen3_list[[length(pmen3_list) + 1]] <- current_graphed_res
    current_iter <- current_iter + 1
  }

  pmen3_mega_total <- series_display(pmen3_list, graph_name, flanks_veccy, insert_name, region = "whole")
  toc()
  ## before blast res 
  tic("Before flanks sum up")
  pmen3_list_bef <- list()
  flanks_veccy <- as.integer(sub("_before_flank_blast_res","",basename(before_blast_res)))
  current_iter <- 1
  for(k in before_blast_res){
    #next
    
    current_graphed_res <- graph_getter(k, graph_name, flanks_veccy[current_iter], insert_name)
    pmen3_list_bef[[length(pmen3_list_bef) + 1]] <- current_graphed_res
    current_iter <- current_iter + 1
  }
  
  before_name <- paste(graph_name, "before flanks")
  pmen3_mega_before <- series_display(pmen3_list_bef, before_name, flanks_veccy, insert_name, region = "before")
  toc()
  
  ## after blast res 
  tic("After flanks sum up")
  pmen3_list_aft <- list()
  flanks_veccy <- as.integer(sub("_after_flank_blast_res","",basename(after_blast_res)))
  current_iter <- 1
  for(k in after_blast_res){
    
    current_graphed_res <- graph_getter(k, graph_name, flanks_veccy[current_iter], insert_name)
    pmen3_list_aft[[length(pmen3_list_aft) + 1]] <- current_graphed_res
    current_iter <- current_iter + 1
  }
  
  after_name <- paste(graph_name, "after flanks")
  pmen3_mega_after <- series_display(pmen3_list_aft, after_name, flanks_veccy, insert_name, region = "after")
  toc()
  return(list(whole_res = pmen3_mega_total, before_res = pmen3_mega_before, after_res = pmen3_mega_after))
  
}

parse_args <- function(){
  input_args <- commandArgs(trailingOnly = TRUE)
  
  if(length(input_args) != 3){
    cat("Not enough arguments, need 3, you have:", length(input_args), "\n")
    cat("\n")
    cat("Usage: Rscript --vanilla ./gaph_results.R <Graph name> <folder_of_flank_res> <pdf_out_file>")
    cat("\n")
    cat("\n")
    stop("Not enough input files")
  }
  
  return(input_args)
}

###############################################################################
## system arguments: 1: Graph name 
##                   2: Directory of blast results 
##                   3: pdf out 


input_args <- parse_args()

graph_name <- input_args[1]
graph_name <- "GPS MEGA"

input_args[2] <- "~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/mega_res/gps_mega_run4/"

pmen_mega_res_total <- folder_to_res(results_folder = input_args[2],
                               graph_name = graph_name,
                               insert_name = "total")

## Get some info out from the gamma scores 
mega_whole_res <- pmen_mega_res_total$whole_res$total_df

## MGE insert sites 
inserts_whole <- mega_whole_res[mega_whole_res$control == "Actual",]
mean(inserts_whole$bit_pneumo)
reference_whole <- mega_whole_res[mega_whole_res$control == "Control",]
mean(reference_whole$bit_pneumo)

wilcox.test(inserts_whole$bit_pneumo, reference_whole$bit_pneumo)
t.test(inserts_whole$bit_pneumo, reference_whole$bit_pneumo)

ggplot(data = mega_whole_res, aes(x = control, y = bit_pneumo, fill = control, colour = control)) +
  geom_violin()

whole_df <- bind_rows(pmen_mega_res_total$whole_res$total_df , bind_rows(pmen_mega_res_total$before_res$total_df, pmen_mega_res_total$after_res$total_df))

whole_df_insert <- whole_df[whole_df$control == "Actual",]
whole_df_ref <- whole_df[whole_df$control == "Control",]
wilcox.test(whole_df_insert$bit_pneumo, whole_df_ref$bit_pneumo)
median(whole_df_insert$bit_pneumo)
median(whole_df_ref$bit_pneumo)

boxplot_whole <- ggplot(data = whole_df,aes(x = control, y = bit_pneumo)) + geom_boxplot() 
boxplot_whole

wilcox.test()

whole_df <- whole_df %>% mutate(control = ifelse(control == "Actual", "Tag","Control"))
whole_df$control <- factor(whole_df$control, levels = c("Tag","Control"))

violin_plot_whole <- ggplot(data = whole_df, aes(x = control, y = bit_pneumo)) + geom_violin(aes(fill = control)) +
  labs(x = "Isolate", y = "Score") + scale_fill_discrete(breaks = c("Actual","Control")) + 
  guides(fill = guide_legend(title = "Isolate"))

violin_plot_whole

whole_plot <- pmen_mega_res$whole_res$pneumo_plot

###############################################################################
## Run it through for Tn916 ###################################################
###############################################################################

pmen_tn916_res_total <- folder_to_res(results_folder = "~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/tn916_res/gps_tn916_run_free/",
                                     graph_name = graph_name,
                                     insert_name = "total")

gps_tn916_total <- pmen_tn916_res_total$whole_res$total_df
inserts_whole <- gps_tn916_total[gps_tn916_total$control == "Actual",]
mean(inserts_whole$bit_pneumo)
reference_whole <- gps_tn916_total[gps_tn916_total$control == "Control",]
mean(reference_whole$bit_pneumo)

wilcox.test(inserts_whole$bit_pneumo, reference_whole$bit_pneumo)

t.test(inserts_whole$bit_pneumo, reference_whole$bit_pneumo)


whole_df <- bind_rows(pmen_tn916_res_total$whole_res$total_df , bind_rows(pmen_tn916_res_total$before_res$total_df, pmen_tn916_res_total$after_res$total_df))

whole_df_insert <- whole_df[whole_df$control == "Actual",]
whole_df_ref <- whole_df[whole_df$control == "Control",]

## Lets pair it up 
pairing_df <- whole_df %>% mutate(iso = ifelse(grepl("!",isolate),sub(".*!","",isolate), isolate))

wilcox.test(whole_df_insert$bit_pneumo, whole_df_ref$bit_pneumo)
median(whole_df_insert$bit_pneumo)
median(whole_df_ref$bit_pneumo)

boxplot_whole <- ggplot(data = whole_df,aes(x = control, y = bit_pneumo)) + geom_boxplot() 
boxplot_whole

wilcox.test()

whole_df <- whole_df %>% mutate(control = ifelse(control == "Actual", "Tag","Control"))
whole_df$control <- factor(whole_df$control, levels = c("Tag","Control"))

violin_plot_whole <- ggplot(data = whole_df, aes(x = control, y = bit_pneumo)) + geom_violin(aes(fill = control)) +
  labs(x = "Isolate", y = "Score") + scale_fill_discrete(breaks = c("Actual","Control")) + 
  guides(fill = guide_legend(title = "Isolate"))

violin_plot_whole

whole_plot <- pmen_mega_res$whole_res$pneumo_plot

###############################################################################

pmen_mega_res_tag <- folder_to_res(results_folder = input_args[2],
                                     graph_name = "GPS Tn1207.1",
                                     insert_name = "28")
## Get some info on the insertions 
## Look into the top species hits for tag at 500bp upstream. 


pmen_mega_res_tag$whole_res$graph_data
whole_df <- bind_rows(pmen_mega_res_tag$before_res$total_df, pmen_mega_res_tag$after_res$total_df)
whole_df <- whole_df %>% mutate(control = ifelse(control == "Actual", "Tag","Control"))
whole_df$control <- factor(whole_df$control, levels = c("Tag","Control"))

violin_plot_whole <- ggplot(data = whole_df, aes(x = control, y = bit_pneumo)) + geom_violin(aes(fill = control)) +
  labs(x = "Isolate", y = "Score") + scale_fill_discrete(labels = c(expression(paste(italic("tag"), " insertion",sep = "")),"Unmodified Locus")) + 
  guides(fill = guide_legend(title = "Isolate")) + geom_point(aes(x = control, y = bit_pneumo, color = control),
                                                              position = position_jitter(width = 0.45, height = 0),
                                                              size = 0.25) +
  scale_color_discrete(labels = c(expression(paste(italic("tag"), " insertion",sep = "")),"Unmodified Locus"), name = "Isolate") + 
  theme(legend.text.align = 0) +
  scale_x_discrete(labels = c(expression(paste(italic("tag"), " insertion",sep = "")),"Unmodified Locus")) +
  theme_bw()

violin_plot_whole

#whole_plot <- pmen_mega_res$whole_res$pneumo_plot


mean_40_60 <- ggdraw()  +
  draw_plot(pmen_mega_res_tag$before_res$pneumo_plot + labs(x = NULL) + theme_bw() + theme(legend.position = "none"), x = 0, y = 0.5, width = .4, height = .5) +
  draw_plot(pmen_mega_res_tag$after_res$pneumo_plot  + labs(x = "Flank Length (bp)")  + theme_bw() + theme(legend.position = "none"), x = 0, y = 0, width = 0.4, height = 0.5) +
  draw_plot(violin_plot_whole + theme(legend.text.align = 0), x = 0.4, y = 0, width = 0.6, height = 1) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0, 0.4), y = c(0.95, 0.45, 0.95)) 

mean_40_60
  mean_95 <- ggdraw()  +
  draw_plot(pmen_mega_res$whole_res$pneumo_plot, x = 0, y = 0.6, width = .5, height = .3) +
  draw_plot(pmen_mega_res$before_res$pneumo_plot, x = 0, y = 0.3, width = .5, height = .3) +
  draw_plot(pmen_mega_res$after_res$pneumo_plot, x = 0, y = 0, width = 0.5, height = 0.3) +
  draw_plot(violin_plot_whole, x = 0.5, y = 0, width = 0.5, height = 0.9) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0, 0, 0.5), y = c(0.9, 0.6, 0.3, 0.9)) + draw_plot_label("Mean score, 95% range", size = 15,
                                                                                    x = 0.2, y = 1)

mean_95

median_iqr <- ggdraw()  +
  draw_plot(pmen_mega_res$whole_res$pneumo_plot, x = 0, y = 0.6, width = .5, height = .3) +
  draw_plot(pmen_mega_res$before_res$pneumo_plot, x = 0, y = 0.3, width = .5, height = .3) +
  draw_plot(pmen_mega_res$after_res$pneumo_plot, x = 0, y = 0, width = 0.5, height = 0.3) +
  draw_plot(violin_plot_whole, x = 0.5, y = 0, width = 0.5, height = 0.9) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0, 0, 0.5), y = c(0.9, 0.6, 0.3, 0.9)) + draw_plot_label("Median score, IQR range", size = 15,
                                                                                    x = 0.2, y = 1)

median_iqr

median_95 <- ggdraw()  +
  draw_plot(pmen_mega_res$whole_res$pneumo_plot, x = 0, y = 0.6, width = .5, height = .3) +
  draw_plot(pmen_mega_res$before_res$pneumo_plot, x = 0, y = 0.3, width = .5, height = .3) +
  draw_plot(pmen_mega_res$after_res$pneumo_plot, x = 0, y = 0, width = 0.5, height = 0.3) +
  draw_plot(violin_plot_whole, x = 0.5, y = 0, width = 0.5, height = 0.9) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0, 0, 0.5), y = c(0.9, 0.6, 0.3, 0.9)) + draw_plot_label("Median score, 95% range", size = 15,
                                                                                    x = 0.2, y = 1)

median_95

pdf(file = "~/Dropbox/phd/LSA/figure_4_options.pdf", paper = "a4r", width = 12, height = 7)

print(mean_iqr)

print(mean_95)

print(median_iqr)

print(median_95)

dev.off()




mean_95 <- plot_grid(draw_label("Mean score, 95% range") , mean_95, ncol = 1, rel_heights = c(0.1,1))


median_95 <- ggdraw() +
   draw_plot(pmen_mega_res$whole_res$pneumo_plot, x = 0, y = 0.67, width = .5, height = .33) +
   draw_plot(pmen_mega_res$before_res$pneumo_plot, x = 0, y = 0.34, width = .5, height = .33) +
   draw_plot(pmen_mega_res$after_res$pneumo_plot, x = 0, y = 0, width = 0.5, height = 0.33) +
   draw_plot(violin_plot_whole, x = 0.5, y = 0, width = 0.5, height = 1) +
   draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                   x = c(0, 0, 0, 0.5), y = c(1, 0.67, 0.33, 1))

pdf(file = input_args[3], paper = "a4r", width = 12, height = 7)

print(pmen_mega_res$whole_res$pneumo_plot)

print(pmen_mega_res$before_res$pneumo_plot)

print(pmen_mega_res$after_res$pneumo_plot)


dev.off()



pmen_mega_res_tag <- folder_to_res(results_folder = input_args[2],
                                   graph_name = "GPS Tn1207.1",
                                   insert_name = "28")


tag_res <- graph_getter(dir_to_files = "~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/mega_res/gps_mega_run4/500_before_flank_blast_res/",
                        graph_name = "before_tag", insert_name = "28",current_flanks = "500")


###############################################################################
## Check for Tn916 ############################################################
###############################################################################






















