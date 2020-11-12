###############################################################################
## Recombination density plots for GPS ########################################
###############################################################################

## Input: - The gubbins results locations
##        - The reccy hits csv
##        - The mge name
##        - Outfile name 

require(ggpubr, quietly = TRUE)
require(stringr, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(tictoc, quietly = TRUE)
require(tools, quietly = TRUE)
require("cowplot", quietly = TRUE)
###############################################################################
## FUNCTIONS ##################################################################
###############################################################################


recombination_gff_cleaner <- function(in_gff_pmen3){

  pmen3_reccy_data <- data.frame(data = (matrix(data = NA, nrow = nrow(in_gff_pmen3),
                                                ncol = 7)))
  colnames(pmen3_reccy_data) <- c("isolate","number_isolates", "snp_count",
                                  "length","density", "start", "end")
  
  for( k in 1:nrow(in_gff_pmen3)){
    current_row <- in_gff_pmen3[k,]
    current_atty <- current_row$attribute
    current_atty_taxa_snp <- str_split_fixed(current_atty,"taxa=",2)[,2]
    taxa <- str_split_fixed(current_atty_taxa_snp,";snp_count=",2)[,1]
    taxa_first <- str_split_fixed(taxa, "_.57195_E01.1",2)[1,1]
    taxa_first <- sub("\\s*","",taxa_first)
    snp_count <- str_split_fixed(current_atty, ";snp_count=",2)[,2]
    snp_count <- as.numeric(sub(";","",snp_count))
    num_taxa <- str_count(taxa, pattern = "_.57195_E01.1")
    length_of_event <- current_row$end - current_row$start + 1
    snp_density <- snp_count / length_of_event
    
    pmen3_reccy_data$isolate[k] <- taxa_first
    pmen3_reccy_data$number_isolates[k] <- num_taxa
    pmen3_reccy_data$snp_count[k] <- snp_count
    pmen3_reccy_data$length[k] <- length_of_event
    pmen3_reccy_data$density[k] <- snp_density
    pmen3_reccy_data$start[k] <- current_row$start
    pmen3_reccy_data$end[k] <- current_row$end 
    
  }
  return(pmen3_reccy_data)
}

reccy_hits_cleaner <- function(reccy_hits, mge_name){

  dims <- nrow(reccy_hits)
  df_mges <- data.frame(matrix(data = NA, nrow = dims, ncol = 6))
  colnames(df_mges) <- c("length", "density","lineage","MGE", "start", "end")
  
  df_mges$lineage <- reccy_hits$cluster_name[1]
  
  
  df_mges$length <- reccy_hits$gub_length
  df_mges$density <- reccy_hits$density
  df_mges$start <- reccy_hits$start_gubb
  df_mges$end <- reccy_hits$end_gubb
  df_mges$MGE <- mge_name
  
  return(df_mges)
  
}

reccy_combiner <- function(whole_csv, mge_csv){

  df_both <- whole_csv[,colnames(whole_csv) %in% c("length", "density", "start", "end")]
  
  df_both$MGE <- rep("No", nrow(df_both))
  
  lengers <- NULL
  counter <- 0
  
  for(k in 1:nrow(df_both)){
    current_row <- df_both[k,]
    subset_mge <- mge_csv[mge_csv$start == current_row$start,]
    if(nrow(subset_mge) > 0){
      counter <- counter + 1
      further_subset <- subset_mge[subset_mge$end == current_row$end,]
      if(nrow(further_subset) >= 1){
        even_further <- further_subset[further_subset$density == current_row$density,]
        if(nrow(even_further) >= 1){
          df_both$MGE[k] <- "Yes"
          lengers <- append(lengers,current_row$bases)
          
        }
        
      }
    }
    
  }
  
  return(df_both)
}

delim_reader <- function(delim_file){
  reccy_csv <- "No"
  try({
    reccy_csv <- read.delim(reccy_gff,header = FALSE, comment.char = "#")
    colnames(reccy_csv) <- c("type","prog","class","start","end","trent","alexander","arnold","attribute")
  }, silent = TRUE)
  return(reccy_csv)
}

input_checker <- function(){
  input_args <- commandArgs(trailingOnly = TRUE)
  
  if(length(input_args) != 4){
    cat("Incorrect number of files, this requires 4, you have:", length(input_args))
    cat("\n")
    cat("Usage: Rscript --vanilla ./recombination_length_plots.R <gubbins_dirs> <reccy_csv> <mge_name> <out_prefix>")
    cat("\n")
    cat("\n")
    stop("Not enough input files")
  }else{
    if(!(file.exists(input_args[2]))){
      cat("The input reccy csv doesn't exist")
      cat("\n")
      cat("Usage: Rscript --vanilla ./recombination_length_plots.R <gubbins_dirs> <reccy_csv> <mge_name> <out_prefix>")
      cat("\n")
      cat("\n")
      stop("Reccy file doesn't exist")
    }
    
    
    return(input_args)
  }
}

###############################################################################
## BEGIN ######################################################################
###############################################################################

input_args <- input_checker()

gubbins_res <- input_args[1]#"~/Dropbox/phd/insertion_site_analysis/data/pmen_run/"
if(tools::file_ext(input_args[2]) == "txt"){
  cat("Using multiple reccy csv files", "\n")
  reccy_files <- readLines(input_args[2])
  reccy_hits <- NULL
  for(file in reccy_files){
    current_hit <- read.csv(file, stringsAsFactors = FALSE)
    reccy_hits <- bind_rows(reccy_hits, current_hit)
  }
}else{
  reccy_hits <- read.csv(input_args[2],
                         stringsAsFactors = FALSE)
}
mge_name <- input_args[3]
out_name <- input_args[4]#"~/Dropbox/phd/insertion_site_analysis/pmen_mega_lib_updated_flanks/pmen_mega_recombination_length.pdf"

out_pdf_file <- paste(out_name, ".pdf", sep = "")
out_png_name <- paste(out_name, ".png", sep = "")
out_csv_name <- paste(out_name, ".csv", sep = "")

## If input args one is just a directory run through the hits with the 
## element in them, if its a file run through all the listed gubbins dirs 
tot_csv <- NULL

if(!(dir.exists(gubbins_res))){
  gubbins_locs <- readLines(gubbins_res)
  cat("This many clusters to work through:", length(gubbins_locs), "\n")
  start_time <- Sys.time()
  
  for(current_dir in gubbins_locs){
    cluster <- basename(current_dir)
    cluster <- sub("_run_data","",cluster)
    tic(paste("curating cluster:", cluster))
    gubb_files <- list.files(current_dir, full.names = TRUE)
    reccy_gff <- gubb_files[grep(".recombination_predictions.gff", gubb_files)]
    narrowed_hits_df <- reccy_hits[reccy_hits$cluster_name == cluster,]
    if(nrow(narrowed_hits_df) > 0){
      narrowed_hits_df$density <- narrowed_hits_df$snp_count / (narrowed_hits_df$gub_length)
      narrowed_hits <- reccy_hits_cleaner(narrowed_hits_df, mge_name)
      
      reccy_csv <- delim_reader(reccy_gff)
      if(class(reccy_csv) == "character"){
        next
      }
      reccy_csv <- recombination_gff_cleaner(reccy_csv)
      
      both_df <- reccy_combiner(reccy_csv, narrowed_hits)
      both_df$cluster <- cluster
      
      tot_csv <- bind_rows(tot_csv, both_df)
    }else{
      
      reccy_csv <- delim_reader(reccy_gff)
      if(class(reccy_csv) == "character"){
        next
      }
      reccy_csv <- recombination_gff_cleaner(reccy_csv)
      df_both <- reccy_csv[,colnames(reccy_csv) %in% c("length", "density", "start", "end")]
      
      df_both$MGE <- rep("No", nrow(df_both))
      df_both$cluster <- cluster
      tot_csv <- bind_rows(tot_csv, df_both)
    }
    toc()
    
  }  
  
}else{

  unique_clusters <- unique(reccy_hits$cluster_name)
  cat("This many clusters to work through:", length(unique_clusters))
  start_time <- Sys.time()
  
  for(cluster in unique_clusters){
    tic(paste("curating cluster:", cluster))
    current_dir <- paste(gubbins_res,cluster, "_run_data", sep = "")
    gubb_files <- list.files(current_dir, full.names = TRUE)
    reccy_gff <- gubb_files[grep(".recombination_predictions.gff", gubb_files)]
    narrowed_hits_df <- reccy_hits[reccy_hits$cluster_name == cluster,]
    narrowed_hits_df$density <- narrowed_hits_df$snp_count / (narrowed_hits_df$gub_length)
    narrowed_hits <- reccy_hits_cleaner(narrowed_hits_df, mge_name)
    
    reccy_csv <- delim_reader(reccy_gff)
    if(class(reccy_csv) == "character"){
      next
    }
    reccy_csv <- recombination_gff_cleaner(reccy_csv)
    
    both_df <- reccy_combiner(reccy_csv, narrowed_hits)
    both_df$cluster <- cluster
    
    tot_csv <- bind_rows(tot_csv, both_df)
    toc()
  }

}

tic("Plotting out results")
histo_mges <- ggplot(data = tot_csv) + geom_histogram(aes(tot_csv$density,
                                                          fill = MGE), alpha  = 1,
                                                      colour = "black") +
  scale_x_log10(limits = c(0.0001, 1)) + scale_y_log10()+
  labs(title = paste(mge_name,"SNP density of event"), y = "Count", x = "Density")

histo_mges_length <- ggplot(data = tot_csv) + geom_histogram(aes(tot_csv$length,
                                                                 fill = MGE), alpha = 1,
                                                             colour = "black") +
  scale_x_log10(limits = c(1, 1e5)) + scale_y_log10() +
  labs(title = paste(mge_name,"Length of event"), y = "Count", x = "Bases")

compo_dot_plot_mges <- ggplot(data = tot_csv) + geom_point(aes(x = length, y = density, colour = MGE, alpha = MGE))

compo_dot_plot_log_mges <- compo_dot_plot_mges + scale_x_log10(limits = c(1,1e5)) + scale_y_log10(limits = c(0.001, 1)) +
  labs(title = paste(mge_name,"SNP density against size of insert"), y = "Density", x = "Bases") + 
  stat_density2d(aes(x = length, y = density))


sum_up_by_MGE <- ggarrange(histo_mges_length, histo_mges, compo_dot_plot_log_mges,
                           ncol = 1, nrow = 3, align = "hv")


pdf(file = out_pdf_file, paper = "a4r", width = 11, height = 7)
print(compo_dot_plot_log_mges)
print(sum_up_by_MGE)

dev.off()

histo_mges <- ggplot(data = tot_csv) + geom_histogram(aes(tot_csv$density,
                                                          fill = MGE), alpha  = 1,
                                                      colour = "black") +
  scale_x_log10(limits = c(0.0001, 1)) + scale_y_log10()+
  labs(y = "Count", x = "Density")

histo_mges_length <- ggplot(data = tot_csv) + geom_histogram(aes(tot_csv$length,
                                                                 fill = MGE), alpha = 1,
                                                             colour = "black") +
  scale_x_log10(limits = c(1, 1e5)) + scale_y_log10() +
  labs( y = "Count", x = "Bases")

compo_dot_plot_mges <- ggplot(data = tot_csv) + geom_point(aes(x = length, y = density, colour = MGE, alpha = MGE))

compo_dot_plot_log_mges <- compo_dot_plot_mges + scale_x_log10(limits = c(1,1e5)) + scale_y_log10(limits = c(0.001, 1)) +
  labs( y = "Density", x = "Bases") + 
  stat_density2d(aes(x = length, y = density))



png(filename = out_png_name, width = 17, height = 15, units = "cm", res = 1000)
ggdraw() +
  draw_plot(histo_mges_length, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(histo_mges, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(compo_dot_plot_log_mges, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))

dev.off()


write.csv(tot_csv, file = out_csv_name, row.names = FALSE, quote = FALSE)

toc()
end_time <- Sys.time()
cat("Finished in", (end_time - start_time))

