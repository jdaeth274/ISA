###############################################################################
## This is our script to output a merged blast csv from an initial blast csv ##
###############################################################################

require(stringr, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(gtools, quietly = TRUE, warn.conflicts = FALSE)

###############################################################################
## First lets read in the csv from the blast search ###########################
###############################################################################

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


in_files <- commandArgs(trailingOnly = TRUE)

if(length(in_files) != 3){
  print("Not correct number of inputs, need 2: <blast_res> <contig_dir> <out_csv_name>")
  quit(save = "no", status = 1, runLast = FALSE)
}

in_csv_file <- commandArgs(trailingOnly = TRUE)[1]
#in_csv_file <- "~/Documents/phd_downloads/gps_downloads/cluster_100_list_gffs/tmp_dna_dir/tmp_blast_csv.csv"
contig_dir <- commandArgs(trailingOnly = TRUE)[2]
#contig_dir <- "~/Documents/phd_downloads/gps_downloads/cluster_100_list_gffs/contig_bounds/"



if(substrRight(contig_dir, 1) != "/")
  contig_dir <- paste(contig_dir, "/", sep = "")


out_csv_file_name <- commandArgs(trailingOnly = TRUE)[3]
in_csv <- read.csv(in_csv_file, stringsAsFactors = FALSE, header = FALSE)
colnames(in_csv) <- c("Reference", "subject",	"pident",	"align",
                      "mismatch",	"gapopen",	"qstart",	"qend",
                      "sstart",	"ssend",	"eval",	"bitscore")

###############################################################################
## Lets now input the function for mergin our in_csv ##########################
###############################################################################

contig_checker <- function(isolate_blast_rows, contig_table){
  isolate_blast_rows$contig <- rep(0, nrow(isolate_blast_rows))
  for(k in 1:nrow(isolate_blast_rows)){
    subject_locs <- as.numeric(as.vector(isolate_blast_rows[k, c(6,7)]))
    subject_locs <- subject_locs[order(subject_locs)]
    current_contig <- which((contig_table[,1] - 15) <= subject_locs[1] & (contig_table[,2] + 15) >= subject_locs[2])
    isolate_blast_rows$contig[k] <- current_contig 
  }
  return(isolate_blast_rows)
}


duplicated_hit_checker <- function(merged_blast_rows){
  dup_row <- NULL
  for(k in 1:nrow(merged_blast_rows)){
    
    current_locs <- as.vector(as.integer(merged_blast_rows[k, c(4,5)]))
    ident_test <- which((merged_blast_rows$qstart == current_locs[1]) & (merged_blast_rows$qend == current_locs[2]))
    
    if(length(ident_test) > 1){
      if(any(sapply(dup_row, `%in%`, ident_test)) == FALSE)   
        dup_row[[length(dup_row) + 1]] <- ident_test
    }
    
    
  }
  
  if(length(dup_row) > 0){
    losers <- NULL
    for(k in 1:length(dup_row)){
      current_similar <- dup_row[[k]]
      largest_contig_hit <- which.min(merged_blast_rows$sstart[current_similar])
      current_lose <- current_similar[-largest_contig_hit]
      losers <- append(losers, current_lose)
    }
    
    merged_blast_rows <- merged_blast_rows[-losers,]
  }
  
  return(merged_blast_rows)
  
}

contig_tries <- function(accession, contig_dir, try_type){
  try_type <- 1
  contig_tab <- NULL
  print("About to try try_type 1")
  try(contig_tab <- contig_file_read(accession = accession, contig_dir = contig_dir,
                                     try_type = 1), silent = TRUE)
  if(length(contig_tab) == 0){
    print("About to try try_type 2")
    try(contig_tab <- contig_file_read(accession = accession, contig_dir = contig_dir,
                                       try_type = 2), silent = TRUE)
    if(length(contig_tab) == 0){
      stop("Can't recognise contig file name format")
    }
  }
  
  return(try_type)
  
}

contig_file_read <- function(accession, contig_dir, try_type){
  
  if(try_type == 1){
  
  contig_accession <- sub(pattern = "#",replacement = "_",
                          accession)
  contig_file <- paste(contig_dir, contig_accession, "#contig_bounds.csv",sep = "")
  contig_tab <- read.csv(contig_file)
  }else if(try_type == 2){
    print("In try type 2")
    print(accession)
    contig_accession <- sub(pattern = "\\.contigs_velvet.fa",replacement = "",
                            accession)
    contig_accession <- sub(pattern = "#",replacement = "_",
                            contig_accession)
    print(contig_accession)
    contig_file <- paste(contig_dir, contig_accession, "#contig_bounds.csv",sep = "")
    contig_tab <- read.csv(contig_file)
  }
  
  return(contig_tab)
}



delete_dups_blast <- function(hit_csv, contig_dir){

  subject_ids <- hit_csv$subject
  # subject_ids <- str_split_fixed(subject_ids, ".c",2)[, 1]
  # subject_ids <- str_split_fixed(subject_ids, "\\.",2)[, 2]
  # subject_ids <- str_split_fixed(subject_ids, "\\.",2)[,1]
  # 
  # subject_ids <- sub("\\./","",subject_ids)
  
  subject_ids <- str_split_fixed(subject_ids,"\\.", 3)[,1]
  print(head(subject_ids))
  
  for( i in 1:length(subject_ids)){
    current_id <- subject_ids[i]
    if((str_count(current_id, "_")) == 2){
      second_undo <- str_locate_all(current_id, "_")[[1]][2,1]
      substr(current_id, second_undo, second_undo) <- "#"
      subject_ids[i] <- current_id
    }
  }
  print(head(subject_ids))
  subject_count <- dplyr::count(as.data.frame(subject_ids), subject_ids)
  max_accession_frequency <- max(dplyr::count(as.data.frame(subject_ids), subject_ids)[,2])
  
  dups_df <- hit_csv[,c(2:4,7:10,12)]
  dups_df$subject <- subject_ids
  dups_df$file <- hit_csv$subject
  ###############################################################################
  ## Now we remove any nested sequences from the csv ############################
  ###############################################################################
  
  loss_row <- NULL
  loss_row_many <- NULL
  del_row <- NULL
  for(i in 1:nrow(subject_count)){
    duplicate <- dups_df[dups_df$subject==subject_count[i,1],] 
    number_hits <- subject_count[i, 2]
    
    if(number_hits > 1){
      perms_mat <- gtools::permutations(n = number_hits,
                                        r = 2,
                                        v = c(1:number_hits))
      
      for(ii in 1:nrow(perms_mat)){
        if((duplicate[perms_mat[ii,1],4]>=duplicate[perms_mat[ii,2],4])
           &(duplicate[perms_mat[ii,1],5]<=duplicate[perms_mat[ii,2],5])&
           (duplicate[perms_mat[ii,1],3]!=duplicate[perms_mat[ii,2],3])){
          
          del_row <- duplicate[perms_mat[ii,1],]
          loss_row_many <- rbind(loss_row_many,del_row)
        }
        
      }
      loss_row <- unique(loss_row_many)
    }
    
    #print(i/nrow(subject_count)*100)
  }
  
  ###############################################################################
  ## Now we will get rid of these nested values from the csv ####################
  ###############################################################################
  if(length(loss_row) > 0){
    row_nums <- as.numeric(rownames(loss_row))
    data_no_nest <- dups_df[-c(row_nums),]
    nest_accession <- dplyr::count(data_no_nest, subject)
  }else{
    data_no_nest <- dups_df
    nest_accession <- dplyr::count(data_no_nest, subject)
  }
  
  ###############################################################################
  ## Now lets make a note of the exact match seqeunces that are dotted around ###
  ## the genome #################################################################
  ###############################################################################
  ident_row=NULL
  ident_row_many=NULL
  double_hits_df <- NULL
  for (i in 1:nrow(nest_accession)){
    
    
    duplicate=data_no_nest[data_no_nest$subject==nest_accession[i,1],]
    hit_number <- nest_accession[i, 2]
    if (nest_accession[i,2]>1) {                                             ###This Isolates duplicate data sets
      #creating the combination matrix
      perms_mat <- gtools::permutations(n = hit_number,
                                        r = 2,
                                        v = c(1:hit_number))
      
      #creating the data frame to contain all the row to lose
      ident_row=NULL
      ident_row_many <- NULL
      #creating the fopr loop to compare all the different combinations against each other
      
      for (j in 1:nrow(perms_mat)){
        if((duplicate[perms_mat[j,1],4]==duplicate[perms_mat[j,2],4])
           &(duplicate[perms_mat[j,1],5]==duplicate[perms_mat[j,2],5]))
        {ident_row=duplicate[perms_mat[j,1],]
        ident_row_many=rbind(ident_row_many,ident_row)
        }
      }
      
      ident_row=unique(ident_row_many)
      if(length(ident_row) > 0){
        counters <- dplyr::count(ident_row,align)
        if(nrow(counters) > 1){
          repeat_seq <- c(1:counters[1,2])
          segment_num <- nrow(counters)
          nrow_seq <- c(0, (1:(nrow(counters)-1))*length(repeat_seq))
          new_df <- ident_row[order(ident_row$qstart),]
          a <- sort(rep(repeat_seq,segment_num)) + rep(nrow_seq,length(repeat_seq))
          new_df <- new_df[a, ]
          new_df$clust_num <- sort(rep(repeat_seq,segment_num))
        }else{
          ident_row$clust_num <- c(1:counters[1,2])
          new_df <- ident_row
        }
        double_hits_df <- rbind.data.frame(double_hits_df,
                                           new_df)
        
      }
    }
    #print((i/nrow(nest_accession))*100,digits=3)
  }
  
  ###############################################################################
  ## Lets join up the different dfs #############################################
  ###############################################################################
  data_no_nest$clust_num <- rep("ND",nrow(data_no_nest))
  repped_ids <- as.character(dplyr::count(double_hits_df,subject)[,1])
  if(length(repped_ids > 0)){
    new_double_df <- NULL
    for(k in 1:length(repped_ids)){
      no_nest_data <- data_no_nest[data_no_nest$subject == repped_ids[k], ]
      double_df_data <- double_hits_df[double_hits_df$subject == repped_ids[k], ]
      
      if(nrow(no_nest_data) != nrow(double_df_data)){
        if(no_nest_data$sstart[1] > no_nest_data$ssend[1])
          no_nest_data <- no_nest_data[order(no_nest_data$sstart, decreasing = TRUE),]
        else{
          no_nest_data <- no_nest_data[order(no_nest_data$sstart),]
        }
        no_nest_data$clust_num[1] <- 1
        for(g in 2:nrow(no_nest_data)){
          if(no_nest_data[g,4] > no_nest_data[g-1,5])
            no_nest_data$clust_num[g] <- no_nest_data$clust_num[g-1]
          else
            no_nest_data$clust_num[g] <- as.numeric(no_nest_data$clust_num[g-1]) + 1
        }
        double_df_data <- no_nest_data
      }
      
      new_double_df <- rbind.data.frame(new_double_df, double_df_data)
    }
    
    data_no_nest[data_no_nest$subject %in% repped_ids, ] <- new_double_df
    print("There are some isolates with more than one hit in them")
  }
  data_orientation <- transform(data_no_nest,
                                orientation= ifelse(sstart<ssend,
                                                    "forward", "reverse"))
  
  ###############################################################################
  ## Lets do some merging! ######################################################
  ###############################################################################
  new_row <- NULL
  data_merged <- NULL
  
  for (i in 1:nrow(nest_accession)){
     # if( i == 168){
     #   browser()
     # }
     # 
    if(i == 1){
      try_type <- contig_tries(accession = nest_accession[i, 1],
                               contig_dir = contig_dir,
                               try_type = 1)
    }
    
    contig_tab <- contig_file_read(accession = nest_accession[i, 1],
                                   contig_dir = contig_dir,
                                   try_type = try_type)
    
    clust_rows <- NULL
    new_row <- NULL
          
    
    
    if(nest_accession[i,2]==1){
      ## If there's only one hit in the genome we can just use this row as the 
      ## entry for our db 
      new_row <- data_orientation[data_orientation$subject==nest_accession[i,1],]
      cols_to_lose <- c("clust_num", "contig")
      new_row <- new_row[, -which(colnames(new_row) %in% cols_to_lose)]
      colnames(new_row)[9] <- "file_loc"
    }
    
    if(nest_accession[i,2]>1){
      data_set <- data_orientation[data_orientation$subject==nest_accession[i,1],]
      
      data_set <- contig_checker(data_set, contig_tab)
      
      if(data_set$clust_num[1] != "ND"){
        ## If there are multiple repeats of the same fragments throughout the genome
        ## we have to join those clusters together and then put each as a new entry
        ## into our db 
        clust_rows <- NULL
        new_row <- NULL
        
        contig_nums <- dplyr::count(data_set, contig)
        for(k in 1:nrow(contig_nums)){
          current_contig <- contig_nums[k, 1]
          data_set_contig <- data_set[data_set$contig == current_contig, ]
          
          if(length(unique(data_set_contig$clust_num)) > 1){
            
            data_set_contig <- duplicated_hit_checker(data_set_contig)
            
          }
          
          subject <- paste(data_set_contig[1,1], "_", as.character(k), sep = "")                                             ###outcome from a column
          pident <- (sum(data_set_contig$pident*data_set_contig$align))/(sum(data_set_contig$align))
          align <- sum(data_set_contig$align)
          if(data_set_contig$orientation[1] == "forward"){
                      qstart <- min(data_set_contig$qstart)
                      qend <- max(data_set_contig$qend)
                      sstart <- min(data_set_contig$sstart)
                      ssend <- max(data_set_contig$ssend)
                    }else{
                      qstart <- min(data_set_contig$qstart)
                      qend <- max(data_set_contig$qend)
                      sstart <- max(data_set_contig$sstart)
                      ssend <- min(data_set_contig$ssend)
                    }
          bitscore <- sum(data_set_contig$bitscore)
          #ifelse("y" %in% data_set$mefA == TRUE, mefA<-"y", mefA<- "n")
          orientation <- data_set_contig$orientation[1]
          file_loc <- data_set_contig$file[1]
          clust_rows <- data.frame(subject, pident, align, qstart, qend, sstart,
                                             ssend, bitscore, file_loc, orientation)
                    
          new_row <- rbind.data.frame(new_row, clust_rows)
    
        }
        
        
        
        new_row <- duplicated_hit_checker(new_row)
        
        # 
        # 
        # for(j in 1:max(data_set$clust_num)){
        #   narrow_data <- data_set[data_set$clust_num == j,]
        #   subject <- paste(narrow_data[1,1],"_",j,sep = "")
        #   pident <- (sum(narrow_data$pident*narrow_data$align))/(sum(narrow_data$align))
        #   align <- sum(narrow_data$align)
        #   if(narrow_data$orientation[1] == "forward"){
        #     qstart <- min(narrow_data$qstart)
        #     qend <- max(narrow_data$qend)
        #     sstart <- min(narrow_data$sstart)
        #     ssend <- max(narrow_data$ssend)
        #   }else{
        #     qstart <-min(narrow_data$qstart)
        #     qend <-max(narrow_data$qend)
        #     sstart <- max(narrow_data$sstart)
        #     ssend <- min(narrow_data$ssend)
        #   }
        #   bitscore <- sum(narrow_data$bitscore)
        #   #ifelse("y" %in% data_set$mefA == TRUE, mefA<-"y", mefA<- "n")
        #   orientation <- narrow_data$orientation[1]
        #   file_loc <- narrow_data$file[1]
        #   clust_rows <- data.frame(subject, pident, align, qstart, qend, sstart,
        #                            ssend, bitscore, file_loc, orientation)
        #   if(abs(clust_rows$sstart - clust_rows$ssend) >= 40000){
        #     ## Sometimes though we can join up hits that are too far apart to 
        #     ## realistically represent a single entry into the genome. Now 
        #     ## we have to look at the individual hits and re cluster them 
        #     ## together if they are close enough for us to do this
        #     
        #     clust_rows <- NULL
        #     # clust_rows <- narrow_data[,-which(colnames(narrow_data)=="clust_num")]
        #     # colnames(clust_rows)[which(colnames(clust_rows) == "file")] <- "file_loc"
        #     data_to_work_with <- narrow_data[order(narrow_data$sstart),]
        #     if(nrow(data_to_work_with) == 2){
        #       ## If there is just two separate hits we can just use 
        #       ## these as two separate entries into the db
        #       new_row <- data_to_work_with
        #       sub_vector <- c(1:2)
        #       new_row$subject <- paste(data_to_work_with$subject,
        #                                sub_vector, sep = "_")
        #       cols_to_lose <- c("clust_num", "contig")
        #       new_row <- new_row[, -which(colnames(new_row) %in% cols_to_lose)]
        #       colnames(new_row)[9] <- "file_loc"
        #       
        #     }else{
        #       ## If there's more than two we have to compare each of the hits 
        #       ## against each other to see if we can consider them close enough
        #       ## to merge.
        #       combin_mat <- combinat::combn(c(1:nrow(data_to_work_with)),2)
        #       join_rows <- matrix(NA,nrow = nrow(data_to_work_with), ncol = ncol(combin_mat))
        #       ## This is our matrix to store the sequences that do line up 
        #       ## together, can start a run from each of the combinations used 
        #       ## above. 
        #       join_rows_size <- dim(join_rows)[1] * dim(join_rows)[2]
        #       for(q in 1:ncol(combin_mat)){
        #         ## This for loop compares each of the sequences against each 
        #         ## other. Starting n keeps track of how many seqs we can manage
        #         ## to pair with another seq. while end_n is used to mark the
        #         ## column we get the end seq value from. 
        #         if(abs(data_to_work_with$ssend[combin_mat[1,q]] - data_to_work_with$sstart[combin_mat[2,q]]) <= 20000){
        #           starting_col <- q
        #           starting_n <- 1
        #           diff_len <- 1
        #           current_column <- q
        #           if(starting_col == ncol(combin_mat)){
        #             rows_to_be_joined <- c(combin_mat[1,starting_col]:combin_mat[2,starting_col])
        #           }else{
        #             while (diff_len <= 20000 & current_column < ncol(combin_mat)){
        #               current_column <- starting_col + starting_n
        #               diff_len <- abs(data_to_work_with$ssend[combin_mat[2, (current_column-1)]] - 
        #                                 data_to_work_with$sstart[combin_mat[2, current_column]])
        #               starting_n <- starting_n + 1
        #               end_n <- starting_n - 1 + starting_col
        #             }
        #             rows_to_be_joined <- c(combin_mat[1,starting_col]:combin_mat[2,end_n - 1])
        #           }
        #           join_rows[c(rows_to_be_joined[1]:max(rows_to_be_joined)),starting_col] <- rows_to_be_joined
        #         }
        #       }
        #       if(length(which(is.na(join_rows))) == join_rows_size){
        #         ## This occurs when none of the individual seqs are close 
        #         ## enough together to pair up, so we just use the individual
        #         ## seqs as the rows for the db
        #         new_row <- data_to_work_with
        #         sub_vector <- c(1:nrow(data_to_work_with))
        #         new_row$subject <- paste(data_to_work_with$subject,
        #                                  sub_vector, sep = "_")
        #         cols_to_lose <- c("clust_num", "contig")
        #         new_row <- new_row[, -which(colnames(new_row) %in% cols_to_lose)]
        #         colnames(new_row)[9] <- "file_loc"
        #       }else{
        #         rows_to_check <- 1
        #         nrow_checker <- nrow(data_to_work_with)-1
        #         ## From the join matrix the only rows we're actually 
        #         ## interested in are those were we have compared a 
        #         ## new value in the top row against values from the
        #         ## bottom row. This while loop gets those column ids
        #         ## it is a bit misleading with the name rows_to_check 
        #         ## but we are actually getting the important column numbers
        #         while(nrow_checker > 1){
        #           rows_to_check <- c(rows_to_check,
        #                              rows_to_check[length(rows_to_check)] + nrow_checker)
        #           nrow_checker <- nrow_checker - 1
        #         }
        #         na_rose <- NULL
        #         ## This loop below then determines which of these columns has
        #         ## the longest chain in it. This is where we will then start 
        #         ## our analysis of if the chains overlap.
        #         for (rosa in 1:length(rows_to_check)){
        #           nas <- length(which(!is.na(join_rows[,rows_to_check[rosa]])))
        #           na_rose  <- c(na_rose,nas)
        #         }
        #         max_row <- which.max(na_rose)
        #         
        #         tot_removed <- 0
        #         ## Columns that just contain NAs can occur if our new base
        #         ## comparator is not close to another hit. These are therefore
        #         ## useless to compare against so we will get rid of them.
        #         for(p in 1:length(rows_to_check)){
        #           index <- p
        #           if(p > 1){
        #             index <- p - tot_removed
        #           }
        #           no_nas <- length(which(is.na(join_rows[,rows_to_check[index]])))
        #           if(no_nas == nrow(join_rows)){
        #             rows_to_check <- rows_to_check[-index]
        #             na_rose <- na_rose[-index]
        #             tot_removed <- tot_removed + 1
        #           }
        #         }
        #         max_row <- which.max(na_rose)
        #         ## if after geting rid of these values we now only have one 
        #         ## column from which to get values from, we can just merge 
        #         ## these values together and enter them into our db 
        #         if(length(rows_to_check) == 1){
        #           cols_to_keep <- join_rows[,rows_to_check]
        #           rows_of_data <- join_rows[which(!is.na(cols_to_keep)), rows_to_check]
        #           narrow_data <- data_to_work_with[rows_of_data, ]
        #           subject <- paste(narrow_data[1,1],"_",1,sep = "")
        #           pident <- (sum(narrow_data$pident*narrow_data$align))/(sum(narrow_data$align))
        #           align <- sum(narrow_data$align)
        #           if(narrow_data$orientation[1] == "forward"){
        #             qstart <- min(narrow_data$qstart)
        #             qend <- max(narrow_data$qend)
        #             sstart <- min(narrow_data$sstart)
        #             ssend <- max(narrow_data$ssend)
        #           }else{
        #           qstart <- min(narrow_data$qstart)
        #           qend <- max(narrow_data$qend)
        #           sstart <- max(narrow_data$sstart)
        #           ssend <- min(narrow_data$ssend)
        #           }
        #           bitscore <- sum(narrow_data$bitscore)
        #           #ifelse("y" %in% data_set$mefA == TRUE, mefA<-"y", mefA<- "n")
        #           orientation <- narrow_data$orientation[1]
        #           file_loc <- narrow_data$file[1]
        #           clust_rows <- data.frame(subject, pident, align, qstart, qend, sstart,
        #                                    ssend, bitscore, file_loc, orientation)
        #           new_row <- rbind.data.frame(new_row, clust_rows)
        #           
        #         }else{
        #           ## Now we will test if there's any overlap in the values that
        #           ## we have seen string together. We will then keep only the 
        #           ## longest chains that overlap 
        #           cols_to_keep <- NULL
        #           removed_last_row <- NULL
        #           tot_removed <- 0
        #           while(length(rows_to_check) > 1){
        #             current_max_col <- rows_to_check[max_row]
        #             rows_to_check <- rows_to_check[-max_row]
        #             tot_removed <- 0
        #             
        #             for(k in 1:length(rows_to_check)){
        #               max_col_vals <- which(!is.na(join_rows[,current_max_col]))
        #               if(k > 1)
        #                 index <- k - tot_removed
        #               else 
        #                 index <- k
        #               
        #               overlap_test <- join_rows[,rows_to_check[index]] %in% join_rows[max_col_vals, current_max_col]
        #               any_trues <- length(grep(TRUE, overlap_test))
        #               if(any_trues > 0){
        #                 rows_to_check <- rows_to_check[-index]
        #                 na_rose <- na_rose[-index]
        #                 removed_last_row <- "Y"
        #                 tot_removed <- tot_removed + 1
        #               }
        #             }
        #             cols_to_keep <- cbind(cols_to_keep, join_rows[,current_max_col])
        #             na_rose <- na_rose[-max_row]
        #             max_row <- which.max(na_rose)
        #             if(length(rows_to_check) == 1)
        #               cols_to_keep <- cbind(cols_to_keep, join_rows[, rows_to_check])
        #           }
        #           newer_row <- NULL
        #           for(bnbn in 1: ncol(cols_to_keep)){
        #             ## Now we merge each of these longer chains into one entry 
        #             ## for our db 
        #             rows_of_data <- cols_to_keep[!is.na(cols_to_keep[,bnbn]),bnbn]
        #             narrow_data <- data_to_work_with[rows_of_data, ]
        #             subject <- paste(narrow_data[1,1],"_",bnbn,sep = "")
        #             pident <- (sum(narrow_data$pident*narrow_data$align))/(sum(narrow_data$align))
        #             align <- sum(narrow_data$align)
        #             qstart <- min(narrow_data$qstart)
        #             qend <- max(narrow_data$qend)
        #             if(narrow_data$orientation[1] == "forward"){
        #               sstart <- min(narrow_data$sstart)
        #               ssend <- max(narrow_data$ssend)
        #             } else{
        #               sstart <- max(narrow_data$sstart)
        #               ssend <- min(narrow_data$ssend)
        #             }
        #             bitscore <- sum(narrow_data$bitscore)
        #             #ifelse("y" %in% data_set$mefA == TRUE, mefA<-"y", mefA<- "n")
        #             orientation <- narrow_data$orientation[1]
        #             file_loc <- narrow_data$file[1]
        #             clust_rows <- data.frame(subject, pident, align, qstart, qend, sstart,
        #                                      ssend, bitscore, file_loc, orientation)
        #             newer_row <- rbind.data.frame(newer_row, clust_rows)
        #             
        #           }
        #         }
        #         clust_rows <- newer_row
        #         ## Some of the entries will not have formed a chain with 
        #         ## others, here we also enter these into the db 
        #         combed_rows <- cols_to_keep[which(!is.na(cols_to_keep))]
        #         if(length(combed_rows) != nrow(data_to_work_with)){
        #           tot_seqs <- c(1:nrow(data_to_work_with))
        #           tester <- tot_seqs %in% combed_rows
        #           vals_to_put_in <- tot_seqs[grep(FALSE, tester)]
        #           seq_add <- ncol(cols_to_keep)
        #           if(is.null(ncol(cols_to_keep)) == TRUE){
        #             seq_add <- 1
        #           }
        #           
        #           for(r in 1:length(vals_to_put_in)){
        #             new_row_add <- data_to_work_with[vals_to_put_in[r],]
        #             sub_vector <- seq_add + r
        #             new_row_add$subject <- paste(new_row_add$subject,
        #                                          sub_vector, sep = "_")
        #             cols_to_lose <- c("clust_num", "contig")
        #             new_row <- new_row[, -which(colnames(new_row) %in% cols_to_lose)]
        #             colnames(new_row_add)[9] <- "file_loc"
        #             newer_row <- rbind.data.frame(newer_row, new_row_add)
        #           }
        #           clust_rows <- newer_row
        #         }
        #       }
        #     }
        #   }
        #   new_row <- rbind(new_row, clust_rows)
        # }
        
      }else{
        contig_nums <- dplyr::count(data_set, contig)
        for(k in 1:nrow(contig_nums)){
          current_contig <- contig_nums[k, 1]
          data_set_contig <- data_set[data_set$contig == current_contig, ]
          
          subject <- paste(data_set_contig[1,1], "_", as.character(k), sep = "")                                             ###outcome from a column
          pident <- (sum(data_set_contig$pident*data_set_contig$align))/(sum(data_set_contig$align))
          align <- sum(data_set_contig$align)
          if(data_set_contig$orientation[1] == "forward"){
                      qstart <- min(data_set_contig$qstart)
                      qend <- max(data_set_contig$qend)
                      sstart <- min(data_set_contig$sstart)
                      ssend <- max(data_set_contig$ssend)
                    }else{
                      qstart <- min(data_set_contig$qstart)
                      qend <- max(data_set_contig$qend)
                      sstart <- max(data_set_contig$sstart)
                      ssend <- min(data_set_contig$ssend)
                    }
          bitscore <- sum(data_set_contig$bitscore)
          #ifelse("y" %in% data_set$mefA == TRUE, mefA<-"y", mefA<- "n")
          orientation <- data_set_contig$orientation[1]
          file_loc <- data_set_contig$file[1]
          clust_rows <- data.frame(subject, pident, align, qstart, qend, sstart,
                                             ssend, bitscore, file_loc, orientation)
                    
          new_row <- rbind.data.frame(new_row, clust_rows)
    
        }
        
        
        
        new_row <- duplicated_hit_checker(new_row)
        
        
          
        # if(data_set$orientation[1] =="forward"){
        #   subject <- data_set[1,1]                                             ###outcome from a column
        #   pident <- (sum(data_set$pident*data_set$align))/(sum(data_set$align))
        #   align <- sum(data_set$align)
        #   qstart <-min(data_set$qstart)
        #   qend <- max(data_set$qend)
        #   sstart <- min(data_set$sstart)
        #   ssend <- max(data_set$ssend)
        #   bitscore <- sum(data_set$bitscore)
        #   #ifelse("y" %in% data_set$mefA == TRUE, mefA<-"y", mefA<- "n")
        #   orientation <- data_set$orientation[1]
        #   file_loc <- data_set$file[1]
        #   new_row <- data.frame(subject,pident,align,qstart,qend,sstart,ssend
        #                         ,bitscore, file_loc,orientation)
        #   if(abs(new_row$sstart - new_row$ssend) >= 40000){
        #     
        #     new_row <- NULL
        #     # clust_rows <- narrow_data[,-which(colnames(narrow_data)=="clust_num")]
        #     # colnames(clust_rows)[which(colnames(clust_rows) == "file")] <- "file_loc"
        #     data_to_work_with <- data_set[order(data_set$sstart),]
        #     if(nrow(data_to_work_with) == 2){
        #       new_row <- data_to_work_with
        #       sub_vector <- c(1:2)
        #       new_row$subject <- paste(data_to_work_with$subject,
        #                                sub_vector, sep = "_")
        #       new_row <- new_row[, -which(colnames(new_row) == "clust_num")]
        #       colnames(new_row)[9] <- "file_loc"
        #       
        #     }else{
        #       combin_mat <- combinat::combn(c(1:nrow(data_to_work_with)),2)
        #       join_rows <- matrix(NA,nrow = nrow(data_to_work_with), ncol = ncol(combin_mat))
        #       join_rows_size <- dim(join_rows)[1] * dim(join_rows)[2]
        #       for(q in 1:ncol(combin_mat)){
        #         if(abs(data_to_work_with$ssend[combin_mat[1,q]] - data_to_work_with$sstart[combin_mat[2,q]]) <= 20000){
        #           starting_col <- q
        #           starting_n <- 1
        #           diff_len <- 1
        #           current_column <- q
        #           if(starting_col == ncol(combin_mat)){
        #             rows_to_be_joined <- c(combin_mat[1,starting_col]:combin_mat[2,starting_col])
        #           }else{
        #             while (diff_len <= 20000 & current_column < ncol(combin_mat)){
        #               current_column <- starting_col + starting_n
        #               diff_len <- abs(data_to_work_with$ssend[combin_mat[2, (current_column-1)]] - 
        #                                 data_to_work_with$sstart[combin_mat[2, current_column]])
        #               starting_n <- starting_n + 1
        #               end_n <- starting_n - 1 + starting_col
        #             }
        #             rows_to_be_joined <- c(combin_mat[1,starting_col]:combin_mat[2,end_n - 1])
        #           }
        #           join_rows[c(rows_to_be_joined[1]:max(rows_to_be_joined)),starting_col] <- rows_to_be_joined
        #         }
        #       }
        #       if(length(which(is.na(join_rows))) == join_rows_size){
        #         new_row <- data_to_work_with
        #         sub_vector <- c(1:nrow(data_to_work_with))
        #         new_row$subject <- paste(data_to_work_with$subject,
        #                                  sub_vector, sep = "_")
        #         new_row <- new_row[, -which(colnames(new_row) == "clust_num")]
        #         colnames(new_row)[9] <- "file_loc"
        #       }else{
        #         rows_to_check <- 1
        #         nrow_checker <- nrow(data_to_work_with)-1
        #         while(nrow_checker > 1){
        #           rows_to_check <- c(rows_to_check,
        #                              rows_to_check[length(rows_to_check)] + nrow_checker)
        #           nrow_checker <- nrow_checker - 1
        #         }
        #         na_rose <- NULL
        #         for (rosa in 1:length(rows_to_check)){
        #           nas <- length(which(!is.na(join_rows[,rows_to_check[rosa]])))
        #           na_rose  <- c(na_rose,nas)
        #         }
        #         max_row <- which.max(na_rose)
        #         
        #         tot_removed <- 0
        #         for(p in 1:length(rows_to_check)){
        #           index <- p
        #           if(p > 1){
        #             index <- p - tot_removed
        #           }
        #           no_nas <- length(which(is.na(join_rows[,rows_to_check[index]])))
        #           if(no_nas == nrow(join_rows)){
        #             rows_to_check <- rows_to_check[-index]
        #             na_rose <- na_rose[-index]
        #             tot_removed <- tot_removed + 1
        #           }
        #         }
        #         max_row <- which.max(na_rose)
        #         if(length(rows_to_check) == 1){
        #           cols_to_keep <- join_rows[,rows_to_check]
        #           rows_of_data <- join_rows[which(!is.na(cols_to_keep)), rows_to_check]
        #           narrow_data <- data_to_work_with[rows_of_data, ]
        #           subject <- paste(narrow_data[1,1],"_",1,sep = "")
        #           pident <- (sum(narrow_data$pident*narrow_data$align))/(sum(narrow_data$align))
        #           align <- sum(narrow_data$align)
        #           if(narrow_data$orientation[1] == "forward"){
        #             qstart <- min(narrow_data$qstart)
        #             qend <- max(narrow_data$qend)
        #             sstart <- min(narrow_data$sstart)
        #             ssend <- max(narrow_data$ssend)
        #           }else{
        #             qstart <- min(narrow_data$qstart)
        #             qend <- max(narrow_data$qend)
        #             sstart <- max(narrow_data$sstart)
        #             ssend <- min(narrow_data$ssend)
        #           }
        #           bitscore <- sum(narrow_data$bitscore)
        #           #ifelse("y" %in% data_set$mefA == TRUE, mefA<-"y", mefA<- "n")
        #           orientation <- narrow_data$orientation[1]
        #           file_loc <- narrow_data$file[1]
        #           clust_rows <- data.frame(subject, pident, align, qstart, qend, sstart,
        #                                    ssend, bitscore, file_loc, orientation)
        #           new_row <- rbind.data.frame(new_row, clust_rows)
        #           
        #         }else{
        #           
        #           cols_to_keep <- NULL
        #           removed_last_row <- NULL
        #           tot_removed <- 0
        #           while(length(rows_to_check) > 1){
        #             current_max_col <- rows_to_check[max_row]
        #             rows_to_check <- rows_to_check[-max_row]
        #             tot_removed <- 0
        #             for(k in 1:length(rows_to_check)){
        #               max_col_vals <- which(!is.na(join_rows[,current_max_col]))
        #               if(k > 1)
        #                 index <- k - tot_removed
        #               else 
        #                 index <- k
        #               
        #               overlap_test <- join_rows[,rows_to_check[index]] %in% join_rows[max_col_vals, current_max_col]
        #               any_trues <- length(grep(TRUE, overlap_test))
        #               if(any_trues > 0){
        #                 rows_to_check <- rows_to_check[-index]
        #                 na_rose <- na_rose[-index]
        #                 removed_last_row <- "Y"
        #                 tot_removed <- tot_removed + 1
        #               }
        #             }
        #             cols_to_keep <- cbind(cols_to_keep, join_rows[,current_max_col])
        #             na_rose <- na_rose[-max_row]
        #             max_row <- which.max(na_rose)
        #             if(length(rows_to_check) == 1)
        #               cols_to_keep <- cbind(cols_to_keep, join_rows[, rows_to_check])
        #           }
        #           for(bnbn in 1: ncol(cols_to_keep)){
        #             rows_of_data <- cols_to_keep[!is.na(cols_to_keep[,bnbn]),bnbn]
        #             narrow_data <- data_to_work_with[rows_of_data, ]
        #             subject <- paste(narrow_data[1,1],"_",bnbn,sep = "")
        #             pident <- (sum(narrow_data$pident*narrow_data$align))/(sum(narrow_data$align))
        #             align <- sum(narrow_data$align)
        #             if(narrow_data$orientation[1] == "forward"){
        #               qstart <- min(narrow_data$qstart)
        #               qend <- max(narrow_data$qend)
        #               sstart <- min(narrow_data$sstart)
        #               ssend <- max(narrow_data$ssend)
        #             }else{
        #               qstart <- min(narrow_data$qstart)
        #               qend <- max(narrow_data$qend)
        #               sstart <- max(narrow_data$sstart)
        #               ssend <- min(narrow_data$ssend)
        #             }
        #             bitscore <- sum(narrow_data$bitscore)
        #             #ifelse("y" %in% data_set$mefA == TRUE, mefA<-"y", mefA<- "n")
        #             orientation <- narrow_data$orientation[1]
        #             file_loc <- narrow_data$file[1]
        #             clust_rows <- data.frame(subject, pident, align, qstart, qend, sstart,
        #                                      ssend, bitscore, file_loc, orientation)
        #             new_row <- rbind.data.frame(new_row, clust_rows)
        #             
        #           }
        #         }
        #         ## checking for rows that weren't combined at all 
        #         combed_rows <- cols_to_keep[which(!is.na(cols_to_keep))]
        #         if(length(combed_rows) != nrow(data_to_work_with)){
        #           tot_seqs <- c(1:nrow(data_to_work_with))
        #           tester <- tot_seqs %in% combed_rows
        #           vals_to_put_in <- tot_seqs[grep(FALSE, tester)]
        #           seq_add <- ncol(cols_to_keep)
        #           if(is.null(ncol(cols_to_keep)) == TRUE){
        #             seq_add <- 1
        #           }
        #           for(r in 1:length(vals_to_put_in)){
        #             new_row_add <- data_to_work_with[vals_to_put_in[r],]
        #             sub_vector <- seq_add + r
        #             new_row_add$subject <- paste(new_row_add$subject,
        #                                          sub_vector, sep = "_")
        #             new_row_add <- new_row_add[, -which(colnames(new_row_add) == "clust_num")]
        #             colnames(new_row_add)[9] <- "file_loc"
        #             new_row <- rbind.data.frame(new_row, new_row_add)
        #           }
        #         }
        #         print(paste0("This has some spread out hits as well: ", new_row$subject[1]))
        #       }
        #     }
        #   }
        #   
        #   #if(abs(new_row$sstart - new_row$ssend) >= 40000){
        #   # new_row <- data_set[, -which(colnames(data_set) == "clust_num")]
        #   #colnames(new_row)[9] <- "file_loc"
        #   
        # }
        # if(data_set$orientation[1]=="reverse"){
        #   subject <- data_set[1,1]
        #   pident <- (sum(data_set$pident*data_set$align))/(sum(data_set$align))
        #   align <- sum(data_set$align)
        #   qstart <- min(data_set$qstart)
        #   qend <- max(data_set$qend)
        #   sstart <- max(data_set$sstart)
        #   ssend <- min(data_set$ssend)
        #   bitscore <- sum(data_set$bitscore)
        #   #ifelse("y" %in% data_set$mefA == TRUE, mefA<-"y", mefA<- "n")
        #   orientation <- data_set$orientation[1]
        #   file_loc <- data_set$file[1]
        #   new_row <- data.frame(subject,pident,align,qstart,qend,sstart,ssend
        #                         ,bitscore, file_loc,orientation)
        #   if(abs(new_row$sstart - new_row$ssend) >= 40000){
        #     
        #     new_row <- NULL
        #     # clust_rows <- narrow_data[,-which(colnames(narrow_data)=="clust_num")]
        #     # colnames(clust_rows)[which(colnames(clust_rows) == "file")] <- "file_loc"
        #     data_to_work_with <- data_set[order(data_set$sstart),]
        #     if(nrow(data_to_work_with) == 2){
        #       new_row <- data_to_work_with
        #       sub_vector <- c(1:2)
        #       new_row$subject <- paste(data_to_work_with$subject,
        #                                sub_vector, sep = "_")
        #       new_row <- new_row[, -which(colnames(new_row) == "clust_num")]
        #       colnames(new_row)[9] <- "file_loc"
        #       
        #     }else{
        #       combin_mat <- combinat::combn(c(1:nrow(data_to_work_with)),2)
        #       join_rows <- matrix(NA,nrow = nrow(data_to_work_with), ncol = ncol(combin_mat))
        #       join_rows_size <- dim(join_rows)[1] * dim(join_rows)[2]
        #       for(q in 1:ncol(combin_mat)){
        #         if(abs(data_to_work_with$ssend[combin_mat[1,q]] - data_to_work_with$sstart[combin_mat[2,q]]) <= 20000){
        #           starting_col <- q
        #           starting_n <- 1
        #           diff_len <- 1
        #           current_column <- q
        #           if(starting_col == ncol(combin_mat)){
        #             rows_to_be_joined <- c(combin_mat[1,starting_col]:combin_mat[2,starting_col])
        #           }else{
        #             while (diff_len <= 20000 & current_column < ncol(combin_mat)){
        #               current_column <- starting_col + starting_n
        #               diff_len <- abs(data_to_work_with$ssend[combin_mat[2, (current_column-1)]] - 
        #                                 data_to_work_with$sstart[combin_mat[2, current_column]])
        #               starting_n <- starting_n + 1
        #               end_n <- starting_n - 1 + starting_col
        #             }
        #             rows_to_be_joined <- c(combin_mat[1,starting_col]:combin_mat[2,end_n - 1])
        #           }
        #           join_rows[c(rows_to_be_joined[1]:max(rows_to_be_joined)),starting_col] <- rows_to_be_joined
        #           
        #         }
        #       }
        #       if(length(which(is.na(join_rows))) == join_rows_size){
        #         new_row <- data_to_work_with
        #         sub_vector <- c(1:nrow(data_to_work_with))
        #         new_row$subject <- paste(data_to_work_with$subject,
        #                                  sub_vector, sep = "_")
        #         new_row <- new_row[, -which(colnames(new_row) == "clust_num")]
        #         colnames(new_row)[9] <- "file_loc"
        #       }else{
        #         rows_to_check <- 1
        #         nrow_checker <- nrow(data_to_work_with)-1
        #         while(nrow_checker > 1){
        #           rows_to_check <- c(rows_to_check,
        #                              rows_to_check[length(rows_to_check)] + nrow_checker)
        #           nrow_checker <- nrow_checker - 1
        #         }
        #         na_rose <- NULL
        #         for (rosa in 1:length(rows_to_check)){
        #           nas <- length(which(!is.na(join_rows[,rows_to_check[rosa]])))
        #           na_rose  <- c(na_rose,nas)
        #         }
        #         max_row <- which.max(na_rose)
        #         
        #         tot_removed <- 0
        #         for(p in 1:length(rows_to_check)){
        #           index <- p
        #           if(p > 1){
        #             index <- p - tot_removed
        #           }
        #           no_nas <- length(which(is.na(join_rows[,rows_to_check[index]])))
        #           if(no_nas == nrow(join_rows)){
        #             rows_to_check <- rows_to_check[-index]
        #             na_rose <- na_rose[-index]
        #             tot_removed <- tot_removed + 1
        #           }
        #         }
        #         max_row <- which.max(na_rose)
        #         if(length(rows_to_check) == 1){
        #           cols_to_keep <- join_rows[,rows_to_check]
        #           rows_of_data <- join_rows[which(!is.na(cols_to_keep)), rows_to_check]
        #           narrow_data <- data_to_work_with[rows_of_data, ]
        #           subject <- paste(narrow_data[1,1],"_",1,sep = "")
        #           pident <- (sum(narrow_data$pident*narrow_data$align))/(sum(narrow_data$align))
        #           align <- sum(narrow_data$align)
        #           if(narrow_data$orientation[1] == "forward"){
        #             qstart <- min(narrow_data$qstart)
        #             qend <- max(narrow_data$qend)
        #             sstart <- min(narrow_data$sstart)
        #             ssend <- max(narrow_data$ssend)
        #           }else{
        #             qstart <- min(narrow_data$qstart)
        #             qend <- max(narrow_data$qend)
        #             sstart <- max(narrow_data$sstart)
        #             ssend <- min(narrow_data$ssend)
        #           }
        #           bitscore <- sum(narrow_data$bitscore)
        #           #ifelse("y" %in% data_set$mefA == TRUE, mefA<-"y", mefA<- "n")
        #           orientation <- narrow_data$orientation[1]
        #           file_loc <- narrow_data$file[1]
        #           clust_rows <- data.frame(subject, pident, align, qstart, qend, sstart,
        #                                    ssend, bitscore, file_loc, orientation)
        #           new_row <- rbind.data.frame(new_row, clust_rows)
        #           
        #         }else{
        #           cols_to_keep <- NULL
        #           removed_last_row <- NULL
        #           tot_removed <- 0
        #           while(length(rows_to_check) > 1){
        #             
        #             tot_removed <- 0
        #             current_max_col <- rows_to_check[max_row]
        #             rows_to_check <- rows_to_check[-max_row]
        #             for(k in 1:length(rows_to_check)){
        #               max_col_vals <- which(!is.na(join_rows[,current_max_col]))
        #               if(k > 1)
        #                 index <- k - tot_removed
        #               else 
        #                 index <- k
        #               
        #               overlap_test <- join_rows[,rows_to_check[index]] %in% join_rows[max_col_vals, current_max_col]
        #               any_trues <- length(grep(TRUE, overlap_test))
        #               if(any_trues > 0){
        #                 rows_to_check <- rows_to_check[-index]
        #                 na_rose <- na_rose[-index]
        #                 removed_last_row <- "Y"
        #                 tot_removed <- tot_removed + 1
        #               }
        #             }
        #             cols_to_keep <- cbind(cols_to_keep, join_rows[,current_max_col])
        #             na_rose <- na_rose[-max_row]
        #             max_row <- which.max(na_rose)
        #             if(length(rows_to_check) == 1){
        #               cols_to_keep <- cbind(cols_to_keep, join_rows[, rows_to_check])
        #             }
        #           }
        #           for(bnbn in 1: ncol(cols_to_keep)){
        #             rows_of_data <- cols_to_keep[!is.na(cols_to_keep[,bnbn]),bnbn]
        #             narrow_data <- data_to_work_with[rows_of_data, ]
        #             subject <- paste(narrow_data[1,1],"_",bnbn,sep = "")
        #             pident <- (sum(narrow_data$pident*narrow_data$align))/(sum(narrow_data$align))
        #             align <- sum(narrow_data$align)
        #             if(narrow_data$orientation[1] == "forward"){
        #               qstart <- min(narrow_data$qstart)
        #               qend <- max(narrow_data$qend)
        #               sstart <- min(narrow_data$sstart)
        #               ssend <- max(narrow_data$ssend)
        #             }else{
        #               qstart <- min(narrow_data$qstart)
        #               qend <- max(narrow_data$qend)
        #               sstart <- max(narrow_data$sstart)
        #               ssend <- min(narrow_data$ssend)
        #             }
        #             bitscore <- sum(narrow_data$bitscore)
        #             #ifelse("y" %in% data_set$mefA == TRUE, mefA<-"y", mefA<- "n")
        #             orientation <- narrow_data$orientation[1]
        #             file_loc <- narrow_data$file[1]
        #             clust_rows <- data.frame(subject, pident, align, qstart, qend, sstart,
        #                                      ssend, bitscore, file_loc, orientation)
        #             new_row <- rbind.data.frame(new_row, clust_rows)
        #             
        #           }
        #         }
        #         ## checking for rows that weren't combined at all 
        #         combed_rows <- cols_to_keep[which(!is.na(cols_to_keep))]
        #         if(length(combed_rows) != nrow(data_to_work_with)){
        #           tot_seqs <- c(1:nrow(data_to_work_with))
        #           tester <- tot_seqs %in% combed_rows
        #           vals_to_put_in <- tot_seqs[grep(FALSE, tester)]
        #           seq_add <- ncol(cols_to_keep)
        #           if(is.null(ncol(cols_to_keep)) == TRUE){
        #             seq_add <- 1
        #           }
        #           for(r in 1:length(vals_to_put_in)){
        #             new_row_add <- data_to_work_with[vals_to_put_in[r],]
        #             sub_vector <- seq_add + r
        #             new_row_add$subject <- paste(new_row_add$subject,
        #                                          sub_vector, sep = "_")
        #             new_row_add <- new_row_add[, -which(colnames(new_row_add) == "clust_num")]
        #             colnames(new_row_add)[9] <- "file_loc"
        #             new_row <- rbind.data.frame(new_row, new_row_add)
        #           }
        #         }
        #         print(paste0("This has some spread out hits as well: ", new_row$subject[1]))
        #       }
        #     }
        #   }
        #   
        # }
      }
    }
    
    data_merged <- rbind(data_merged,new_row)
    cat("Done", i , "-", nrow(nest_accession), "\n")
    
  }
  
  
  return(data_merged)
}


###############################################################################
## Now lets use this function on our in_csv ###################################
###############################################################################

merged_csv <- delete_dups_blast(in_csv, contig_dir)

write.csv(merged_csv, out_csv_file_name, row.names = FALSE)
