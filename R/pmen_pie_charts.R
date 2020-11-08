require(dplyr)
require(ggplot2)

load("~/Dropbox/phd/PMEN3/whole_tree_relaxed_gamma_bactdating",
     verbose = TRUE)

pmen3_dated_tree <- central_res$tree


current_pmen3_fasta_list <- readLines("~/Dropbox/phd/gubbins_project/data/pmen3/pmen3_run_data/pmen3_fastas_list")

dated_fasta_list <- NULL
dated_isolates <- NULL
for(k in current_pmen3_fasta_list){
  isolate_id <- sub("\\.[a-z,A-Z].*$","",basename(k))

  if(isolate_id %in% pmen3_dated_tree$tip.label){
    dated_fasta_list <- append(dated_fasta_list, k)
    dated_isolates <- append(dated_isolates, isolate_id)
  }
  
}

pmen3_dated_tree <- drop.tip(pmen3_dated_tree, c("6678_3#1","6678_3#12"))

write.tree(pmen3_dated_tree, file = "~/Dropbox/phd/PMEN3/dated_tree_gubbins/dated_661_tree.nwk")

writeLines(dated_fasta_list, con = "~/Dropbox/phd/PMEN3/dated_tree_gubbins/dated_pmen_fasta_list.txt")

saved_tree <- read.tree("~/Dropbox/phd/PMEN3/dated_tree_gubbins/dated_661_tree.nwk")



###############################################################################
## PMEN3 pie charts ###########################################################
###############################################################################

pmen_3_epi_data <- read.csv("~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_22_08_2020.csv", stringsAsFactors = FALSE)
pmen_9_epi_data <- read.csv("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_22_08_2020.csv", stringsAsFactors = FALSE)

## serotype pie charts
pmen_3_epi_data$Serotype <- pmen_3_epi_data$serotype__autocolour
pmen_9_epi_data$Serotype <- pmen_9_epi_data$serotype__autocolour

pmen3_epi_serotype <- dplyr::count(pmen_3_epi_data, Serotype) %>%
  mutate(under10 = ifelse(n < 10, 5,n)) %>% mutate(Sero = ifelse(under10 == 5,"misc",Serotype)) %>%
  group_by(Sero) %>% summarise(count = sum(n))  %>% 
  arrange(desc(Sero)) %>% mutate(prop = count / sum(pmen3_epi_serotype$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>% rename(Serotype = Sero) %>% mutate(prop_lab = paste(as.character(round(prop)), "%", sep = "")) 

  
  
pmen9_epi_serotype <- dplyr::count(pmen_9_epi_data, Serotype)  %>%
  mutate(under10 = ifelse(n < 10, 5,n)) %>% mutate(Sero = ifelse(under10 == 5,"misc",Serotype)) %>%
  group_by(Sero) %>% summarise(count = sum(n))  %>% 
  arrange(desc(Sero)) %>% mutate(prop = count / sum(pmen3_epi_serotype$count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>% rename(Serotype = Sero) %>% mutate(prop_lab = paste(as.character(round(prop)), "%", sep = "")) 



pmen3_serotype <- ggplot(data = pmen3_epi_serotype, aes(x="",y = prop, fill = Serotype)) +
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) + 
  theme_void() +
  geom_text(aes(y = ypos, label = count), color = "White")+
  theme_void() + theme(legend.text = element_text(size = 12))
pmen3_serotype

pmen9_serotype <- ggplot(data = pmen9_epi_serotype, aes(x="",y = prop, fill = Serotype)) +
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) + 
  theme_void() +
  geom_text(aes(y = ypos, label = count), color = "White")+
  theme_void() + theme(legend.text = element_text(size = 12))
pmen9_serotype


length(pmen3_tree$tip.label)

pmen3_tree <- keep.tip(pmen3_tree, pmen3_dated_tree$tip.label)

pmen3_tree$tip.label <- paste(pmen3_tree$tip.label, ".contigs_velvet.fa", sep = "")
write.tree(pmen3_tree, file = "~/Dropbox/phd/PMEN3/dated_tree_gubbins/starting_gubbins_661_tree.nwk")



