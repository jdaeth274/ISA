
pmen_mega_hits <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_mega_updated_lib/pmen_mega_hits_df.csv",
                           stringsAsFactors = FALSE)
pmen_tn916_hits <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_updated_lib/pmen_tn916_hits_df.csv",
                            stringsAsFactors = FALSE)
pmen3_trimethoprim <- read.csv("~/Dropbox/phd/PMEN3/test_pen/pmen3_dhfr_res.csv", stringsAsFactors = FALSE)

pmen3_sulfa <- read.csv("~/Dropbox/phd/PMEN3/test_pen/pmen3_folP_res.csv", stringsAsFactors = FALSE)

pmen9_trimethoprim <- read.csv("~/Dropbox/phd/PMEN3/test_pen/pmen9_dhfr_res2.csv", stringsAsFactors = FALSE,
                               header = FALSE, col.names = c("isolate_id","Resistance"))

pmen9_sulfa <- read.csv("~/Dropbox/phd/PMEN3/test_pen/pmen_folP_res.csv", stringsAsFactors = FALSE)

pmen9_epi_csv <- read.csv("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_14_10_2020.csv",
                          stringsAsFactors = FALSE)

pmen3_epi_csv <- read.csv("~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_22_08_2020.csv",
                          stringsAsFactors = FALSE)

pmen_mega_narrow <- select(pmen_mega_hits, c(id, insert_name)) %>% rename(insert_name_mega = insert_name)
pmen_tn916_narrow <- select(pmen_tn916_hits, c(id, insert_name)) %>% rename(insert_name_tn916 = insert_name)

pmen_trim <- bind_rows(pmen3_trimethoprim, pmen9_trimethoprim) %>% rename(id = isolate_id, trimethoprim__autocolour = Resistance)
pmen_sulfa <- bind_rows(pmen3_sulfa, pmen9_sulfa) %>% rename(id = isolate_id, sulfa__autocolour = Resistance)

pmen9_epi_csv <- left_join(x = pmen9_epi_csv,y = pmen_mega_narrow, by = "id") %>% 
  mutate(insert_name_mega = ifelse(is.na(insert_name_mega), 0, insert_name_mega)) %>% rename(insert_name_mega__autocolour = insert_name_mega) %>% left_join(y = pmen_tn916_narrow, by = "id") %>% 
  mutate(insert_name_tn916 = ifelse(is.na(insert_name_tn916), 0, insert_name_tn916)) %>% rename(insert_name_tn916__autocolour = insert_name_tn916) %>%
  left_join(pmen_sulfa, by = "id") %>% mutate(sulfa__autocolour = ifelse(is.na(sulfa__autocolour), "ND", sulfa__autocolour))
                                                                                             
pmen3_epi_csv <- left_join(x = pmen3_epi_csv,y = pmen_mega_narrow, by = "id") %>% 
  mutate(insert_name_mega = ifelse(is.na(insert_name_mega), 0, insert_name_mega)) %>% rename(insert_name_mega__autocolour = insert_name_mega) %>% left_join(y = pmen_tn916_narrow, by = "id") %>% 
  mutate(insert_name_tn916 = ifelse(is.na(insert_name_tn916), 0, insert_name_tn916)) %>% rename(insert_name_tn916__autocolour = insert_name_tn916) %>%
  left_join(pmen_sulfa, by = "id") %>% mutate(sulfa__autocolour = ifelse(is.na(sulfa__autocolour), "ND", sulfa__autocolour)) %>%
  left_join(pmen_trim, by = "id") %>% mutate(trimethoprim__autocolour = ifelse(is.na(trimethoprim__autocolour), "ND", trimethoprim__autocolour))

write.csv(pmen9_epi_csv, file = "~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_04_11_2020.csv",
          row.names = FALSE, quote = FALSE)

write.csv(pmen3_epi_csv, file = "~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_04_11_2022.csv",
          row.names = FALSE, quote = FALSE)













