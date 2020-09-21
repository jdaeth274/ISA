## gff list creator 

pmen3_gff <- readLines("~/Dropbox/phd/gubbins_project/data/pmen3/pmen3_run_data/pmen3_gff_list")
pmen9_gff <- readLines("~/Dropbox/phd/gubbins_project/data/pmen9/pmen9_run_data/pmen9_gff_list")
pmen9_gff <- sub(".dna_prokka",".fa_prokka", pmen9_gff)

pmen3_ref <- "~/Dropbox/phd/gubbins_project/data/pmen3/pmen3_run_data/pmen3_reference.gff"
pmen9_ref <- "~/Dropbox/phd/gubbins_project/data/pmen9/pmen9_run_data/pmen9_reference.gff"

pmens <- c(pmen3_gff, pmen9_gff)

pmens_ref <- c(rep(pmen3_ref, length(pmen3_gff)), rep(pmen9_ref, length(pmen9_gff)))

pmen_ref_gff <- cbind.data.frame(pmens, pmens_ref)
colnames(pmen_ref_gff) <- c("isolate","reference")

write.csv(pmen_ref_gff, file = "~/Dropbox/phd/insertion_site_analysis/data/reference_isolate_gff.csv", row.names = FALSE)

shortened_csv <- pmen_ref_gff[c(625:650, 725:750),]
write.csv(shortened_csv, file = "~/Dropbox/phd/insertion_site_analysis/data/short_reference_isolate.csv", row.names = FALSE)

pmen9_only <- pmen_ref_gff[710:715,]
write.csv(pmen9_only, file = "~/Dropbox/phd/insertion_site_analysis/data/pmen9_reference_isolate.csv", row.names = FALSE)

## fasta gff creator
pmen3_fasta <- readLines("~/Dropbox/phd/gubbins_project/data/pmen3/pmen3_run_data/pmen3_fastas_list")
pmen9_fasta <- readLines("~/Dropbox/phd/gubbins_project/data/pmen9/pmen9_run_data/pmen9_fastas_list")

pmen3_ref <- "~/Dropbox/phd/gubbins_project/data/pmen3/pmen3_run_data/pmen3_reference.gff"
pmen9_ref <- "~/Dropbox/phd/gubbins_project/data/pmen9/pmen9_run_data/pmen9_reference.gff"

pmens <- c(pmen3_fasta, pmen9_fasta)

pmens_ref <- c(rep(pmen3_ref, length(pmen3_fasta)), rep(pmen9_ref, length(pmen9_fasta)))

pmen_ref_fasta <- cbind.data.frame(pmens, pmens_ref)
colnames(pmen_ref_gff) <- c("isolate","reference")

write.csv(pmen_ref_fasta, file = "~/Dropbox/phd/insertion_site_analysis/data/reference_isolate_fasta.csv", row.names = FALSE)

shortened_csv <- pmen_ref_fasta[c(625:650, 725:750),]
write.csv(shortened_csv, file = "~/Dropbox/phd/insertion_site_analysis/data/short_reference_isolate_fasta.csv", row.names = FALSE)

pmen9_only <- pmen_ref_gff[710:715,]
write.csv(pmen9_only, file = "~/Dropbox/phd/insertion_site_analysis/data/pmen9_reference_isolate.csv", row.names = FALSE)
