#!/bin/bash

set -e

echo "This is the folder $1"
echo "This is the prefix $2"

ls -d "$PWD/$4/isolate_seqs_for_blast/"*_whole_blast_seq.fasta* > "$4/$5_isolate_blast_list"
ls -d "$PWD/$4/isolate_seqs_for_blast/"*_before_flank.fasta* > "$4/$5_before_flanks_list"
ls -d "$PWD/$4/isolate_seqs_for_blast/"*_after_flank.fasta* > "$4/$5_after_flanks_list"


  blast_length=$(ls -l "$4/isolate_seqs_for_blast/"*_whole_blast_seq.fasta | wc -l)

  el_nuevo_fazza="$4/$5_isolate_blast_list"
  before_flanks_fazza="$4/$5_before_flanks_list"
  after_flanks_fazza="$4/$5_after_flanks_list"

  if [ ! -d "$4/blast_results" ]
    then
        mkdir "$4/blast_results"
    fi

  if [ ! -d "$4/before_flank_blast_res" ]
    then
        mkdir "$4/before_flank_blast_res"
    fi

  if [ ! -d "$4/after_flank_blast_res" ]
    then
        mkdir "$4/after_flank_blast_res"
    fi


  counter_2=0

  echo "This many isolates to BLAST: $blast_length"

  if [ $9 == "r" ]
  then

    while read coffee
    do
      namo=$(basename $coffee)
      echo "Starting BLAST on this isolate: $namo"
      blastn -db nt -query $coffee -out "$4/blast_results/$namo".csv \
      -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore" \
      -remote

      python "${pythondir}eutils_test.py" "$4/blast_results/$namo".csv \
      "$4/blast_results/"



    done < $el_nuevo_fazza
  else
    DB_FILE=~/Dropbox/phd/strep_reference_collection/dna_fasta_files/strep_ref_database.nin
    db_loc=~/Dropbox/phd/strep_reference_collection/dna_fasta_files/strep_ref_database
    if [ ! -f $DB_FILE ]
    then
      makeblastdb -dbtype nucl -out strep_ref_dna_db -max_file_sz 2GB \
      -in ~/Dropbox/phd/strep_reference_collection/dna_fasta_files/strep_ref_dna_db
    fi
    while read coffee
    do
      namo=$(basename $coffee)
      echo "Starting BLAST on this isolate: $namo"
      blastn -db $db_loc -query $coffee -out "$4/blast_results/$namo".csv \
      -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"

      python "${pythondir}blast_local_results_interpreter.py" \
      "$4/blast_results/$namo".csv "$4/blast_results/"


    done < $el_nuevo_fazza

    ## Just doing the before flanks now

    while read coffee
    do
      namo=$(basename $coffee)
      echo "Starting before BLAST on this isolate: $namo"
      blastn -db $db_loc -query $coffee -out "$4/before_flank_blast_res/$namo".csv \
      -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"

      python "${pythondir}blast_local_results_interpreter.py" \
      "$4/before_flank_blast_res/$namo".csv "$4/before_flank_blast_res/"


    done < $before_flanks_fazza


    ## Now for the after flanks

    while read coffee
    do
      namo=$(basename $coffee)
      echo "Starting after BLAST on this isolate: $namo"
      blastn -db $db_loc -query $coffee -out "$4/after_flank_blast_res/$namo".csv \
      -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"

      python "${pythondir}blast_local_results_interpreter.py" \
      "$4/after_flank_blast_res/$namo".csv "$4/after_flank_blast_res/"



    done < $after_flanks_fazza


  fi

  echo "Starting species list formation"

  ls -d "$PWD/$4/blast_results/"*_species_list.csv > "$4/blast_results/$5_species_list_list"
  ls -d "$PWD/$4/before_flank_blast_res/"*_species_list.csv > "$4/before_flank_blast_res/$5_species_list_before"
  ls -d "$PWD/$4/after_flank_blast_res/"*_species_list.csv > "$4/after_flank_blast_res/$5_species_list_after"

  echo "Done"

  python "${pythondir}blast_results_to_info.py" \
  "$4/blast_results/$5_species_list_list" "$4$pyt_csv" "$4/blast_results/$5_species_compo.csv"

  python "${pythondir}blast_results_to_info.py" \
  "$4/before_flank_blast_res/$5_species_list_before" "$4$pyt_csv" "$4/before_flank_blast_res/$5_species_compo_before.csv"

  python "${pythondir}blast_results_to_info.py" \
  "$4/after_flank_blast_res/$5_species_list_after" "$4$pyt_csv" "$4/after_flank_blast_res/$5_species_compo.csv"
