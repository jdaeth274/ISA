#!/bin/bash

set -e

if [ $# != 9 ]
then
  echo "This need 9 args you have $#."
  echo "Requires: <run_output_folder/prefix> <tmp_dna_dir> <actcompos_dir/prefix> <outdir> <prefix> <flank length file> <gubbins_res_loc> <strep_ref_db_loc> <contig_bound_location>"
  exit
fi

hit_csv="$1_hits_df.csv"

reccy_csv="$1_reccy_hits.csv"


echo "This is the hit_locs csv: $hit_csv"
echo "This is the reccy_hits csv: $reccy_csv"
echo "This is the dna directory: $2"

echo "These are the actcompos with prefix: $3"
echo "This is the outdir: $4"
echo "This is the prefix: $5"
echo "Flanking region length file: $6"
echo "This is the gubbins_res_loc $7"
echo "This is the reference db base_file $8"
echo "This is the contig bound location $9"


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
parentdir="$(dirname "$DIR")"
pythondir="${parentdir}/python/"
perldir="${parentdir}/perl/"
rdir="${parentdir}/R/"

if [ ! -f $hit_csv ] || [ ! -f $reccy_csv ] || [ ! -f $6 ] || [ ! -d $2 ]
then
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  echo "Missing files check the usage above"
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  exit
else
  if [ ! -d $4 ]
  then
  mkdir $4
  fi

  less $6 | while read line; do \
    current_flanks=$line

    echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
    echo "On flank length: ${current_flanks}"

  if [ ! -d "$4/$5_${current_flanks}_flanks" ]
      then
      mkdir "$4/$5_${current_flanks}_flanks"
    else
      rm -r "$4/$5_${current_flanks}_flanks"
      mkdir "$4/$5_${current_flanks}_flanks"
    fi

  python "${pythondir}blast_hits.py" --gubbins_res $7 --reccy_hits $reccy_csv \
      --hit_csv $hit_csv --act_compos $3 --flank_length $current_flanks --dna_dir $2 \
      --out_dir "$4/$5_${current_flanks}_flanks" --out_name "$4/$5_${current_flanks}_flanks_extracted.csv" --contig_bounds $9 \
    --proper_hits "$4/$5_proper_hits.csv"


    ls -d "$PWD/$4/$5_${current_flanks}_flanks/"*_whole_blast_seq.fasta* > "$4/$5_${current_flanks}_isolate_blast_list"
    ls -d "$PWD/$4/$5_${current_flanks}_flanks/"*_before_flank.fasta* > "$4/$5_${current_flanks}_before_flanks_list"
    ls -d "$PWD/$4/$5_${current_flanks}_flanks/"*_after_flank.fasta* > "$4/$5_${current_flanks}_after_flanks_list"


    blast_length=$(ls -l "$4/$5_${current_flanks}_flanks/"*_whole_blast_seq.fasta | wc -l)

    el_nuevo_fazza="$4/$5_${current_flanks}_isolate_blast_list"
    before_flanks_fazza="$4/$5_${current_flanks}_before_flanks_list"
    after_flanks_fazza="$4/$5_${current_flanks}_after_flanks_list"
    if [ ! -d "$4/${current_flanks}_blast_results" ]
      then
      mkdir "$4/${current_flanks}_blast_results"
    else
      rm -r "$4/${current_flanks}_blast_results"
      mkdir "$4/${current_flanks}_blast_results"
    fi

    if [ ! -d "$4/${current_flanks}_before_flank_blast_res" ]
      then
      mkdir "$4/${current_flanks}_before_flank_blast_res"
    else
      rm -r "$4/${current_flanks}_before_flank_blast_res"
      mkdir "$4/${current_flanks}_before_flank_blast_res"
    fi

    if [ ! -d "$4/${current_flanks}_after_flank_blast_res" ]
      then
      mkdir "$4/${current_flanks}_after_flank_blast_res"
    else
      rm -r "$4/${current_flanks}_after_flank_blast_res"
      mkdir "$4/${current_flanks}_after_flank_blast_res"
    fi

    counter_2=0

    echo "This many isolates to BLAST: $blast_length"


      DB_FILE="$8.nin"
      db_loc=$8
      if [ ! -f $DB_FILE ]
      then
        makeblastdb -dbtype nucl -out strep_ref_dna_db -max_file_sz 2GB \
        -in $8
      fi
      while read coffee
      do
        namo=$(basename $coffee)
        echo "Starting BLAST on this isolate: $namo"
        blastn -db $db_loc -query $coffee -out "$4/${current_flanks}_blast_results/$namo".csv \
        -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"

        python "${pythondir}blast_local_results_interpreter.py" \
        --results_csv "$4/${current_flanks}_blast_results/$namo".csv --out_dir "$4/${current_flanks}_blast_results/"


      done < $el_nuevo_fazza

      ## Just doing the before flanks now

      while read coffee
      do
        namo=$(basename $coffee)
        echo "Starting before BLAST on this isolate: $namo"
        blastn -db $db_loc -query $coffee -out "$4/${current_flanks}_before_flank_blast_res/$namo".csv \
        -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"

        python "${pythondir}blast_local_results_interpreter.py" \
        --results_csv "$4/${current_flanks}_before_flank_blast_res/$namo".csv --out_dir "$4/${current_flanks}_before_flank_blast_res/"


      done < $before_flanks_fazza


      ## Now for the after flanks

      while read coffee
      do
        namo=$(basename $coffee)
        echo "Starting after BLAST on this isolate: $namo"
        blastn -db $db_loc -query $coffee -out "$4/${current_flanks}_after_flank_blast_res/$namo".csv \
        -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"

        python "${pythondir}blast_local_results_interpreter.py" \
        --results_csv "$4/${current_flanks}_after_flank_blast_res/$namo".csv --out_dir "$4/${current_flanks}_after_flank_blast_res/"


      done < $after_flanks_fazza


    ls -d "$PWD/$4/${current_flanks}_blast_results/"*_species_list.csv > "$4/${current_flanks}_blast_results/$5_species_list_list"
    ls -d "$PWD/$4/${current_flanks}_before_flank_blast_res/"*_species_list.csv > "$4/${current_flanks}_before_flank_blast_res/$5_species_list_before"
    ls -d "$PWD/$4/${current_flanks}_after_flank_blast_res/"*_species_list.csv > "$4/${current_flanks}_after_flank_blast_res/$5_species_list_after"

    python "${pythondir}blast_to_info.py" \
    --list_file "$4/${current_flanks}_blast_results/$5_species_list_list" --hit_locs $hit_csv \
     --out_name "$4/${current_flanks}_blast_results/$5_species_compo.csv"

    python "${pythondir}blast_to_info.py" \
    --list_file "$4/${current_flanks}_before_flank_blast_res/$5_species_list_before" --hit_locs $hit_csv \
     --out_name "$4/${current_flanks}_before_flank_blast_res/$5_species_compo_before.csv"

    python "${pythondir}blast_to_info.py" \
    --list_file "$4/${current_flanks}_after_flank_blast_res/$5_species_list_after" --hit_locs $hit_csv \
    --out_name "$4/${current_flanks}_after_flank_blast_res/$5_species_compo_after.csv"



    echo "Done"
  done
  echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
fi
