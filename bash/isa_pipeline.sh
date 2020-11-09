#!/bin/bash

set -e
ref_isolate_gff="$1_reference_isolate_gff.csv"
ref_isolate_fasta="$1_reference_isolate_fasta.csv"
fasta_list="$1_fastas_list"

echo "This is the reference isolate gff $ref_isolate_gff"
echo "This is the reference isolate fasta $ref_isolate_fasta"
echo "This is the fasta list $fasta_list"

#echo "This is the reference gff: $1"
#echo "This is the location of the collection of fastas: $2"
echo "This is the MGE/gene fasta: $2"
#echo "This is the gubbins recombination gff: $4"
echo "This is the alignment length cutoff: $3"
#echo "This is the node labelled tree to use: $6"
echo "This is the outfolder name: $4"
echo "This is the prefix: $5"
echo "This is where the gubbins_res are stored $6"
echo "Flank length to extract: $7"
echo "Go through Act compos (yes/loc_of_act_compos): $8"
echo "Just blast res (yes/no): $9"
echo "Reference db loc : ${10}"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
parentdir="$(dirname "$DIR")"
pythondir="${parentdir}/python/"
perldir="${parentdir}/perl/"
rdir="${parentdir}/R/"


echo "This is number 9:$9"
if [ $9 == "no" ]
then
  echo "This is working here"
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
fi
#read -p "Are you sure? " -n 1 -r
#echo    # (optional) move to a new line
#if [[ $REPLY =~ ^[Nn]$ ]]
if [ ! -f $ref_isolate_gff ] || [ ! -f $ref_isolate_fasta ] || [ ! -f $fasta_list ] || [ ! -f $2 ]
then
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "Input files don't exist"
    echo "exiting script now"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  else
    if [ $9 == "no" ] || [ $9 == "NO" ] || [ $9 == "No" ]
    then
      Lines=$(cat $fasta_list | wc -l)

      counter=0
      if [ $8 == "yes" ]
      then

        if [ -d ./act_compos ]
        then
          rm -r ./act_compos
          mkdir ./act_compos
        else
          mkdir ./act_compos
        fi

          python "${pythondir}running_act_comparisons.py" \
          --csv $ref_isolate_fasta --perl_dir $perldir --act_dir ./act_compos/

          act_compos_loc="./act_compos/"
      else
        act_compos_loc=$8
      fi
      if [ -d ./tmp_dna_dir ]
      then
        echo "Using already made tmp_dna_dir"
        cd ./tmp_dna_dir
        if [ ! -f output.mfa ]
        then
          cat *.dna > output.mfa
        fi
      else
        mkdir ./tmp_dna_dir
        cd ./tmp_dna_dir
        less "../${fasta_list}" | while read line; do \
          perl "${perldir}converting_velvet_contigs_to_dna.pl" $line;

          done
        if [ ! -f output.mfa ]
        then
          cat *.dna > output.mfa
        fi

      fi




      if [ ! -f temp_blast_db.nin ] && [ ! -f temp_blast_db.nal ]
      then
      makeblastdb -dbtype nucl -out temp_blast_db -max_file_sz 2GB \
      -in output.mfa
      fi
      #rm output.mfa

      blastn -db temp_blast_db -query "../$2" -outfmt 10 -num_alignments 1000000 -evalue 0.05 -out "$5_tmp_blast_csv.csv"

      cd ../




      if [ ! -d $4 ]
      then
          mkdir $4
      fi

      if [ ! -d ./contig_bounds ]
      then
      exec 3<&0
      exec 0<$fasta_list
      while read line
       do
        grep -n ">" $line > "grep_output_temp"

        python "${pythondir}contig_bound_check.py" \
        --grep_file ./grep_output_temp --fasta_file $line

         counter=$(expr $counter + 1)

         if [ $counter != $Lines ]
         then
         lines_left=$(expr $Lines - $counter)
         echo -ne "This many files left for contigs check: $lines_left "\\r
         fi

       done
      echo "Done contig bounds checker"

      if [ ! -d "./contig_bounds" ]
      then
           mkdir ./contig_bounds
      fi


      mv *contig_bounds.csv ./contig_bounds
      fi


    Rscript --vanilla "${rdir}merging_blast_hits.R" \
    "./tmp_dna_dir/$5_tmp_blast_csv.csv" ./contig_bounds "$4/$5_merged_blast_file"

    Fazza=$fasta_list
    Lines=$(cat $fasta_list | wc -l)
    counter=0

    ## Now we run the library creation step
    python "${pythondir}library_creator.py" \
    --hit_csv "$4/$5_merged_blast_file" --reference_csv $ref_isolate_gff --align_cutoff $3 \
    --act_loc "${act_compos_loc}referoo.fasta." --contig_loc ./contig_bounds/ --output "$4/$5" \
    --fasta_csv $ref_isolate_fasta

    ## Now for the hit allocation step
    python "${pythondir}hit_allocator.py" \
    --blast_csv "$4/$5_merged_blast_file" --lib_csv "$4/$5_library.csv" --reference_csv $ref_isolate_gff \
    --fasta_csv $ref_isolate_fasta --act_loc "${act_compos_loc}referoo.fasta." --contig_loc ./contig_bounds/ \
    --output "$4/$5" --align $3

    ## Now for the edge list
      python "${pythondir}edge_list.py" --hit_csv "$4/$5_hits_df.csv" --out_name "$4/$5_edge_list.tsv"

    ## Now lets find the reccy hits for the isolates
      python "${pythondir}reccy_detector.py" --gubbins_res $6 --hit_locs_csv "$4/$5_hits_df.csv" \
      --out_name "$4/$5" --contig_bounds ./contig_bounds/

  fi

    python "${pythondir}node_mutation.py" --gubbins_res $6 --hit_locs_csv "$4/$5_hits_df.csv" \
    --out_name $5

     if [ ! -d "$4/$5_flanks" ]
      then
        mkdir "$4/$5_flanks"
      else
      rm -r "$4/$5_flanks"
      mkdir "$4/$5_flanks"
      fi


    ## Now to find out the blast locations
    python "${pythondir}blast_hits.py" --gubbins_res $6 --reccy_hits "$4/$5_reccy_hits.csv" \
    --hit_csv "$4/$5_hits_df.csv" --act_compos "${act_compos_loc}referoo.fasta." --flank_length $7 --dna_dir "./tmp_dna_dir/" \
    --out_dir "$4/$5_flanks" --out_name "$4/$5_flanks_extracted.csv" --contig_bounds ./contig_bounds/ \
    --proper_hits "$4/$5_proper_hits.csv"



    ls -d "$PWD/$4/$5_flanks/"*_whole_blast_seq.fasta* > "$4/$5_isolate_blast_list"
    ls -d "$PWD/$4/$5_flanks/"*_before_flank.fasta* > "$4/$5_before_flanks_list"
    ls -d "$PWD/$4/$5_flanks/"*_after_flank.fasta* > "$4/$5_after_flanks_list"


  blast_length=$(ls -l "$4/$5_flanks/"*_whole_blast_seq.fasta | wc -l)

  el_nuevo_fazza="$4/$5_isolate_blast_list"
  before_flanks_fazza="$4/$5_before_flanks_list"
  after_flanks_fazza="$4/$5_after_flanks_list"

  if [ ! -d "$4/blast_results" ]
    then
        mkdir "$4/blast_results"
  else
      rm -r "$4/blast_results"
      mkdir "$4/blast_results"
    fi

  if [ ! -d "$4/before_flank_blast_res" ]
    then
        mkdir "$4/before_flank_blast_res"
  else
    rm -r "$4/before_flank_blast_res"
    mkdir "$4/before_flank_blast_res"

    fi

  if [ ! -d "$4/after_flank_blast_res" ]
    then
        mkdir "$4/after_flank_blast_res"

  else
    rm -r "$4/after_flank_blast_res"
    mkdir "$4/after_flank_blast_res"
    fi


  counter_2=0

  echo "This many isolates to BLAST: $blast_length"



    DB_FILE=${10}.nin
    db_loc=${10}
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
      --results_csv "$4/blast_results/$namo".csv --out_dir "$4/blast_results/"


    done < $el_nuevo_fazza

    ## Just doing the before flanks now

    while read coffee
    do
      namo=$(basename $coffee)
      echo "Starting before BLAST on this isolate: $namo"
      blastn -db $db_loc -query $coffee -out "$4/before_flank_blast_res/$namo".csv \
      -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"

      python "${pythondir}blast_local_results_interpreter.py" \
      --results_csv "$4/before_flank_blast_res/$namo".csv --out_dir "$4/before_flank_blast_res/"


    done < $before_flanks_fazza


    ## Now for the after flanks

    while read coffee
    do
      namo=$(basename $coffee)
      echo "Starting after BLAST on this isolate: $namo"
      blastn -db $db_loc -query $coffee -out "$4/after_flank_blast_res/$namo".csv \
      -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"

      python "${pythondir}blast_local_results_interpreter.py" \
      --results_csv "$4/after_flank_blast_res/$namo".csv --out_dir "$4/after_flank_blast_res/"



    done < $after_flanks_fazza



  echo "Starting species list formation"

  ls -d "$PWD/$4/blast_results/"*_species_list.csv > "$4/blast_results/$5_species_list_list"
  ls -d "$PWD/$4/before_flank_blast_res/"*_species_list.csv > "$4/before_flank_blast_res/$5_species_list_before"
  ls -d "$PWD/$4/after_flank_blast_res/"*_species_list.csv > "$4/after_flank_blast_res/$5_species_list_after"

  python "${pythondir}blast_to_info.py" \
  --list_file "$4/blast_results/$5_species_list_list" --hit_locs "$4/$5_hits_df.csv" \
   --out_name "$4/blast_results/$5_species_compo.csv"

  python "${pythondir}blast_to_info.py" \
  --list_file "$4/before_flank_blast_res/$5_species_list_before" --hit_locs "$4/$5_hits_df.csv" \
   --out_name "$4/before_flank_blast_res/$5_species_compo_before.csv"

  python "${pythondir}blast_to_info.py" \
  --list_file "$4/after_flank_blast_res/$5_species_list_after" --hit_locs "$4/$5_hits_df.csv" \
  --out_name "$4/after_flank_blast_res/$5_species_compo_after.csv"



  echo "Done"


fi