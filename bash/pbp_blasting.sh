#!/bin/bash

set -e

if [[ $# == 0 ]]
then
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  echo "Missing files check the usage below"
  echo "Requires: [-p|--pbp_csv] [-d|--dna_dir] [-a|--act_compos] [-o|--out_dir] [-pr|--prefix] [-f|--flanks_file] [-g|--gubbins_res] [-r|--reference_db] [-c|--contig_bounds] [-gff|--gff_loc_csv]"
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  exit
fi


POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
  -p|--pbp_csv)
  PBP="$2"
  shift
  shift
  ;;
  -d|--dna_dir)
  DNA_DIR="$2"
  shift
  shift
  ;;
  -a|--act_compos)
  ACT_COMPOS="$2"
  shift
  shift
  ;;
  -o|--out_dir)
  OUT_DIR="$2"
  shift
  shift
  ;;
  -pr|--prefix)
  PREFIX="$2"
  shift
  shift
  ;;
  -f|--flanks_file)
  FLANKS="$2"
  shift
  shift
  ;;
  -g|--gubbins_res)
  GUBBINS="$2"
  shift
  shift
  ;;
  -r|--reference_db)
  REF_DB="$2"
  shift
  shift
  ;;
  -c|--contig_bounds)
  CONTIGS="$2"
  shift
  shift
  ;;
  -gff|--gff_loc_csv)
  GFF_LOC="$2"
  shift
  shift
  ;;
esac
done

set -- "S{POSITIONAL[@]}"

echo "This is the pbp profile csv: ${PBP}"
echo "This is the dna directory: ${DNA_DIR}"
echo "These are the actcompos with prefix: ${ACT_COMPOS}"
echo "This is the outdir: ${OUT_DIR}"
echo "This is the prefix: ${PREFIX}"
echo "Flanking region length file: ${FLANKS}"
echo "This is the gubbins_res_loc ${GUBBINS}"
echo "This is the reference db base_file ${REF_DB}"
echo "This is the contig bound location ${CONTIGS}"
echo "This is the reference_gff_loc csv ${GFF_LOC}"



DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
parentdir="$(dirname "$DIR")"
pythondir="${parentdir}/python/"
perldir="${parentdir}/perl/"
rdir="${parentdir}/R/"

if [ ! -f $PBP ] || [ ! -d $DNA_DIR ] || [ ! -f $FLANKS ] || [ ! -d $CONTIGS ]
then
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  echo "Missing files check the usage below"
  echo "Requires: [-p|--pbp_csv] [-d|--dna_dir] [-a|--act_compos] [-o|--out_dir] [-pr|--prefix] [-f|--flanks_file] [-g|--gubbins_res] [-r|--reference_db] [-c|--contig_bounds] [-gff|--gff_loc_csv]"
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  exit
else
  if [ ! -d $OUT_DIR ]
  then
  mkdir $OUT_DIR
  fi

  python "${pythondir}pbp_reccy_finder.py" --gubbins_res $GUBBINS --pbp_profiles $PBP\
   --gff_csv $GFF_LOC --contig_bounds $CONTIGS --out_name "${OUT_DIR}/${PREFIX}"

  hit_csv="${OUT_DIR}/${PREFIX}_hits_df.csv"

  less $FLANKS | while read line; do \
    current_flanks=$line

    echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
    echo "On flank length: ${current_flanks}"

  if [ ! -d "${OUT_DIR}/${PREFIX}_${current_flanks}_flanks" ]
      then
      mkdir "${OUT_DIR}/${PREFIX}_${current_flanks}_flanks"
    else
      rm -r "${OUT_DIR}/${PREFIX}_${current_flanks}_flanks"
      mkdir "${OUT_DIR}/${PREFIX}_${current_flanks}_flanks"
    fi

  python "${pythondir}pbp_blast_hits.py" --gubbins_res $GUBBINS --reccy_hits "${OUT_DIR}/${PREFIX}_reccy_hits.csv" \
      --hit_csv "${OUT_DIR}/${PREFIX}_hits_df.csv" --act_compos "${ACT_COMPOS}/referoo.fasta." --flank_length $current_flanks --dna_dir ${DNA_DIR} \
      --out_dir "${OUT_DIR}/${PREFIX}_${current_flanks}_flanks" --out_name "${OUT_DIR}/${PREFIX}_${current_flanks}_flanks_extracted.csv" --contig_bounds $CONTIGS \
    --gff_csv $GFF_LOC


    ls -d "$PWD/${OUT_DIR}/${PREFIX}_${current_flanks}_flanks/"*_whole_blast_seq.fasta* > "${OUT_DIR}/${PREFIX}_${current_flanks}_isolate_blast_list"
    ls -d "$PWD/${OUT_DIR}/${PREFIX}_${current_flanks}_flanks/"*_before_flank.fasta* > "${OUT_DIR}/${PREFIX}_${current_flanks}_before_flanks_list"
    ls -d "$PWD/${OUT_DIR}/${PREFIX}_${current_flanks}_flanks/"*_after_flank.fasta* > "${OUT_DIR}/${PREFIX}_${current_flanks}_after_flanks_list"
    ls -d "$PWD/${OUT_DIR}/${PREFIX}_${current_flanks}_flanks/"*_pbp_sequence.fasta* > "${OUT_DIR}/${PREFIX}_${current_flanks}_pbp_seq_list"

    blast_length=$(ls -l "${OUT_DIR}/${PREFIX}_${current_flanks}_flanks/"*_whole_blast_seq.fasta | wc -l)

    el_nuevo_fazza="${OUT_DIR}/${PREFIX}_${current_flanks}_isolate_blast_list"
    before_flanks_fazza="${OUT_DIR}/${PREFIX}_${current_flanks}_before_flanks_list"
    after_flanks_fazza="${OUT_DIR}/${PREFIX}_${current_flanks}_after_flanks_list"
    pbp_seqs="${OUT_DIR}/${PREFIX}_${current_flanks}_pbp_seq_list"
    if [ ! -d "${OUT_DIR}/${current_flanks}_blast_results" ]
      then
      mkdir "${OUT_DIR}/${current_flanks}_blast_results"
    else
      rm -r "${OUT_DIR}/${current_flanks}_blast_results"
      mkdir "${OUT_DIR}/${current_flanks}_blast_results"
    fi

    if [ ! -d "${OUT_DIR}/${current_flanks}_before_flank_blast_res" ]
      then
      mkdir "${OUT_DIR}/${current_flanks}_before_flank_blast_res"
    else
      rm -r "${OUT_DIR}/${current_flanks}_before_flank_blast_res"
      mkdir "${OUT_DIR}/${current_flanks}_before_flank_blast_res"
    fi

    if [ ! -d "${OUT_DIR}/${current_flanks}_after_flank_blast_res" ]
      then
      mkdir "${OUT_DIR}/${current_flanks}_after_flank_blast_res"
    else
      rm -r "${OUT_DIR}/${current_flanks}_after_flank_blast_res"
      mkdir "${OUT_DIR}/${current_flanks}_after_flank_blast_res"
    fi

    if [ ! -d "${OUT_DIR}/${current_flanks}_pbp_seqs_res" ]
    then
      mkdir "${OUT_DIR}/${current_flanks}_pbp_seqs_res"
    else
      rm -r "${OUT_DIR}/${current_flanks}_pbp_seqs_res"
      mkdir "${OUT_DIR}/${current_flanks}_pbp_seqs_res"
    fi


    counter_2=0

    echo "This many isolates to BLAST: $blast_length"


      DB_FILE="${REF_DB}.nin"
      db_loc=${REF_DB}
      if [ ! -f $DB_FILE ]
      then
        makeblastdb -dbtype nucl -out strep_ref_dna_db -max_file_sz 2GB \
        -in ${REF_DB}
      fi
      while read coffee
      do
        namo=$(basename $coffee)
        echo "Starting BLAST on this isolate: $namo"
        blastn -db $db_loc -query $coffee -out "${OUT_DIR}/${current_flanks}_blast_results/$namo".csv \
        -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"




      done < $el_nuevo_fazza

      ## Just doing the before flanks now

      while read coffee
      do
        namo=$(basename $coffee)
        echo "Starting before BLAST on this isolate: $namo"
        blastn -db $db_loc -query $coffee -out "${OUT_DIR}/${current_flanks}_before_flank_blast_res/$namo".csv \
        -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"




      done < $before_flanks_fazza


      ## Now for the after flanks

      while read coffee
      do
        namo=$(basename $coffee)
        echo "Starting after BLAST on this isolate: $namo"
        blastn -db $db_loc -query $coffee -out "${OUT_DIR}/${current_flanks}_after_flank_blast_res/$namo".csv \
        -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"




      done < $after_flanks_fazza

      while read coffee
      do
        namo=$(basename $coffee)
        echo "Starting PBP seq BLAST on this isolate: $namo"
        blastn -db $db_loc -query $coffee -out "${OUT_DIR}/${current_flanks}_pbp_seqs_res/$namo".csv \
        -outfmt "10 qesqid sseqid qstart qend sstart send pident length evalue bitscore"

      done < $pbp_seqs

    ls -d "$PWD/${OUT_DIR}/${current_flanks}_blast_results/"*_whole_blast_seq.fasta* > "${OUT_DIR}/${current_flanks}_whole_run_through"
    ls -d "$PWD/${OUT_DIR}/${current_flanks}_before_flank_blast_res/"*_before_flank.fasta* > "${OUT_DIR}/${current_flanks}_before_run_through"
    ls -d "$PWD/${OUT_DIR}/${current_flanks}_after_flank_blast_res/"*_after_flank.fasta* > "${OUT_DIR}/${current_flanks}_after_run_through"
    ls -d "$PWD/${OUT_DIR}/${current_flanks}_pbp_seqs_res/"*_pbp_sequence.fasta* > "${OUT_DIR}/${current_flanks}_pbp_run_through"

    python "${pythondir}blast_pbp_local_results_interpreter.py" \
    --results_list "${OUT_DIR}/${current_flanks}_whole_run_through" --out_dir "${OUT_DIR}/${current_flanks}_blast_results/"
    python "${pythondir}blast_pbp_local_results_interpreter.py" \
    --results_list "${OUT_DIR}/${current_flanks}_before_run_through" --out_dir "${OUT_DIR}/${current_flanks}_before_flank_blast_res/"
    python "${pythondir}blast_pbp_local_results_interpreter.py" \
    --results_list "${OUT_DIR}/${current_flanks}_after_run_through" --out_dir "${OUT_DIR}/${current_flanks}_after_flank_blast_res/"
    python "${pythondir}blast_pbp_local_results_interpreter.py" \
    --results_list "${OUT_DIR}/${current_flanks}_pbp_run_through" --out_dir "${OUT_DIR}/${current_flanks}_pbp_seqs_res/"

    ls -d "$PWD/${OUT_DIR}/${current_flanks}_blast_results/"*_species_list.csv > "${OUT_DIR}/${current_flanks}_blast_results/${PREFIX}_species_list_list"
    ls -d "$PWD/${OUT_DIR}/${current_flanks}_before_flank_blast_res/"*_species_list.csv > "${OUT_DIR}/${current_flanks}_before_flank_blast_res/${PREFIX}_species_list_before"
    ls -d "$PWD/${OUT_DIR}/${current_flanks}_after_flank_blast_res/"*_species_list.csv > "${OUT_DIR}/${current_flanks}_after_flank_blast_res/${PREFIX}_species_list_after"
    ls -d "$PWD/${OUT_DIR}/${current_flanks}_pbp_seqs_res/"*_species_list.csv > "${OUT_DIR}/${current_flanks}_pbp_seqs_res/${PREFIX}_species_list_pbp"

    python "${pythondir}pbp_blast_to_info.py" \
    --list_file "${OUT_DIR}/${current_flanks}_blast_results/${PREFIX}_species_list_list" --hit_locs $hit_csv \
     --out_name "${OUT_DIR}/${current_flanks}_blast_results/${PREFIX}_species_compo.csv"

    python "${pythondir}pbp_blast_to_info.py" \
    --list_file "${OUT_DIR}/${current_flanks}_before_flank_blast_res/${PREFIX}_species_list_before" --hit_locs $hit_csv \
     --out_name "${OUT_DIR}/${current_flanks}_before_flank_blast_res/${PREFIX}_species_compo_before.csv"

    python "${pythondir}pbp_blast_to_info.py" \
    --list_file "${OUT_DIR}/${current_flanks}_after_flank_blast_res/${PREFIX}_species_list_after" --hit_locs $hit_csv \
    --out_name "${OUT_DIR}/${current_flanks}_after_flank_blast_res/${PREFIX}_species_compo_after.csv"

    python "${pythondir}pbp_blast_to_info.py" \
    --list_file "${OUT_DIR}/${current_flanks}_pbp_seqs_res/${PREFIX}_species_list_pbp" --hit_locs $hit_csv \
    --out_name "${OUT_DIR}/${current_flanks}_pbp_seqs_res/${PREFIX}_species_compo_pbp.csv"


    echo "Done"
  done
  echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
fi
