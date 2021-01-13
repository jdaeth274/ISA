#!/bin/bash

set -e

if [ $# != 4 ]
then
  echo "This need 4 args you have $#."
  echo "Requires: <outdir> <prefix> <flank length file> <hit_csv>"
  exit
fi

echo "This is the outdir: $1" # $4
echo "This is the prefix: $2" # $5
echo "Flanking region length file: $3" # $6
echo "This is the hit_csv: $4" # $hit_csv

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
parentdir="$(dirname "$DIR")"
pythondir="${parentdir}/python/"
perldir="${parentdir}/perl/"
rdir="${parentdir}/R/"

if [ ! -f $4 ] || [ ! -f $3 ] || [ ! -d $1 ]
then
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  echo "Missing files check the usage above"
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  exit
else
  if [ ! -d $1 ]
  then
  mkdir $1
  fi



  less $3 | while read line; do \
    current_flanks=$line

    echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
    echo "On flank length: ${current_flanks}"

  if [ ! -d "$1/$2_${current_flanks}_flanks" ]
      then
      mkdir "$1/$2_${current_flanks}_flanks"
    else
      rm -r "$1/$2_${current_flanks}_flanks"
      mkdir "$1/$2_${current_flanks}_flanks"
    fi

    ls -d "$PWD/$1/${current_flanks}_blast_results/"*_whole_blast_seq.fasta* > "$1/${current_flanks}_whole_run_through"
    ls -d "$PWD/$1/${current_flanks}_before_flank_blast_res/"*_before_flank.fasta* > "$1/${current_flanks}_before_run_through"
    ls -d "$PWD/$1/${current_flanks}_after_flank_blast_res/"*_after_flank.fasta* > "$1/${current_flanks}_after_run_through"

    python "${pythondir}blast_local_results_interpreter.py" \
    --results_list "$1/${current_flanks}_whole_run_through" --out_dir "$1/${current_flanks}_blast_results/"
    python "${pythondir}blast_local_results_interpreter.py" \
    --results_list "$1/${current_flanks}_before_run_through" --out_dir "$1/${current_flanks}_before_flank_blast_res/"
    python "${pythondir}blast_local_results_interpreter.py" \
    --results_list "$1/${current_flanks}_after_run_through" --out_dir "$1/${current_flanks}_after_flank_blast_res/"

    ls -d "$PWD/$1/${current_flanks}_blast_results/"*_species_list.csv > "$1/${current_flanks}_blast_results/$2_species_list_list"
    ls -d "$PWD/$1/${current_flanks}_before_flank_blast_res/"*_species_list.csv > "$1/${current_flanks}_before_flank_blast_res/$2_species_list_before"
    ls -d "$PWD/$1/${current_flanks}_after_flank_blast_res/"*_species_list.csv > "$1/${current_flanks}_after_flank_blast_res/$2_species_list_after"

    python "${pythondir}blast_to_info.py" \
    --list_file "$1/${current_flanks}_blast_results/$2_species_list_list" --hit_locs $4\
     --out_name "$1/${current_flanks}_blast_results/$2_species_compo.csv"

    python "${pythondir}blast_to_info.py" \
    --list_file "$1/${current_flanks}_before_flank_blast_res/$2_species_list_before" --hit_locs $4 \
     --out_name "$1/${current_flanks}_before_flank_blast_res/$2_species_compo_before.csv"

    python "${pythondir}blast_to_info.py" \
    --list_file "$1/${current_flanks}_after_flank_blast_res/$2_species_list_after" --hit_locs $4 \
    --out_name "$1/${current_flanks}_after_flank_blast_res/$2_species_compo_after.csv"



    echo "Done"
  done
  echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
fi
