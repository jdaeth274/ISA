#!/bin/bash
set -e
echo "This is the flanks length file: $1"
echo "This is the outfolder: $2"
echo "This is the prefix: $3"
hit_csv="$2/$3_hits_df.csv"
echo "Thi is the hit_csv"


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
parentdir="$(dirname "$DIR")"
pythondir="${parentdir}/python/"
perldir="${parentdir}/perl/"
rdir="${parentdir}/R/"


if [ ! -f $1 ] || [ ! -d $2 ] || [ ! -f $hit_csv ]
then
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  echo "Missing files check the usage above"
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  exit
else


  less $1 | while read line; do \
    current_flanks=$line

    echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
    echo "On flank length: ${current_flanks}"



    blast_length=$(ls -l "$2/$3_${current_flanks}_flanks/"*_whole_blast_seq.fasta | wc -l)

    el_nuevo_fazza="$2/${current_flanks}_isolate_blast_list"
    before_flanks_fazza="$2/$3_${current_flanks}_before_flanks_list"
    after_flanks_fazza="$2/$3_${current_flanks}_after_flanks_list"

    ls -d "$PWD/$2/${current_flanks}_blast_results/"*_whole_blast_seq.fasta.csv > "$2/${current_flanks}_whole_run_through"
    ls -d "$PWD/$2/${current_flanks}_before_flank_blast_res/"*_before_flank.fasta.csv > "$2/${current_flanks}_before_run_through"
    ls -d "$PWD/$2/${current_flanks}_after_flank_blast_res/"*_after_flank.fasta.csv > "$2/${current_flanks}_after_run_through"

    rm -f "$2/${current_flanks}_blast_results"/*species_list.csv
    rm -f "$2/${current_flanks}_before_flank_blast_res"/*species_list.csv
    rm -f "$2/${current_flanks}_after_flank_blast_res"/*species_list.csv

    python "${pythondir}blast_local_results_interpreter.py" \
    --results_list "$2/${current_flanks}_whole_run_through" --out_dir "$2/${current_flanks}_blast_results/"
    python "${pythondir}blast_local_results_interpreter.py" \
    --results_list "$2/${current_flanks}_before_run_through" --out_dir "$2/${current_flanks}_before_flank_blast_res/"
    python "${pythondir}blast_local_results_interpreter.py" \
    --results_list "$2/${current_flanks}_after_run_through" --out_dir "$2/${current_flanks}_after_flank_blast_res/"



    ls -d "$PWD/$2/${current_flanks}_blast_results/"*_species_list.csv > "$2/${current_flanks}_blast_results/$3_species_list_list"
    ls -d "$PWD/$2/${current_flanks}_before_flank_blast_res/"*_species_list.csv > "$2/${current_flanks}_before_flank_blast_res/$3_species_list_before"
    ls -d "$PWD/$2/${current_flanks}_after_flank_blast_res/"*_species_list.csv > "$2/${current_flanks}_after_flank_blast_res/$3_species_list_after"

    python "${pythondir}blast_to_info.py" \
    --list_file "$2/${current_flanks}_blast_results/$3_species_list_list" --hit_locs $hit_csv \
     --out_name "$2/${current_flanks}_blast_results/$3_species_compo.csv"

    python "${pythondir}blast_to_info.py" \
    --list_file "$2/${current_flanks}_before_flank_blast_res/$3_species_list_before" --hit_locs $hit_csv \
     --out_name "$2/${current_flanks}_before_flank_blast_res/$3_species_compo_before.csv"

    python "${pythondir}blast_to_info.py" \
    --list_file "$2/${current_flanks}_after_flank_blast_res/$3_species_list_after" --hit_locs $hit_csv \
    --out_name "$2/${current_flanks}_after_flank_blast_res/$3_species_compo_after.csv"



    echo "Done"
  done
  echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
fi
