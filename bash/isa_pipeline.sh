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
echo "Go through Act compos (yes/no): $6"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
parentdir="$(dirname "$DIR")"
pythondir="${parentdir}/python/"
perldir="${parentdir}/perl/"
rdir="${parentdir}/R/"

read -p "Are you sure? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Nn]$ ]]
then
    echo "exiting script now"
  else
    Lines=$(cat $fasta_list | wc -l)

    counter=0
    if [ $6 == "yes" ]
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

    fi
    if [ -d ./tmp_dna_dir ]
    then
      rm -r ./tmp_dna_dir
      mkdir ./tmp_dna_dir
    else
      mkdir ./tmp_dna_dir
    fi


    cd ./tmp_dna_dir
    less "../${fasta_list}" | while read line; do \
        perl "${perldir}converting_velvet_contigs_to_dna.pl" $line;

        done
    cat *.dna > output.mfa





    makeblastdb -dbtype nucl -out temp_blast_db -max_file_sz 2GB \
    -in output.mfa

    rm output.mfa

    blastn -db temp_blast_db -query "../$2" -outfmt 10 -out tmp_blast_csv.csv

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
      ./grep_output_temp $line

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
     ./tmp_dna_dir/tmp_blast_csv.csv ./contig_bounds "$4/$5_merged_blast_file"

  Fazza=$fasta_list
  Lines=$(cat $fasta_list | wc -l)
  counter=0

  ## Now we run the library creation step
  python "${pythondir}library_creator.py" \
  --hit_csv "$4/$5_merged_blast_file" --reference_csv $ref_isolate_gff --align_cutoff $3 \
  --act_loc ./act_compos/referoo.fasta. --contig_loc ./contig_bounds/ --output "$4/$5_library.csv"

fi