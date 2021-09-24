#!/bin/bash/ -e
#####################################################################################

#####################################################################################

##Setting fixed variables
fref=/kakapo-data/References/kakapo_full_ref.fa
mref=/kakapo-data/References/kakapo_no_Wchromosome.fa
repout=/kakapo-data/delly/replicate

for j in {2..3}
do
    for i in {01..11}
    do
        for female in /kakapo-data/alignments_female/batch${i}/*.sorted.bam
        do
            id=$(basename ${female} .sorted.bam)
            echo "Running Delly call for ${id} in batch ${i} for replicate ${j}..."
            delly call -g ${fref} -o ${repout}${j}/female_SVcalls/${id}_SVcalls${j}.bcf ${female}   
        done &
    done &
    for k in {01..14}
    do
        for male in /kakapo-data/alignments_male/batch${k}/*.sorted.bam
        do
            id=$(basename ${male} .sorted.bam)
            echo "Running Delly call for ${id} in batch ${k} for replicate ${j}..."
            delly call -g ${mref} -o ${repout}${j}/male_SVcalls/${id}_SVcalls${j}.bcf ${male}
        done &
     done &
done
wait