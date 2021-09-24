#!/bin/sh -e
################################################################################
# Brief trial script for running flye.
################################################################################

data=/kakapo-data/ONT/fastq/
out=/kakapo-data/ONT/flye_out/

for fasta in ${data}*.fastq.gz
do
base=$(basename ${fasta} .fastq.gz)
#mkdir ${out}${base}
printf "\nRunning flye for $base...\n"
flye --nano-raw ${fasta} --out-dir ${out}${base} --genome-size 1.2g --threads 16 --debug
done
