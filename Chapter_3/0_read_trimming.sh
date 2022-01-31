#!/bin/sh
##############################################################################
# Scripts for processing tara iti and Australian fairy tern nova-seq samples
# recieved 12th April 2021.
## NOTE: All fairy tern samples have a double peak in their fastqc GC curves.
##############################################################################

# First began with 2 colour chemistry specific trimmming.
raw=/data/tara_iti_shortreads/raw_reads/2021_august/
trim_out=/data/tara_iti_shortreads/trimmed_reads2/lib2/
reports=/data/tara_iti_shortreads/trim_reports2/lib2/

for samp in ${raw}*R1.fastq.gz
do
base=$(basename ${samp} _R1.fastq.gz)
echo "Running Trim_galore for ${base}..."
trim_galore --paired \
    --nextseq 28 \
    --2colour 20 \
    --cores 8 \
    --fastqc_args "--nogroup --outdir trim_reports --threads 32" \
    --length 50 \
    --output_dir ${trim_out}/ \
    --clip_R1 20 \
    --clip_R2 20 \
    --three_prime_clip_R1 5 \
    --three_prime_clip_R2 5 \
    --retain_unpaired \
    --length_1 55 \
    --length_2 55 \
    ${raw}${base}_R1.fastq.gz \
    ${raw}${base}_R2.fastq.gz
mv ${trim_out}*.txt ${reports}
mv ${trim_out}*fastqc* ${reports}
done