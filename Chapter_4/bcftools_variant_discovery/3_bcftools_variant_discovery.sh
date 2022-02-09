#!/bin/sh -e
#######################################################################################################
# Scripts provided by Roger Moraga, pipeline augmented by Jana Wold.
#######################################################################################################
ref=/data/reference/bSteHir1.pri.cur.20190820.fasta #Reference genome for alignment
scriptdir=/data/Moraga_bcftools_pipeline/ # Directory holding bcftools pipeline
bamdir=/data/common_tern/alignments/merged_bam/ # BAM file directory
bcf_dir=/data/common_tern/bcftools_HQvariant_calls/ # BCF output

printf "\nChunking bam files for mpileup...\n"
ls ${bamdir}*.bam > ${bcf_dir}bam_list.txt
mkdir ${bcf_dir}chunks
perl ${scriptdir}split_bamfiles_tasks.pl \
    -b ${bcf_dir}bam_list.txt \
    -g ${ref} -n 24 -o ${bcf_dir}chunks/ | parallel -j 24 {}

printf "\nRunning mpileup on chunks of bam files...\n"
for (( i=1; i<=24; i++ ))
    do
    bcftools mpileup \
        --threads 24 \
        -f $ref\
        -a AD,ADF,ADR,DP,SP \
        -o ${bcf_dir}fairy_tern_${i}_raw.bcf \
        ${bcf_dir}chunks/${i}/* &
done
wait
printf "\nMPILEUP is complete. Beginning variant calling...\n"
for file in ${bcf_dir}*.bcf
do
    base=$(basename $file .bcf)
    bcftools call $file -a INFO/PV4,FORMAT/GQ,FORMAT/GP -mv -O v -o ${bcf_dir}${base}_VariantCalls.vcf &
done
wait

ls ${bcf_dir}*_VariantCalls.vcf > ${bcf_dir}list_of_VCFs.txt

printf "Completed file preparation with bgzip and indexing, concatenating files now...\n"

bcftools concat --threads 24 --file-list ${bcf_dir}list_of_VCFs.txt \
    -O v -o ${bcf_dir}Fairy_tern_VariantCalls.vcf --threads 24
printf "\nVCF is ready for filtering!\n"
