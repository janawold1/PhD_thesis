#!/bin/sh -e
#######################################################################################################
# Scripts provided by Roger Moraga, pipeline augmented by Jana Wold.
#######################################################################################################
ref=/data/tara_iti_shortreads/reference/bSteHir1.pri.cur.20190820.fasta #Reference genome for alignment
scriptdir=/data/tara_iti_shortreads/Moraga_bcftools_pipeline/ # Directory holding bcftools pipeline
bamdir=/data/tara_iti_shortreads/alignments/merged_bam/ # BAM file directory
bcf_dir=/data/tara_iti_shortreads/bcftools_variants/ # BCF output

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
        --threads 32 \
        -f $ref\
        -a AD,ADF,ADR,DP,SP \
        --threads 64 \
        -o ${bcf_dir}fairy_tern_${i}_raw.bcf \
        ${bcf_dir}chunks/${i}/* &
done
wait
printf "\nMPILEUP is complete. Beginning variant calling...\n"
for file in ${bcf_dir}*.bcf
do
    base=$(basename $file .bcf)
    bcftools call $file -f GQ -mv -O v -o ${bcf_dir}${base}_VariantCalls.vcf &
done
wait

#printf "\nVariant calling is complete.\nPreparing files for filtering...\n"
#for file in ${bcf_dir}*.vcf
#    do
#    base=$(basename $file .vcf)
#    bcftools reheader \
#        -s ${bamdir}tara_iti_bam_list.txt \
#        -O b -o ${bcf_dir}${base}_reheader.bcf \
#        ${file}
#    bcftools index ${bcf_dir}${base}_reheader.bcf
#    ls ${bcf_dir}${base}_reheader.bcf >> ${bcf_dir}list_of_bcf.txt
#done

printf "Completed file preparation with bgzip and indexing, concatenating files now...\n"

bcftools concat --file-list ${bcf_dir}list_of_bcf.txt \
    -O v -o ${bcf_dir}Fairy_tern_VariantCalls.vcf --threads 16
printf "\nVCF is ready for filtering!\n"
