#!/bin/bash -e
################################################################################
# Running the Delly SV discovery and genotyping pipeline. Delly is a programme
# for structural variant discovery with paired-end sequence data. Delly requires
# sorted bam files with marked duplicates (as in the script 2_align_stat.sh)
# as input. 
################################################################################
fref=/kakapo-data/References/kakapo_full_ref.fa
fbam=/kakapo-data/bwa_female/markdup/
fcall=/kakapo-data/delly/bwa_aligned/female_SVcalls/
mref=/kakapo-data/References/kakapo_no_Wchromosome.fa
mbam=/kakapo-data/bwa_male/markdup/
mcall=/kakapo-data/delly/bwa_aligned/male_SVcalls/
genos=/kakapo-data/delly/bwa_aligned/genotypes/

printf "\nBEGIN DELLY CALL FOR FEMALES\n"
for i in {01..11}
    do
    for female in ${fbam}batch${i}/*.sorted.bam
    do
        fbase=$(basename ${female} .sorted.bam)
        printf "\nRunning Delly for ${fbase} in batch ${i}...\n"
        delly call -g ${fref} -o ${fcall}${fbase} ${female}
    done &
done &
printf "\nBEGIN DELLY CALL FOR MALES\n"
for j in {01..14}
    do
    for male in ${mbam}batch${j}/*.sorted.bam
    do
        mbase=$(basename ${male} .sorted.bam)
        printf "\nRunning Delly for ${mbase} in batch ${j}...\n"
        delly call -g ${mref} -o ${mcall}${mbase} ${male}
    done &
done
wait
printf "\nREMOVING THE W CHROMOSOME FROM FEMALE SV CALLS...\n"
for female in ${fcall}*bcf
    do
        fbase=$(basename ${female} .bcf)
        printf "\nRemoving the W scaffold for ${fbase}..."
        bcftools view -t ^NC_044301.2 -O b -o ${fcall}${fbase}_noW.bcf ${female} &
done
wait
printf "\nMERGING ALL SV CALLS...\n"
delly merge -o /kakapo-data/delly/bwa_aligned/sites.bcf ${fcall}*_noW.bcf ${mcall}.bcf
printf "\nGENOTYPING INDIVIDUALS...\n"
for i in {01..11}
    do
    for female in ${fbam}batch${i}/*.sorted.bam
    do
        fbase=$(basename ${female} .sorted.bam)
        printf "\nRunning Delly genotyping for ${fbase} in batch ${i}...\n"
        delly call -g ${fref} \
            -v /kakapo-data/delly/bwa_aligned/sites.bcf \
            -o ${genos}${fbase}.geno.bcf ${female}
    done &
done &
printf "\nBEGIN DELLY CALL FOR MALES\n"
for j in {01..14}
    do
    for male in ${mbam}batch${j}/*.sorted.bam
    do
        mbase=$(basename ${male} .sorted.bam)
        printf "\nRunning Delly genotyping for ${mbase} in batch ${j}...\n"
        delly call -g ${mref} \
            -v /kakapo-data/delly/bwa_aligned/sites.bcf \
            -o ${genos}${mbase}.geno.bcf ${male}
    done &
done
wait
printf "\nCREATING RAW BCF...\n"
bcftools merge -m id -O b -o /kakapo-data/delly/kakapo_SVs.bcf ${genos}*
delly filter -f germline -o /kakapo-data/delly/kakapo_SVs_filter.bcf /kakapo-data/delly/kakapo_SVs.bcf
        