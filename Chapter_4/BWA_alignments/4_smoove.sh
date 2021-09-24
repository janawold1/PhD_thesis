#!/bin/sh
###################################################################################
# Walkthrough of how I ran the SMOOVE software package. 
###################################################################################
fref=/kakapo-data/References/kakapo_full_ref.fa
mref=/kakapo-data/References/kakapo_no_Wchromosome.fa
annotate=/kakapo-data/GCF_004027225.2_bStrHab1.2.pri_genomic.gff.gz
fbam=/kakapo-data/bwa_female/
fout=/kakapo-data/smoove/bwa_alignments/female_SVcalls/
fvcf=/kakapo-data/smoove_results/interm_VCFs/joint_merged.sites.vcf.gz
mbam=/kakpo-data/bwa_male/
mout=/kakapo-data/smoove/bwa_alignments/male_SVcalls/
mvcf=/kakapo-data/smoove_results/interm_VCFs/joint_noW.vcf.gz


printf "\nRUNNING SMOOVE FOR FEMALES...\n"
for i in {01..11}
    do
    for female in ${fbam}batch${i}/*.sorted.bam
    do
        fsmoove=$(basename ${female} .sorted.bam)
        echo "Running SMOOVE for ${fsmoove} in batch ${i}..."
        smoove call \
            --outdir ${fout} \
            --name ${fsmoove} \
            --fasta ${fref} \
            -p 16 --genotype ${female}
    done &
done &
printf "\nRUNNING SMOOVE FOR MALES...\n"
for j in {01..14}
    do
    for male in ${mbam}batch${j}/*.sorted.bam
        do
        msmoove=$(basename ${male} .sorted.bam)
        echo "Running SMOOVE for ${msmoove} in batch ${j}..."
        smoove call \
            --outdir smoove_results/female/ \
            --name ${msmoove} \
            --fasta ${mref} \
            -p 16 --genotype ${male}
    done &
done
wait
####################################################
# Need to continue from here. Script to refer to in VM is /kakapo-data/smoove/bwa_aligned/10X/scripts/run_smoove10X.sh
printf "\nMERGING ALL CALLED VARIANTS..." 
smoove merge --name joint_merged \
    -f $reference_fasta \
    --outdir /kakapo-data/smoove_results/interm_VCFs/ \
    smoove_results/{female,male}/*.genotyped.vcf.gz 

bcftools view -O z -o interm_VCFs/joint_noW.vcf.gz -t ^NC_044301.2 joint_merged.sites.vcf.gz

for female in alignments_female/*.sorted.bam
    do
    fsmoove=$(basename ${female} .sorted.bam)
    echo "Creating individual genotyped VCF for ${fsmoove}...."
    smoove genotype -d -x -p 1 \
        --name ${fsmoove}-joint \
        --outdir /kakapo-data/smoove_results/female/ \
        --fasta $fref \
        --duphold \
        --vcf ${fvcf} \
        ${female}
done

for female in /kakapo-data/smoove_results/female/*-joint.vcf.gz
    do
    base=$(basename ${female} -joint.vcf.gz)
    echo "Removing W chromosome from ${base} VCF..."
    bcftools view \
        -O z -o joint/${base}.vcf.gz \
        -t ^NC_044301.2 \
        ${female}
done

for male in alignments_male/*.sorted.bam
    do
    msmoove=$(basename ${male} .sorted.bam)
    echo "Creating individual genotyped VCF for ${msmoove}...."
    smoove genotype -d -x -p 1 \
        --name ${msmoove}-joint \
        --outdir /kakapo-data/smoove_results/joint/ \
        --fasta $mref \
        --duphold
        --vcf ${mvcf} \
        ${male}
done

echo "Creating total raw VCF..."
smoove paste --name kakapo_noW \
    /kakapo-data/smoove_results/joint/*.vcf.gz

smoove annotate --gff ${annotate} \
    /kakapo-data/smoove_results/kakapo_noW.smoove.square.vcf.gz \
    | bgzip -c > kakapo_noW.smoove.square.anno.vcf.gz