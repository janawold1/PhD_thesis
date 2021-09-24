#########################################################################################
# Walkthrough of how I ran the BayesTyper software package. Turns out this has a lot of 
# dependencies and steps. 
#########################################################################################

# First installed with bioconda
conda activate bayestyper
#conda install bayestyper

## Before beginning, there's a number of steps to prepare the data. This includes 
## 1) Formatting the candidate SV vcf 
## 2) Running KMC and bayesTyperTools makeBloom for all individuals
## 3) Preparing the Reference genome and 'decoy' scaffolds
## For the most part, these steps can be done independently of one another.

# With a manta vcf, the first step is to identify SV candidates by de novo assembly. 
# This converts all alleles identified the manta vcf to a full sequence. I found 
# that this step works is easier to check if you remove all of the individual genotype
# information prior.
workvcf=/kakapo-data/bayestyper/vcfs

bcftools view -O z -o ${workvcf}/joint_female_diploidSV.vcf.gz \
    -G /kakapo-data/manta/joint_calling/female/results/variants/candidateSV.vcf.gz
bcftools view -O z -o ${workvcf}/male_candidates.vcf.gz \
    -G /manta/joint_calling/male/results/variants/candidateSV.vcf.gz
    
conda activate manta
MANTA=/home/rccuser/anaconda3/envs/manta/share/manta-1.6.0-0/libexec
${MANTA}/convertInversion.py ${MANTA}/samtools \
    /kakapo-data/References/kakapo_full_ref.fa \
    ${workvcf}/joint_female_diploidSV.vcf.gz > ${workvcf}/joint_female_diploidSV_inv.vcf
bcftools norm -c wx -f /kakapo-data/References/kakapo_full_ref.fa \
    -m -any -O v -o ${workvcf}/joint_female_diploidSV_inv_norm.vcf \
    ${workvcf}/joint_female_diploidSV_inv.vcf
    
${MANTA}/convertInversion.py ${MANTA}/samtools \
    /kakapo-data/References/kakapo_full_ref.fa \
    ${workvcf}/joint_male_diploidSV.vcf.gz > ${workvcf}/joint_male_diploidSV_inv.vcf
bcftools norm -c wx -f /kakapo-data/References/kakapo_full_ref.fa \
    -m -any -O v -o ${workvcf}/joint_male_diploidSV_inv_norm.vcf \
    ${workvcf}/joint_male_diploidSV_inv.vcf

##################################################################################################
# This creates 2 VERY large vcfs. The -z flag then gzips these files for the next step.
# Note, bayesTyper does not support bgzip files.These files must be normalised again prior to 
# merging.
##################################################################################################

bayesTyperTools convertAllele --variant-file ${workvcf}/joint_female_diploidSV_inv_norm.vcf \
    --genome-file /kakapo-data/References/kakapo_full_ref.fa \
    --output-prefix joint_female_converted -z

bayesTyperTools convertAllele --variant-file ${workvcf}/joint_male_diploidSV_inv_norm.vcf \
    --genome-file /kakapo-data/References/kakapo_full_ref.fa \
    --output-prefix joint_male_converted -z

bcftools norm -f /kakapo-data/References/kakapo_full_ref.fa \
    -m -any \
    -O v -o ${workvcf}/joint_female_converted_norm.vcf \
    --threads 36 \
    ${workvcf}/joint_female_converted.vcf.gz

bcftools norm -f /kakapo-data/References/kakapo_full_ref.fa \
    -m -any \
    -O v -o ${workvcf}/joint_male_converted_norm.vcf \
    --threads 36 \
    ${workvcf}/joint_male_converted.vcf.gz 

# However, this process only converts the VCFs called from the sex specific variant calling. 
# To get a joint call set, we first merge the raw VCFs from MANTA with bayestyper as per below.

bayesTyperTools combine \
    --variant-files MANTA:${workvcf}/joint_female_converted_norm.vcf, MANTA:${workvcf}/joint_male_converted_norm.vcf \
    --output-prefix ${workvcf}/joint_candidateSVs \
    --filter-ambiguous-all

##################################################################################################
# To create the candidate SV call set from the batched Manta calls, this was repeated as per 
# below. These files were already converted with Manta convertInversion.
##################################################################################################

conda activate bayestyper

for fvcf in /kakapo-data/manta/batch_calling/batches_female/*_manta_INV_conversion_female.vcf.gz
    do
    base=$(basename ${fvcf})
    batch=$(basename ${fvcf} _manta_INV_conversion_female.vcf.gz)
    echo "Copying ${base}..."
    bcftools view -O z -o ${workvcf}/batch/${base} \
        -G ${fvcf}
    echo "Normalising female ${batch}..."
    bcftools norm -c wx -f /kakapo-data/References/kakapo_full_ref.fa \
        -m -any -O v -o ${workvcf}/batch/${batch}_female_norm.vcf \
        ${workvcf}/batch/${base}
    echo "Converting female ${batch}..."
    bayesTyperTools convertAllele \
        --variant-file ${workvcf}/batch/${batch}_female_norm.vcf \
        --genome-file /kakapo-data/References/kakapo_full_ref.fa \
        --output-prefix ${workvcf}/batch/${batch}_female_converted -z
    echo "Normalising converted female ${batch}..."
    bcftools norm -f /kakapo-data/References/kakapo_full_ref.fa \
        -m -any \
        -O v -o ${workvcf}/batch/${batch}_female_converted_norm.vcf \
        --threads 36 \
        ${workvcf}/batch/${batch}_female_converted.vcf.gz 
done
        
for mvcf in /kakapo-data/manta/batch_calling/batches_male/*_manta_INV_conversion_male.vcf.gz
    do
    base=$(basename ${mvcf})
    batch=$(basename ${mvcf} _manta_INV_conversion_male.vcf.gz)
    echo "Copying ${base}..."
    bcftools view -O z -o ${workvcf}/batch/${base} \
        -G ${mvcf}
    echo "Normalising male ${batch}..."
    bcftools norm -c wx -f /kakapo-data/References/kakapo_full_ref.fa \
        -m -any -O v -o ${workvcf}/batch/${batch}_male_norm.vcf \
        ${workvcf}/batch/${base}
    echo "Converting male ${batch}..."
    bayesTyperTools convertAllele \
        --variant-file ${workvcf}/batch/${batch}_male_norm.vcf \
        --genome-file /kakapo-data/References/kakapo_full_ref.fa \
        --output-prefix ${workvcf}/batch/${batch}_male_converted -z
    echo "Normalising converted male ${batch}..."
    bcftools norm -f /kakapo-data/References/kakapo_full_ref.fa \
        -m -any \
        -O v -o ${workvcf}/batch/${batch}_male_converted_norm.vcf \
        --threads 36 \
        ${workvcf}/batch/${batch}_male_converted.vcf.gz 
done

#Then merged all batches with...
bayesTyperTools combine \
    --variant-files MANTA:batch01_INV_female_converted_norm.vcf,MANTA:batch01_INV_male_converted_norm.vcf,MANTA:batch02_INV_female_converted_norm.vcf,MANTA:batch02_INV_male_converted_norm.vcf,MANTA:batch03_INV_female_converted_norm.vcf,MANTA:batch03_INV_male_converted_norm.vcf,MANTA:batch04_INV_female_converted_norm.vcf,MANTA:batch04_INV_male_converted_norm.vcf,MANTA:batch05_INV_female_converted_norm.vcf,MANTA:batch05_INV_male_converted_norm.vcf,MANTA:batch06_INV_female_converted_norm.vcf,MANTA:batch06_INV_male_converted_norm.vcf,MANTA:batch07_INV_female_converted_norm.vcf,MANTA:batch07_INV_male_converted_norm.vcf,MANTA:batch08_INV_male_converted_norm.vcf,MANTA:batch09_INV_male_converted_norm.vcf \
    --output-prefix /kakapo-data/bayestyper/vcfs/batch/batch_candidateSVs \
    --filter-ambiguous-all

##################################################################################################
# Runnning KMC for all bam files can be done concurrently with formatting the VCF as per below...
##################################################################################################

for male in alignments_male/*.sorted.bam
    do
    indiv=$(basename ${male} .sorted.bam)
    echo "Running KMC3 for ${indiv}..."
    kmc -k55 -ci1 -m72 -t288 -fbam ${male} kmc_results/${indiv}_KMC.res kmc/
done

for female in alignments_male/*.sorted.bam
    do
    indiv=$(basename ${female} .sorted.bam)
    echo "Running KMC3 for ${indiv}..."
    kmc -k55 -ci1 -m72 -t288 -fbam ${female} kmc_results/${indiv}_KMC.res kmc/
done

for kmc in kmc_results/*.res.kmc_pre
    do
    bloom=$(basename ${kmc} .kmc_pre)
    echo "Running makeBloom for ${bloom}..."
    bayesTyperTools makeBloom -k kmc_results/${bloom} -p 8
done

# A ploidy file for genotyping must provided in the format below:
# chr    female    male
# 1      2         2
# Z      1         2
# W      1         0
# Removed all unplaced scaffolds where ploidy was unknown in the candidate VCF and 
# reference. Found a script (faSomeRecords) to remove the scaffolds where ploidy is
# unknown at: https://www.biostars.org/p/360658/

##################################################################################################
# Then variant clusters were estimated and genotyped for the joint call set.
# Note, no more than 30 indivs can be genotyped at once...
##################################################################################################

for samp in /kakapo-data/bayestyper/sample_lists/*.tsv
do
    batch=$(basename ${samp} .tsv)
    bayesTyper cluster \
        --variant-file ${workvcf}/joint_candidateSVs.vcf \
        --samples-file ${samp} \
        --genome-file modified_references/kakapo_chr_scaffolds.fa \
        --decoy-file modified_references/kakapo_decoy.fa \
        --output-prefix bayestyper_cluster_data/${batch} \
        --threads 16
done

# Genotyping was then conducted as per
for samp in /kakapo-data/bayestyper/sample_lists/*.tsv
do
    batch=$(basename ${samp} .tsv)
    dir=/kakapo-data/bayestyper
    bayesTyper genotype \
        --variant-clusters-file ${dir}/bayestyper_cluster_data/${batch}_unit_1/variant_clusters.bin \
        --cluster-data-dir ${dir}/bayestyper_cluster_data/${batch}_cluster_data/ \
        --samples-file ${samp}
        --chromosome-ploidy-file ${dir}/intermediate_files/ploidy.tsv \
        --genome-file ${dir}/modified_references/kakapo_chr_scaffolds.fa \
        --decoy-file ${dir}/modified_references/kakapo_decoy.fa \
        --output-prefix ${dir}/genotype_out/${batch} \
        --threads 32
done

##################################################################################################
# Then variant clusters were estimated and genotyped for the batch call set.
# Note, no more than 30 indivs can be genotyped at once...
##################################################################################################

dir=/kakapo-data/bayestyper
for samp in /kakapo-data/bayestyper/sample_lists/*.tsv
do
    batch=$(basename ${samp} .tsv)
    bayesTyper cluster \
        --variant-file ${workvcf}/batch/batch_candidateSVs.vcf \
        --samples-file ${samp} \
        --genome-file ${dir}/modified_references/kakapo_chr_scaffolds.fa \
        --decoy-file ${dir}/modified_references/kakapo_decoy.fa \
        --output-prefix ${dir}/bayestyper_cluster_data/batch/${batch} \
        --threads 16
done

for samp in /kakapo-data/bayestyper/sample_lists/*.tsv
do
    batch=$(basename ${samp} .tsv)
    bayesTyper genotype \
        --variant-clusters-file ${dir}/bayestyper_cluster_data/batch/${batch}_unit_1/variant_clusters.bin \
        --cluster-data-dir ${dir}/bayestyper_cluster_data/batch/${batch}_cluster_data/ \
        --samples-file ${samp}
        --chromosome-ploidy-file ${dir}/intermediate_files/ploidy.tsv \
        --genome-file ${dir}/modified_references/kakapo_chr_scaffolds.fa \
        --decoy-file ${dir}/modified_references/kakapo_decoy.fa \
        --output-prefix ${dir}/genotype_out/${batch} \
        --threads 32
done