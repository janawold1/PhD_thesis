###################################################################################
# Walkthrough of how I ran the BayesTyper software package. Turns out this has a 
# lot of dependencies and steps. 
###################################################################################

# First installed with bioconda
#conda create bayestyper
#conda activate bayestyper
#conda install bayestyper

# With a manta vcf, the first step is to identify SV candidates by de novo assembly. 
# This converts all alleles identified the manta vcf to a full sequence. This 
# creates 2 VERY large vcfs. Must the -z flag gzip's these files for the next step.
# Note, bayesTyper does not support bgzip files. The next step is to normalise and 
# left align the candidate variants with bcftools.

fref=References/kakapo_full_ref.fa
mref=References/kakapo_no_Wchromosome.fa

bcftools merge -m all -O z -o joint_total.vcf.gz \
    female_diploidSV.vcf.gz \
    male_diploidSV.vcf.gz

bcftools norm --check-ref we \
    -f ${fref} \
    -m -any \
    -O z -o joint_total_norm.vcf.gz \
    joint_total.vcf.gz

bayesTyperTools convertAllele \
    --variant-file joint_total_norm.vcf.gz \
    --genome-file ../References/kakapo_full_ref.fa \
    --output-prefix joint_total_converted -z

bcftools norm -f ${fref} \
    -f ${fref}
    -m -any 
    -O v -o bayestyper_variants.vcf
    joint_total_converted.vcf.gz

for fvariants in female_inv_diploidSV.vcf.gz 
    do
    echo "Running bayesTyper convert allele for variants found in females..."
    bayesTyperTools convertAllele \
        --genome-file ${fref} \
        --output-prefix female_converted_SVs \
        --variant-file ${fvariants} -z
    echo "Normalising the converted VCF for females..."
    bcftools norm -f ${fref} \
        -m -any -O z --threads 32 \
        -o femaleSV_norm.vcf.gz \
        female_converted_SVs.vcf.gz
done

for mvariants in male_inv_diploidSV.vcf.gz
    do
    echo "Running bayesTyper convert allele for variants found in males..."
    bayesTyperTools convertAllele \
        --genome-file ${mref} \
        --output-prefix male_converted_SVs \
        --variant-file male_inv_diploidSV.vcf.gz -z
    echo "Normalising the converted VCF for males..."
    bcftools norm -f ${mref} \
        -m -any -O z --threads 32 \
        -o maleSV_norm.vcf.gz \
        male_converted_SVs.vcf.gz
done

# The following steps with KMC and bayesTyperTools makeBloom are necessary to run
# bayesTyper cluster and bayesTyper genotype. It is independent of preparing the
# variant VCFs and can be done concurrently. 

for mkmc in alignments_male/*.sorted.bam
    do
    indiv=$(basename ${mkmc} .sorted.bam)
    echo "Running KMC for ${indiv}..."
    kmc -k55 -ci1 -m24 -t48 -fbam ${mkmc} \
        kmc_results/male/${indiv}_KMC.res \
        kmc_temp/
    echo "Running makeBloom for ${indiv}..."
    bayesTyperTools makeBloom \
        -k kmc_results/male/${indiv}_KMC.res -p 8
done

for fkmc in alignments_female/*.sorted.bam
    do
    indiv=$(basename ${fkmc} .sorted.bam)
    echo "Running KMC for ${indiv}..."
    kmc -k55 -ci1 -m24 -t48 -fbam ${fkmc} \
        kmc_results/female/${indiv}_KMC.res \
        kmc_temp/
    echo "Running makeBloom for ${indiv}..."
    bayesTyperTools makeBloom \
        -k kmc_results/female/${indiv}_KMC.res -p 8
done

# For bayesTyper cluster, a ploidy file in the format below is required
# chr    female    male
# 1      2         2
# Z      1         2
# W      1         0

# Removed all unplaced scaffolds where ploidy was unknown in the candidate VCF and 
# reference. Found a script (faSomeRecords) to remove the scaffolds where ploidy is 
# unknown at: https://www.biostars.org/p/360658/
# Then variant clusters were estimated. Note, no more than 30 indivs can be
# genotyped at once

# The <samples>.tsv file should contain one sample per row with columns
# <sample_id>, <sex> and <kmc_output_prefix> and no header 

mVCF=maleSV_norm.vcf.gz
fVCF=femaleSV_norm.vcf.gz
modMref=modified_references/male_kakapo_ref.fa
modFref=modified_references/female_kakapo_ref.fa
decoy=modified_references/kakapo_decoy.fa

bayesTyper cluster -v ${mVCF} \
    -s samples.tsv \
    -g ${modMref} \
    -d ${decoy} \
    -p 32

# Genotyping was then conducted as per
bayesTyper genotype -v \
    -c \
    -s \
    -g \
    -d \
    -o \
    -z -p 32