#!/bin/bash -e
##################################################################################################################################
# Walkthrough of how I ran the BayesTyper software package. Turns out this has a lot of dependencies and steps. Below is an
# outline of how I filtered the raw Manta SV calls for genotyping with BayesTyper. With a Manta VCF, the first first step is to
# identify SV candidates by de novo assembly. This converts all symbolic alleles to a full sequence. The conversion results in
# large VCFs, so ensure that the output is gzipped iwth the -z flag. Note, bayesTyper does NOT support bgzip files. BayesTyper
# benefits from including SNP data as it helps with initialising the kmers. Therefore it's recommended to conduct SNP discovery
# ith either GATK HaplotypeCaller or FreeBayes. I used FreeBayes since it's fast. The third party program, KMC, and
# bayesTyperTools makeBloom are necessary to run BayesTyper cluster and BayesTyper genotype. It is independent of preparing the
# variant VCFs and can be done concurrrently. 
#################################################################################################################################
chr_ref=/kakapo-data/bwa/manta/bayestyper/kakapo_chromosome_ref.fa
decoy_ref=/kakapo-data/bwa/manta/bayestyper/kakapo_decoy_ref.fa
deepV=DEEPVARIANT:/kakapo-data/bwa/manta/bayestyper/Trained_SV_scaffolds_renamed.sorted.vcf
exclude=/kakapo-data/metadata/kakapo_SVexcluded_scaffolds.bed
female=/kakapo-data/bwa/bwa_female/nodup_bam/
male=/kakapo-data/bwa/bwa_male/nodup_bam/
out=/kakapo-data/bwa/manta/bayestyper/
raw_data=/kakapo-data/bwa/manta/raw_variants/
ref=/kakapo-data/References/kakapo_full_ref.fa
##################################################################################################################################
#  DID NOT USE FREEBAYES OUTPUTS IN THE END
#for i in {01..11}
#   do
#   for fbam in /kakapo-data/bwa/bwa_female/nodup/batch${i}/*_nodup.bam
#       do
#       fbase=$(basename $fbam _nodup.bam)
#        printf "\nRunning FreeBayes for ${fbase}...\n"
#        freebayes -f ${ref} ${fbam} > ${out}bayestyper/freebayes/${fbase}_freebayes.vcf
#        bcftools norm -f ${ref} -T ^${exclude} --check-ref we -m -any -O v -o ${out}bayestyper/freebayes/${fbase}_freebayes_norm.vcf ${out}bayestyper/freebayes/${fbase}_freebayes.vcf
#    done &
#done &
#for j in {01..14}
#    do
#    for mbam in /kakapo-data/bwa/bwa_male/nodup/batch${01..14}/*_nodup.bam
#    do
#        mbase=$(basename $mbam _nodup.bam)
#        printf "\nRunning FreeBayes for ${mbase}...\n"
#        freebayes -f $ref $mbam > /kakapo-data/bwa/manta/bayestyper/freebayes/${mbase}_freebayes.vcf
#        bcftools norm -f ${ref} -T ^${exclude} --check-ref we -m -any -O v -o ${out}bayestyper/freebayes/${mbase}_freebayes_norm.vcf ${out}bayestyper/freebayes/${mbase}_freebayes.vcf
#    done &
#done
#wait
##################################################################################################################################
# Merge sex-specific variant calls and batched calls
bcftools merge -m all -O z -o ${raw_data}joint_total.vcf.gz \
    ${raw_data}joint_calling/female_diploidSV.vcf.gz \
    ${raw_data}joint_calling/male_diploidSV.vcf.gz

bcftools merge -m all -O z -o ${raw_data}batch_total.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch01_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch02_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch03_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch04_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch05_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch06_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch07_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch01_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch02_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch03_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch04_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch05_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch06_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch07_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch08_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch09_manta_INV_conversion_male.vcf.gz

# Implementing SV filter for 'high-quality' SVs
bcftools view -i '(FILTER=="PASS" & FMT/PR >= 10 & FORMAT/FT == "PASS" & SVLEN<=-50 & SVTYPE=="DEL") | (FILTER=="PASS" & PR >= 10 & FMT/FT == "PASS" & SVLEN>=50 & SVTYPE!="BND")' \
    -T ^${exclude} \
    -O z -o ${out}joint_filtered/01_joint_filtered.SVs.vcf.gz \
    ${raw_data}joint_total.vcf.gz

bcftools view -i '(FILTER=="PASS" & FMT/PR >= 10 & FORMAT/FT == "PASS" & SVLEN<=-50 & SVTYPE=="DEL") | (FILTER=="PASS" & PR >= 10 & FMT/FT == "PASS" & SVLEN>=50 & SVTYPE!="BND")' \
    -T ^${exclude} \
    -O z -o ${out}batch_filtered/01_batch_filtered.SVs.vcf.gz \
    ${raw_data}batch_total.vcf.gz
##################################################################################################################################
# For including DeepVariant SNPs, all chromosome sizes were used to match kākāpō125+ and NCBI chromosome IDs. bcftools annotate
# was used to rename chromosomes to be consistent with NCBI scaffold names. SNPs were filtered down to only include regions of interest
# (i.e., scaffolds included for SV analysis). But after converting chromosome names, needed to sort VCF header and body of
# DeepVariant file to match the manta VCFs. Sorted the body with the solution below (https://www.biostars.org/p/84747/), then
# manually edited the header.
bcftools norm -c we -f ${ref} -O v -o ${out}bayestyper/Trained_SV_scaffolds_renamed_norm.vcf \
    ${out}bayestyper/Trained_SV_scaffolds_renamed.vcf
bgzip ${out}Trained_SV_scaffolds_renamed_norm.vcf; tabix -pvcf ${out}Trained_SV_scaffolds_renamed_norm.vcf.gz
cat chr_list.txt | xargs tabix -h ${out}Trained_SV_scaffolds_renamed_norm.vcf.gz > ${out}Trained_SV_scaffolds_renamed.sorted.vcf
##################################################################################################################################
# Convert Manta VCF to remove symbolic alleles
for manta in ${out}{batch,joint}_filtered/
    do
    file=$(basename ${manta})
    echo "Normalising 01_$file.SVs.vcf"
    bcftools norm --check-ref we -f ${ref} -m -any -O v -o ${manta}02_${file}_norm.vcf ${manta}01_$file.SVs.vcf
    echo "Running convertAllele for 02_${file}_norm.vcf..."
    bayesTyperTools convertAllele --variant-file ${manta}02_${file}_norm.vcf --genome-file ${ref} \
        --output-prefix ${manta}03_${file}_converted -z
    echo "Normalising 03_${file}_converted.vcf.gz"
    bcftools norm --threads 24 -f ${ref} -m -any -O z -o ${manta}04_${file}_converted_norm.vcf.gz \
        ${manta}03_${file}_converted.vcf.gz
    echo "Finally combining DeepVariant and Manta alleles into 05_${file}_mantaDV_candidates.vcf.gz"
    bayesTyperTools combine -v ${deepV},MANTA:${manta}04_${file}_converted_norm.vcf.gz \
        -o ${manta}05_${file}_mantaDV_candidates -z
done
# DID NOT USE FREEBAYES RUNS IN THE END.
# Merge freebayes and converted VCFs. Could not run these lines of code simultaneously. files overwrote one another despite being
# output to different files.
#free_vcfs=$(ls /kakapo-data/bwa/manta/bayestyper/freebayes/*_scaffolds_excluded_norm.vcf | sed 's%/kakapo-data/%FREEBAYES:/kakapo-data/%g' | xargs | tr " " ",")
##################################################################################################################################
# Running KMC and makeBloom
# It is really easy to over-resource KMC, that is that if you provide it with too
# much memory/threads the program freezes. I have had the best luck with 24Gb RAM and 48 threads as per below. It also found that
# I needed to run KMC one individual at a time. I had a lot of issues when I attempted to run multiple individuals at once since
# KMC uses the same temp file naming convention each time it runs. Basically temp files from different runs were overwriting each
# other and I needed to designate different locations for the temp files.
ulimit -n 2048
for mkmc in ${male}batch{01..14}/*_nodup.bam
    do
    indiv=$(basename ${mkmc} _nodup.bam)
    echo "Running KMC for ${indiv}..."
    kmc -k55 -ci1 -m24 -t48 -fbam ${mkmc} ${out}kmc/${indiv}_KMC.res ${out}kmc_tmp/
    echo "Running makeBloom for ${indiv}..."
    bayesTyperTools makeBloom -k ${out}kmc/${indiv}_KMC.res -p 8
done
for fkmc in ${female}batch{01..11}/*_nodup.bam
    do
    indiv=$(basename ${fkmc} _nodup.bam)
    echo "Running KMC for ${indiv}..."
    kmc -k55 -ci1 -m24 -t48 -fbam ${fkmc} ${out}kmc/${indiv}_KMC.res ${out}kmc_tmp/
    echo "Running makeBloom for ${indiv}..."
    bayesTyperTools makeBloom -k ${out}kmc/${indiv}_KMC.res -p 8
done
##################################################################################################################################
# For bayesTyper genotype, a ploidy file in the format below is required
# chr    female    male
# 1      2         2
# ...   ...       ...
# Z      1         2
# W      1         0
# Removed all unplaced scaffolds where ploidy was unknown in the candidate VCF and reference. Found a script (faSomeRecords) to
# remove the scaffolds where ploidy is unknown at: https://www.biostars.org/p/360658/ Then variant clusters were estimated.
# Note, no more than 30 indivs can be genotyped at once.
# The <samples>.tsv file should contain one sample per row with columns:
# <sample_id>, <sex> and <kmc_output_prefix> and no header.
for samps in ${out}sample_batches/sample_batch*.tsv
    do
    base=$(basename ${samps} .tsv)
    bayesTyper cluster --variant-file ${out}joint_filtered/05_joint_filtered_mantaDV_candidates.vcf.gz \
        --samples-file ${samps} --genome-file ${chr_ref} --decoy-file ${decoy_ref} \
        --output-prefix ${out}joint_filtered/${base}/ --threads 24
    bayesTyper cluster --variant-file ${out}batch_filtered/05_batch_filtered_mantaDV_candidates.vcf.gz \
        --samples-file ${samps} --genome-file ${chr_ref} --decoy-file ${decoy_ref} \
        --output-prefix ${out}batch_filtered/${base}/ --threads 24
done
for samps in ${out}sample_batches/sample_batch*.tsv
    do
    base=$(basename ${samps} .tsv)
    bayesTyper genotype --variant-clusters-file --cluster-data-dir --samples-file ${samps} \
        --genome-file $chr_ref --decoy-file $decoy_ref --output-prefix ${out}bayestyper/joint_filtered/05_${batch}_genotypes \
        --threads 24 --chromosome-ploidy-file ${out}bayestyper/ploidy.tsv --gzip-output
    bayesTyper genotype --variant-clusters-file --cluster-data-dir --samples-file ${samps} \
        --genome-file $chr_ref --decoy-file $decoy_ref --output-prefix ${out}bayestyper/batch_filtered/05_${batch}_genotypes \
        --threads 24 --chromosome-ploidy-file ${out}bayestyper/ploidy.tsv --gzip-output
done

##################################################################################################################################
# Lineage Comparisons
mkdir -p ${out}lineage_comparisons

bcftools view -s Richard_Henry ${out}bwa_delly_minsize_genotypes_filtered_trio.vcf | \
    bcftools view -i 'GT!="RR" & GT!="mis"' \
    -O z -o ${out}lineage_comparisons/intermediate_bcfs/RH_variants.vcf.gz

bcftools view -s ^Richard_Henry,Kuia,Gulliver,Sinbad,Adelaide,Henry,Marian,Gertrude ${out}bwa_delly_minsize_genotypes_filtered_trio.vcf | \
    bcftools view -i 'GT!="RR" & GT!="mis"' \
    -O z -o ${out}lineage_comparisons/intermediate_bcfs/SI_variants.vcf.gz

tabix ${out}lineage_comparisons/intermediate_bcfs/RH_variants.vcf.gz
tabix ${out}lineage_comparisons/intermediate_bcfs/SI_variants.vcf.gz

bcftools isec ${out}lineage_comparisons/intermediate_bcfs/RH_variants.vcf.gz \
    ${out}lineage_comparisons/intermediate_bcfs/SI_variants.vcf.gz \
    -p ${out}lineage_comparisons/