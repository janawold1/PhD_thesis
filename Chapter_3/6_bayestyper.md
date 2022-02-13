# Runthrough of the BayesTyper software package
BayesTyper is a powerful baysian based genotyping program. However, one draw back is that it is computationally intensive and has multiple intermediate steps to navigate. 

Below I outline of how I filtered the raw Manta SV calls for genotyping with BayesTyper, create kmer counts, cluster variant kmers and finally genotype variants.  

To begin, global variables were defined as below:
```
chr_ref=/kakapo-data/bwa/manta/bayestyper/kakapo_chromosome_ref.fa
decoy_ref=/kakapo-data/bwa/manta/bayestyper/kakapo_decoy_ref.fa
deepV=DEEPVARIANT:/kakapo-data/bwa/manta/bayestyper/Trained_SV_scaffolds_renamed.sorted.vcf
exclude=/kakapo-data/metadata/kakapo_SVexcluded_scaffolds.bed
female=/kakapo-data/bwa/bwa_female/nodup_bam/
male=/kakapo-data/bwa/bwa_male/nodup_bam/
out=/kakapo-data/bwa/manta/bayestyper/
raw_data=/kakapo-data/bwa/manta/raw_variants/
ref=/kakapo-data/References/kakapo_full_ref.fa
trio=/kakapo-data/metadata/sample_trios.csv
````

## Merging sex-specific variant calls and batched calls for SV genotyping
Candidate SVs were called for two Manta data sets 1) the Joint data set; and 2) the Batched data set. 

```
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
```
## Implementing SV filters for 'high-quality' SVs
```
bcftools view -i '(FILTER=="PASS" & FMT/PR >= 10 & FORMAT/FT == "PASS" & SVLEN<=-50 & SVTYPE=="DEL") | (FILTER=="PASS" & PR >= 10 & FMT/FT == "PASS" & SVLEN>=50 & SVTYPE!="BND")' \
    -T ^${exclude} \
    -O z -o ${out}joint_filtered/01_joint_filtered.SVs.vcf.gz \
    ${raw_data}joint_total.vcf.gz

bcftools view -i '(FILTER=="PASS" & FMT/PR >= 10 & FORMAT/FT == "PASS" & SVLEN<=-50 & SVTYPE=="DEL") | (FILTER=="PASS" & PR >= 10 & FMT/FT == "PASS" & SVLEN>=50 & SVTYPE!="BND")' \
    -T ^${exclude} \
    -O z -o ${out}batch_filtered/01_batch_filtered.SVs.vcf.gz \
    ${raw_data}batch_total.vcf.gz
```
For including DeepVariant SNPs, all chromosome sizes were used to match kākāpō125+ and NCBI chromosome IDs. bcftools annotate was used to rename chromosomes to be consistent with NCBI scaffold names. SNPs were filtered down to only include regions of interest (i.e., scaffolds included for SV analysis). But after converting chromosome names, needed to sort VCF header and body of
 DeepVariant file to match the manta VCFs. Sorted the body with [this solution](https://www.biostars.org/p/84747/), then manually edited the header.
 ```
bcftools norm -c we -f ${ref} -O v -o ${out}bayestyper/Trained_SV_scaffolds_renamed_norm.vcf \
    ${out}bayestyper/Trained_SV_scaffolds_renamed.vcf
bgzip ${out}Trained_SV_scaffolds_renamed_norm.vcf; tabix -pvcf ${out}Trained_SV_scaffolds_renamed_norm.vcf.gz
cat chr_list.txt | xargs tabix -h ${out}Trained_SV_scaffolds_renamed_norm.vcf.gz > ${out}Trained_SV_scaffolds_renamed.sorted.vcf
```
## Convert Manta VCF to remove symbolic alleles
This converts all symbolic alleles to a full sequence. The conversion results in large VCFs, so ensure that the output is gzipped iwth the -z flag. Note, bayesTyper does NOT support bgzip files. BayesTyper benefits from including SNP data as it helps with initialising the kmers. Therefore it's recommended to conduct SNP discovery with either GATK HaplotypeCaller or FreeBayes. I used FreeBayes since it's fast. The third party program, KMC, and bayesTyperTools makeBloom are necessary to run BayesTyper cluster and BayesTyper genotype. It is independent of preparing the variant VCFs and can be done concurrrently.
```
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
```
## Running KMC and makeBloom
It is really easy to over-resource KMC, that is that if you provide it with too
much memory/threads the program freezes. I have had the best luck with 24Gb RAM and 48 threads as per below. It also found that I needed to run KMC one individual at a time. I had a lot of issues when I attempted to run multiple individuals at once since KMC uses the same temp file naming convention each time it runs. Basically temp files from different runs were overwriting each other and I needed to designate different locations for the temp files.
```
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
```
For bayesTyper genotype, a ploidy file in the format below is required

|  chr  |  female  |  male  |
|:-----:|:--------:|:------:|
|   1   |     2    |    2   |
|  ...  |    ...   |   ...  |
|   Z   |     1    |    2   |
|   W   |     1    |    0   |

Removed all unplaced scaffolds where ploidy was unknown in the candidate VCF and reference. Found a script (faSomeRecords) to remove the scaffolds where ploidy is unknown at: https://www.biostars.org/p/360658/ Then variant clusters were estimated. Note, no more than 30 indivs can be genotyped at once. The <samples>.tsv file should contain one sample per row with columns:

<sample_id>, <sex> and <kmc_output_prefix> and no header.

```
for samps in ${out}sample_batches/sample_batch*.tsv
    do
    base=$(basename ${samps} .tsv)
    bayesTyper cluster --variant-file ${out}joint_filtered/05_joint_filtered_mantaDV_candidates.vcf.gz \
        --samples-file ${samps} --genome-file ${chr_ref} --decoy-file ${decoy_ref} \
        --output-prefix ${out}joint_filtered/${base}/ --threads 24 &
done
wait
for samps in ${out}sample_batches/sample_batch*.tsv
    do
    base=$(basename ${samps} .tsv)
    bayesTyper cluster --variant-file ${out}batch_filtered/05_batch_filtered_mantaDV_candidates.vcf.gz \
        --samples-file ${samps} --genome-file ${chr_ref} --decoy-file ${decoy_ref} \
        --output-prefix ${out}batch_filtered/${base}/ --threads 24 &
done
wait
for samps in ${out}bayestyper/joint_filtered/sample_batch*
    do
    base=$(basename ${samps})
    bayesTyper genotype --variant-clusters-file ${samps}/${base}_unit_1/variant_clusters.bin \
        --cluster-data-dir ${samps}/${base}_cluster_data/ \
        --samples-file ${out}bayestyper/sample_batches/${base}.tsv \
        --genome-file $chr_ref --decoy-file $decoy_ref \
        --output-prefix ${samps}/${base}_genotypes \
        --threads 24 --chromosome-ploidy-file ${out}bayestyper/ploidy.tsv --gzip-output
done
for samps in ${out}bayestyper/batch_filtered/sample_batch*
    do
    base=$(basename ${samps})
    bayesTyper genotype --variant-clusters-file ${samps}/${base}_unit_1/variant_clusters.bin \
        --cluster-data-dir ${samps}/${base}_cluster_data/ \
        --samples-file ${out}bayestyper/sample_batches/${base}.tsv \
        --genome-file $chr_ref --decoy-file $decoy_ref \
        --output-prefix ${samps}/${base}_genotypes \
        --threads 24 --chromosome-ploidy-file ${out}bayestyper/ploidy.tsv --gzip-output
done
```
## Merging SV genotype batches
```
for vcf in ${out}bayestyper/joint_filtered/sample_batch{1..6}/*_genotypes.vcf.gz
    do
    batch=$(basename $vcf _genotypes.vcf.gz)
    echo "Converting $vcf to bgzip compression..."
    gunzip $vcf
    bgzip ${out}bayestyper/joint_filtered/${batch}/${batch}_genotypes.vcf
    echo "Now indexing $vcf"
    tabix ${out}bayestyper/joint_filtered/${batch}/${batch}_genotypes.vcf.gz
done

bcftools merge -m id -O z -o joint_filtered/06_manta_genotypes.vcf.gz \
    joint_filtered/sample_batch1/sample_batch1_genotypes.vcf.gz \
    joint_filtered/sample_batch2/sample_batch2_genotypes.vcf.gz \
    joint_filtered/sample_batch3/sample_batch3_genotypes.vcf.gz \
    joint_filtered/sample_batch4/sample_batch4_genotypes.vcf.gz \
    joint_filtered/sample_batch5/sample_batch5_genotypes.vcf.gz \
    joint_filtered/sample_batch6/sample_batch6_genotypes.vcf.gz

 bcftools view --threads 24 -i '(ACO="MANTA") & ((N_PASS(GT=="mis") < 17) & (N_PASS(GT="alt")>=1))' \
    -O v -o joint_filtered/07_joint_filtered_genotypes.vcf \
    joint_filtered/06_joint_genotypes.vcf.gz

for vcf in ${out}bayestyper/batch_filtered/sample_batch{1..6}/*_genotypes.vcf.gz
    do
    batch=$(basename $vcf _genotypes.vcf.gz)
    echo "Converting $vcf to bgzip compression..."
    gunzip $vcf
    bgzip ${out}bayestyper/batch_filtered/${batch}/${batch}_genotypes.vcf
    echo "Now indexing $vcf"
    tabix ${out}bayestyper/batch_filtered/${batch}/${batch}_genotypes.vcf.gz
done

bcftools merge -m id -O z -o batch_filtered/06_manta_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch1/sample_batch1_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch2/sample_batch2_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch3/sample_batch3_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch4/sample_batch4_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch5/sample_batch5_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch6/sample_batch6_genotypes.vcf.gz

 bcftools view --threads 24 -i '(ACO="MANTA") & ((N_PASS(GT=="mis") < 17) & (N_PASS(GT="alt")>=1))' \
    -O v -o ${out}bayestyper/batch_filtered/07_batch_filtered_genotypes.vcf \
    ${out}bayestyper/batch_filtered/06_batch_genotypes.vcf.gz

```

## Finding SV type

BayesTyper removes symbolic alleles from the genotype files. To make comparisons among SV tools similar, the SV type called by Manta was resoved as per below.
1) First identify the locations of genotyped variants:
```
bcftools query -f '%CHROM\t%POS\n' ${out}bayestyper/batch_filtered/08_batch_filtered_trios.vcf > ${out}bayestyper/batch_genos
bcftools query -f '%CHROM\t%POS\n' ${out}bayestyper/joint_filtered/08_joint_filtered_trios.vcf > ${out}bayestyper/joint_genos
```

2) Count overlaping SVtypes:
```
bcftools query -T ${out}bayestyper/batch_genos -f '%SVTYPE\n' \
    ${out}bayestyper/batch_filtered/02_batch_filtered_norm.vcf.gz | sort | uniq -c

bcftools query -T ${out}bayestyper/joint_genos -f '%SVTYPE\n' \
    ${out}bayestyper/joint_filtered/02_joint_filtered_norm.vcf.gz | sort | uniq -c
```

Found that 76 SVs don't overlap in the batch data and 126 SVs don't overlap in the joint data. 

3) Find non-overlapping sites (i.e., those present in genotyped output, but absent in call set):
```

```

4) Annotate VCF for individual counts:
First prepare the annotation file by identifying the resolvable SVs and their types:
```
bcftools query -f '%CHROM\t%POS\n' batch_filtered/07_batch_filtered_genotypes.vcf > batch_geno_sites
bcftools query -T batch_geno_sites -f '%CHROM\t%POS\t%SVTYPE\t%SVLEN\n' batch_filtered/02_batch_filtered_norm.vcf.gz > batch_conversion

bcftools query -f '%CHROM\t%POS\n' joint_filtered/07_joint_filtered_genotypes.vcf > joint_geno_sites
bcftools query -T joint_geno_sites -f '%CHROM\t%POS\t%SVTYPE\t%SVLEN\n' joint_filtered/02_joint_filtered_norm.vcf.gz > joint_conversion
```

Then edit the ```batch_conversion``` file so the first line in this file contains ```#CHROM POS SVTYPE  SVLEN``` with nano. 

And create a file that captures the information to be appended into the header of the VCF.

Contents of annots.hdr for example:
```
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant called by Manta">
##INFO=<ID=SVLEN,Number=1,Type=String,Description="Length of structural variant called by Manta">
```
Now time to annotate the file to continue with the genotype based analyses.
```
bgzip batch_conversion
bgzip joint_conversion

tabix -s1 -b2 -e2 batch_conversion.gz
tabix -s1 -b2 -e2 joint_conversion.gz

bcftools annotate -a batch_conversion.gz -h annots.hdr \
    -c CHROM,POS,SVTYPE,SVLEN -O v -o 09_batch_annotated.vcf \
    ${out}bayestyper/batch_filtered/07_batch_filtered_genotypes.vcf

bcftools annotate -a joint_conversion.gz -h annots.hdr \
    -c CHROM,POS,SVTYPE,SVLEN -O v -o 09_joint_annotated.vcf \
    ${out}bayestyper/joint_filtered/07_joint_filtered_genotypes.vcf
```

## Mendelian Inheritance Tests
Tests of Mendelian Inheritance were conducted as per: 

```
bcftools +mendelian -m a -T ${trio} -O v -o ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf \
    ${out}bayestyper/joint_filtered/08_joint_annotated.vcf
bcftools +mendelian -m a -T ${trio} -O v -o ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf \
    ${out}bayestyper/batch_filtered/08_batch_annotated.vcf

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.05' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.05_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.1' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.1_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.2' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.2_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv



bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.05' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.05_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.1' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.1_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.2' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.2_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv
```
## Summarising the number of SVs carried by individuals

```
while read -r line
    do
    indiv=$(echo $line | awk '{print $1}')
    gen=$(echo $line | awk '{print $2}')
    echo "Counting SVs for $indiv..."
    bcftools view -s ${indiv} ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\n' | awk -v var="$gen" '{print $0"\t"var}' >> ${out}bayestyper/summary/manta_batch_generations.tsv
    bcftools view -s ${indiv} ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\n' | awk -v var="$gen" '{print $0"\t"var}' >> ${out}bayestyper/summary/manta_joint_generations.tsv
done < /kakapo-data/metadata/generations.tsv
```


## Lineage Comparisons
```
mkdir -p ${out}bayestyper/lineage_{batch,joint}_comparisons

bcftools view -s Richard_Henry ${out}bayestyper/joint_filtered/09_joint_filtered_trios_annotated.vcf | \
    bcftools view -i 'GT!="RR" & GT!="mis"' -O z -o ${out}bayestyper/lineage_joint_comparisons/RH_variants.vcf.gz

bcftools view -s ^Richard_Henry,Kuia,Gulliver,Sinbad,Adelaide,Henry,Marian,Gertrude \
    ${out}/bayestyper/joint_filtered/09_joint_filtered_trios_annotated.vcf | bcftools view -i 'GT!="RR" & GT!="mis"' \
    -O z -o ${out}bayestyper/lineage_joint_comparisons/SI_variants.vcf.gz

tabix ${out}bayestyper/lineage_joint_comparisons/RH_variants.vcf.gz
tabix ${out}bayestyper/lineage_joint_comparisons/SI_variants.vcf.gz

bcftools isec ${out}bayestyper/lineage_joint_comparisons/RH_variants.vcf.gz \
    ${out}bayestyper/lineage_joint_comparisons/SI_variants.vcf.gz \
    -p ${out}bayestyper/lineage_joint_comparisons/
```