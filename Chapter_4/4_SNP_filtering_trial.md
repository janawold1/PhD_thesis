# SNP data filtering trial
SNP filtering trial for tara iti and Australian fairy tern. Here we tested a number of outcomes for SNP filtering of tara iti specific, Australian fairy tern specific, and global SNPs. The steps taken here are adapted from (this)[https://speciationgenomics.github.io/filtering_vcfs/] tutorial by Mark Ravinet & Joana Meier.

The Pixy software package requires the incorporation of invariant sites for estimates of pi, dxy and Fst. As such, a VCF with SNPs and a separate VCF with invariant sites were assessed as per below.

Global variables were fist defined as:
```
data=/data/common_tern/SNP_filtering_trial/
```
Then created initial VCF for variant and invariant sites, excluding indels, and any sites with more than 2 alleles:
```
printf "\nCreating initial VCF for filtering trial...\n"
bcftools view --threads 24 -m 2 -M 2 -v snps \
    -O z -o ${data}global_variant.vcf.gz \
    /data/common_tern/bcftools_variantCalls/Fairy_tern_AllCalls.sorted.vcf.gz

bcftools view --threads 24 -M 0 \
    -O z -o ${data}global_invariant.vcf.gz \
    /data/common_tern/bcftools_variantCalls/Fairy_tern_AllCalls.sorted.vcf.gz
```
This left 33,188,437 SNPs in the ```global_variant.vcf.gz``` file and 39,109,127 sites in the ```global_invariant.vcf.gz``` file.

# Filtering parameters
Estimating allele frequency, mean site and individual depth, site quality, proportion of missing data per site and individual, heterozygosity and inbreeding per individual. 

```
printf "\nNOW BEGINNING TO CALCULATE STATS...\n"
for vcf in ${data}*.vcf.gz
do
    base=$(basename ${vcf} .vcf.gz)
    echo "Calculating stats for ${base}..."
    vcftools --gzvcf ${vcf} --out ${data}stats/${base} --freq2 &
    vcftools --gzvcf ${vcf} --out ${data}stats/${base} --site-mean-depth &
    vcftools --gzvcf ${vcf} --out ${data}stats/${base} --depth &
    vcftools --gzvcf ${vcf} --out ${data}stats/${base} --site-quality &
    vcftools --gzvcf ${vcf} --out ${data}stats/${base} --missing-site &
    vcftools --gzvcf ${vcf} --out ${data}stats/${base} --missing-indv &
    vcftools --gzvcf ${vcf} --out ${data}stats/${base} --het
done
```