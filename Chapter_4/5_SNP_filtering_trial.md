# SNP data filtering trial
SNP filtering trial for tara iti and Australian fairy tern. Here we tested a number of outcomes for SNP filtering of tara iti specific, Australian fairy tern specific, and global SNPs. The steps taken here are adapted from [this](https://speciationgenomics.github.io/filtering_vcfs/) tutorial by Mark Ravinet & Joana Meier.

The Pixy software package requires the incorporation of invariant sites for estimates of pi, dxy and Fst. As such, a VCF with SNPs and a separate VCF with invariant sites were assessed as per below.

Global variables were fist defined as:
```
data=/data/common_tern/SNP_filtering_trial/
```
Then created initial VCF to explore filtering parameters for variant and invariant sites. As a first pass, all sites had to have <=20% missing data and all individuals had a minimum depth of 4x. The variant file only included biallelic SNPs. 
```
printf "\nCreating initial VCF for filtering trial...\n"
bcftools view -T ${data}chroms --threads 120 -i '(N_PASS(GT="mis")<=5) & (MEAN(FMT/DP>=4))' -m 2 -M 2 -v snps \
    -O z -o ${data}global_biallelic.vcf.gz \
    /data/common_tern/bcftools_variantCalls/Fairy_tern_AllCalls.sorted.vcf.gz

bcftools view -T ./chroms --threads 120 -i '(N_PASS(GT="mis")<=5) & (MEAN(FMT/DP=>4))' -c 0 -C 0 \
    -O z -o ${data}global_invariant.vcf.gz \
    /data/common_tern/bcftools_variantCalls/Fairy_tern_AllCalls.sorted.vcf.gz
```
This left 12,755,476 SNPs in the ```global_biallelic.vcf.gz``` file and 485,166,219 sites in the ```global_invariant.vcf.gz``` file. Prepared the SNP file by including HWE and MAF annotations as per:
```
bcftools +fill-tags ${data}SNP_filtering_trial/global_biallelic.vcf.gz -O z -o ${data}SNP_filtering_trial/global_biallelic_annotated.vcf.gz -- -t HWE,MAF
```

# Filtering parameters
Estimating allele frequency, mean site and individual depth, site quality, proportion of missing data per site and individual, heterozygosity and inbreeding per individual. 

```
annotated=${data}SNP_filtering_trial/global_biallelic_annotated.vcf.gz
invariant=${data}SNP_filtering_trial/global_invariant.vcf.gz
printf "\nNOW BEGINNING TO CALCULATE STATS...\n"
for vcf in ${annotated} ${invariant}
  do
  base=$(basename ${vcf} .vcf.gz)
    if (( "$base" == global_biallelic_annotated ))
      then
      echo "Calculating stats for ${base}..."
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --freq2 &
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --site-mean-depth &
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --depth &
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --site-quality &
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --missing-site &
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --missing-indv &
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --het
    else
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --site-mean-depth &
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --depth &
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --site-quality &
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --missing-site &
      vcftools --gzvcf ${vcf} --out ${data}total/${base} --missing-indv
    fi &
done
```
# Visualing the outputs
To visualise and explore the characteristics of variant and invariant sites, an R script plotting the outputs of VCFtools was as follows:
```
library(tidyverse)
setwd("/data/common_tern/SNP_filtering_trial/")

var_qual <- read_delim("total/global_variant_autosomes.lqual", delim = "\t", 
                   col_names = c("chr", "pos", "qual"), skip = 1)
var_depth <- read_delim("total/global_variant_autosomes.ldepth.mean", delim = "\t", 
                    col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
var_miss <- read_delim("total/global_variant_autosomes.lmiss", delim = "\t", 
                   col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
var_freq <- read_delim("total/global_variant_autosomes.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_ind_depth <- read_delim("total/global_variant.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
var_ind_miss  <- read_delim("total/global_variant.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
var_ind_het <- read_delim("total/global_variant.het", delim = "\t",
           col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)

fix_qual <- read_delim("total/global_invariant_autosomes.lqual", delim = "\t", 
                   col_names = c("chr", "pos", "qual"), skip = 1)
fix_depth <- read_delim("total/global_invariant_autosomes.ldepth.mean", delim = "\t", 
                    col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
fix_miss <- read_delim("total/global_invariant_autosomes.lmiss", delim = "\t", 
                   col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
fix_ind_depth <- read_delim("total/global_invariant.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
fix_ind_miss  <- read_delim("total/global_invariant.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Plot Variant Quality 
qual_plot <- ggplot(var_qual, aes(qual)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(title = "Variant Quality Distribution") +
  theme_light()

# Plot Variant Mean Depth
depth_plot <- ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(title = "Variant Depth Distribution") +
  xlim(0, 50) +
  theme_light()

fix_qual_plot <- ggplot(fix_qual, aes(qual)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(title = "Variant Quality Distribution") +
  theme_light()

# Plot Variant Mean Depth
fix_depth_plot <- ggplot(fix_depth, aes(mean_depth)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(title = "Variant Depth Distribution") +
  xlim(0, 50) +
  theme_light()

print(summary(var_depth$mean_depth))
print(summary(fix_depth$mean_depth))

# Variant missingness
miss_plot <- ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(title = "Variant Missingness Distribution") +
  theme_light()

  fix_miss_plot <- ggplot(fix_miss, aes(fmiss)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(title = "Variant Missingness Distribution") +
  theme_light()

print(summary(var_miss$fmiss))
print(summary(fix_miss$fmiss))

# Minor allele frequency
var_freq$maf <- var_freq %>% 
  select(a1, a2) %>%
  apply(1, function(z) min(z))

maf_plot <- ggplot(var_freq, aes(maf)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(title = "MAF distribution") +
  theme_light()

print(summary(var_freq$maf))

# Saving to PDF
pdf("Variant_characteristics.pdf")
print(qual_plot)
print(depth_plot)
print(miss_plot)
print(maf_plot)
dev.off() 

pdf("Invariant_characteristics.pdf")
print(fix_qual_plot)
print(fix_depth_plot)
print(fix_miss_plot)
dev.off() 
```
# Final SNP Filtering

```
bcftools view --threads 64 \
  -i 'MQ>=20 & (MIN(FMT/DP)>=5) & (MAX(FMT/DP)<=25) & (MEAN(FMT/DP)>=5) & (MEAN(FMT/DP)<=25) & (MAF>=0.05)' \
  ${data}SNP_filtering_trial/global_variant_autosomes.vcf.gz | \
  bcftools +prune -m 0.8 -w 1000 | \
  bcftools view -i 'N_PASS(GT="mis")=0' -O z -o ${data}SNP_filtering_trial/global_MAF.vcf.gz

bcftools view --threads 64 -i 'MQ>=20 & (MIN(FMT/DP)>=5) & (MAX(FMT/DP)<=25) & (MEAN(FMT/DP)>=5) & (MEAN(FMT/DP)<=25) & (MAF>=0.05)' \
  -O z -o global_MAF_noPrune.vcf.gz \
  global_variant_autosomes_annotated.vcf.gz


bcftools view --threads 64 \
  -i '(MQ>=20) & (MEAN(FMT/DP)>=5) & (MEAN(FMT/DP)<=25) & (MIN(FMT/DP)>=5) & (MAX(FMT/DP)<=25) & (N_PASS(GT="mis")=0)' \
  -O z -o /data/common_tern/SNP_filtering_trial/global_invariant_filtered2.vcf.gz /data/common_tern/SNP_filtering_trial/ \
  global_invariant_autosomes.vcf.gz

bcftools view -i 'N_PASS(GT="mis")=0' -O z -o ${data}global_MAFfiltered_LD_nomiss.vcf.gz ${data}global_variant_MAFfiltered_LD.vcf.gz
bcftools view -i 'N_PASS(GT="mis")=0' -O z -o ${data}global_noMAFfiltered_LD_nomiss.vcf.gz ${data}global_variant_noMAFfiltered_LD.vcf.gz
```
The ```global_MAF_noPrune.vcf.gz``` SNPs were used for counting the total number of variable SNPs, private SNPs, fixed SNPs, and the proportion of variable to invariant sites. These data as well as the pruned SNPs were used for estimates of pi, Dxy, and Fst with pixy, estimates of heterozygosity with VCFtools. Finally, admixture was run using the pruned SNPs. 

## Visualising filtered outputs
After filtering, the estimates of 

```
for vcf in ${data}global_MAF_noPrune.vcf.gz
  do
  vcftools --gzvcf ${vcf} --out ${data}filter_stats/noLD_filter_SNPs --freq2 &
  vcftools --gzvcf ${vcf} --out ${data}filter_stats/noLD_filter_SNPs --site-mean-depth &
  vcftools --gzvcf ${vcf} --out ${data}filter_stats/noLD_filter_SNPs --depth &
  vcftools --gzvcf ${vcf} --out ${data}filter_stats/noLD_filter_SNPs --site-quality
done

for vcf in ${data}global_MAF.vcf.gz
  do
  vcftools --gzvcf ${vcf} --out ${data}filter_stats/LD_filtered_SNPs --freq2 &
  vcftools --gzvcf ${vcf} --out ${data}filter_stats/LD_filtered_SNPs --site-mean-depth &
  vcftools --gzvcf ${vcf} --out ${data}filter_stats/LD_filtered_SNPs --depth &
  vcftools --gzvcf ${vcf} --out ${data}filter_stats/LD_filtered_SNPs --site-quality
done
```