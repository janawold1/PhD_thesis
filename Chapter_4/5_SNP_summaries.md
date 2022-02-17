# Overview of Running Pixy
(pixy)[https://pixy.readthedocs.io/en/latest/about.html] is a tool for unbiased estimates of nucleotide diversity from a VCF. It uses both variable and invariant sites for these estimates. 

To begin, SNPs were filtered as per (4_SNP_filtering_trial.md). Filtered biallelic SNPs and filtered invariant sites were merged as per:
```
data=/data/common_tern/

mkdir ${data}pixy

bcftools concat --threads 64 --allow-overlaps \
    ${data}SNP_filtering_trial/global_MAF.vcf.gz \
    ${data}SNP_filtering_trial/global_invariant_filtered.vcf.gz \
    -O z -o ${data}pixy/total_MAF.vcf.gz

bcftools concat --threads 64 --allow-overlaps \
    ${data}SNP_filtering_trial/global_noMAF.vcf.gz \
    ${data}SNP_filtering_trial/global_invariant_filtered.vcf.gz \
    -O z -o ${data}pixy/total_noMAF.vcf.gz

tabix ${data}pixy/total_MAF.vcf.gz
tabix ${data}pixy/total_noMAF.vcf.gz
```

A tab delimited file containing sample name and population was created, and pixy run as per:
```
pixy --stats pi fst dxy \
    --vcf ${data}pixy/total_MAFfiltered.vcf.gz \
    --populations ${data}pixy/pop_map \
    --window_size 10000 \
    --output_folder ${data}pixy \
    --output_prefix total_MAF \
    --n_cores 24

pixy --stats pi fst dxy \
    --vcf ${data}pixy/total_noMAFfiltered.vcf.gz \
    --populations ${data}pixy/pop_map \
    --window_size 10000 \
    --output_folder ${data}pixy \
    --output_prefix total_noMAF \
    --n_cores 24

sed -i 's/Super_scaffold_//g' total_*
```

## Running ADMIXTURE
To run ADMIXTURE, the VCFs must be converted to plink format (I chosed the binary .bed format). For this Chromosome scaffold names were converted as per:
```
gunzip -c ${data}SNP_filtering_trial/global_MAF.vcf.gz | sed 's/Super_scaffold_//g' | grep -v "arrow_ctg1" > ${data}admixture/global_MAF.vcf

gunzip -c ${data}SNP_filtering_trial/global_noMAF.vcf.gz | sed 's/Super_scaffold_//g' | grep -v "arrow_ctg1" > ${data}admixture/global_noMAF.vcf
```
Then VCFs converted to plink format and admixture Run:
```
plink --vcf ${data}admixture/global_noMAF.vcf --recode12 --out admixture/noMAF/global_noMAF

plink --vcf ${data}admixture/global_MAF.vcf --recode12 --out admixture/MAF/global_MAF


for K in 1 2 3 4 5
    do
    admixture --cv ${data}admixture/noMAF/global_noMAF.ped $K | tee log${K}.out
    admixture --cv ${data}admixture/MAF/global_MAF.ped $K | tee log${K}.out
done
```
Cross-validation error was found with:
```
grep -h CV ${data}admixture/MAF/log*.out
grep -h CV ${data}admixture/noMAF/log*.out

```

Figures were plotted in R as per:
```
setwd("/data/common_tern/admixture/")

maf <- read.table("MAF/global_MAF.2.Q")
nomaf <- read.table("noMAF/global_noMAF.2.Q")

pdf("MAF_admixture_plot.pdf")
barplot(t(as.matrix(maf)), col=rainbow(3), xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()

pdf("noMAF_admixture_plot.pdf")
barplot(t(as.matrix(nomaf)), col=rainbow(3), xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()
```

## SNP summaries
Total SNPs, invariant sites, and allele counts per population:
```
bcftools view -S ${data}AU_samples.tsv global_MAF.vcf.gz | bcftools view -H -i 'GT!="RR"' | wc -l # Total SNPs in AU in the SNP file
bcftools view -S ${data}TI_samples.tsv global_MAF.vcf.gz | bcftools view -H -i 'GT!="RR"' | wc -l # Total SNPs in TI in the SNP file

bcftools view -S ${data}AU_samples.tsv global_MAF.vcf.gz | bcftools view -H -i 'N_PASS(GT=="RR")=15' | wc -l # Invariant sites in AU in the SNP file
bcftools view -S ${data}TI_samples.tsv global_MAF.vcf.gz | bcftools view -H -i 'N_PASS(GT=="RR")=11' | wc -l # Invariant sites in TI in the SNP file

bcftools view -S ${data}AU_samples.tsv global_MAF.vcf.gz | bcftools query -f '%AC\n' | awk '{sum += $1}; END {print sum}' # Total AU allele count
bcftools view -S ${data}TI_samples.tsv global_MAF.vcf.gz | bcftools query -f '%AC\n' | awk '{sum += $1}; END {print sum}' # Total TI allele count
```
Private alleles were estimated with the poppr library in R.

Mean Depth per SNP per Individual
```
while read -r line
    do
    echo "Averaging depth for $line..."
    dp=$(bcftools view -s $line ${data}final_SNPs/global_MAF.vcf.gz | bcftools query -f '[%DP]\n' | awk '{sum += $1}; END {print sum/NR}')
    echo "${line},${dp} >> indiv_dp.csv
done < samples.txt
```
Mean SNP depth
```
```
Mean SNPs per Kb
```
```
Ts/Tv Ratios
```
bcftools view -T ${data}metadata/common_tern_macrochromosomes -S ${data}AU_samples.tsv ${data}global_MAF.vcf.gz |\
    bcftools view -e 'N_PASS(GT=="RR")=15' | bcftools stats | grep "TSTV"
bcftools view -T ${data}metadata/common_tern_macrochromosomes -S ${data}TI_samples.tsv ${data}global_MAF.vcf.gz |\
    bcftools view -e 'N_PASS(GT=="RR")=11' | bcftools stats | grep "TSTV"
bcftools view -T ${data}metadata/common_tern_macrochromosomes ${data}global_MAF.vcf.gz | bcftools stats | grep "TSTV"

bcftools view -T ${data}metadata/common_tern_microchromosomes -S ${data}AU_samples.tsv ${data}global_MAF.vcf.gz |\
    bcftools view -e 'N_PASS(GT=="RR")=15' | bcftools stats | grep "TSTV"
bcftools view -T ${data}metadata/common_tern_microchromosomes -S ${data}TI_samples.tsv ${data}global_MAF.vcf.gz |\
    bcftools view -e 'N_PASS(GT=="RR")=11' | bcftools stats | grep "TSTV"
bcftools view -T ${data}metadata/common_tern_microchromosomes ${data}global_MAF.vcf.gz | bcftools stats | grep "TSTV"
```