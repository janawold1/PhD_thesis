# Overview of Running Pixy
[pixy](https://pixy.readthedocs.io/en/latest/about.html) is a tool for unbiased estimates of nucleotide diversity from a VCF. It uses both variable and invariant sites for these estimates. 

To begin, SNPs were filtered as per [5_SNP_filtering_trial.md](https://github.com/janawold1/PhD_thesis/blob/main/Chapter_4/5_SNP_filtering_trial.md). Filtered biallelic SNPs and filtered invariant sites were merged as per:
```
data=/data/common_tern/

mkdir ${data}pixy

bcftools concat --threads 64 --allow-overlaps \
    ${data}SNP_filtering_trial/global_MAF.vcf.gz \
    ${data}SNP_filtering_trial/global_invariant_filtered.vcf.gz \
    -O z -o ${data}pixy/total_MAF.vcf.gz

bcftools concat --threads 64 --allow-overlaps \
    ${data}SNP_filtering_trial/global_MAF_noPrune.vcf.gz \
    ${data}SNP_filtering_trial/global_invariant_filtered.vcf.gz \
    -O z -o ${data}pixy/total_MAF_noPrune.vcf.gz

tabix ${data}pixy/total_MAF.vcf.gz
tabix ${data}pixy/total_MAF_noPrune.vcf.gz
```

A tab delimited file containing sample name and population was created, and pixy run as per:
```
for i in {100000,10000}
    do
    for vcf in total_MAF total_MAF_noPrune
        do
        pixy --stats pi fst dxy \
            --vcf ${data}pixy/${vcf}filtered.vcf.gz \
            --populations ${data}pixy/pop_map \
            --window_size ${i} \
            --output_folder ${data}pixy \
            --output_prefix ${vcf}_${i} \
            --n_cores 24
    done
done

sed -i 's/Super_scaffold_//g' total_*
```
### Post-hoc aggregating
As recommended by the [pixy documentation](https://pixy.readthedocs.io/en/latest/output.html) the sum of the raw counts were used to recompute the differences/comparisons ratios rather than averaging the summary statistics themselves as per:
```
grep AU ${data}pixy/total_MAF_pi.txt | awk '{diff+=$7}; {comp+=$8}; END {print diff/comp}'
grep TI ${data}pixy/total_MAF_pi.txt | awk '{diff+=$7}; {comp+=$8}; END {print diff/comp}'
```
## Running ADMIXTURE
To run ADMIXTURE, the VCFs must be converted to plink format (I chosed the binary .bed format). For this Chromosome scaffold names were converted as per:
```
gunzip -c ${data}SNP_filtering_trial/global_MAF.vcf.gz | sed 's/Super_scaffold_//g' > ${data}admixture/global_MAF.vcf

gunzip -c ${data}SNP_filtering_trial/global_noMAF.vcf.gz | sed 's/Super_scaffold_//g' > ${data}admixture/global_MAF_noPrune.vcf
```
Then VCFs converted to plink format and admixture Run:
```
plink --vcf ${data}admixture/global_MAF.vcf --recode12 --out ${data}admixture/MAF/global_MAF
plink --vcf ${data}admixture/global_MAF_noPrune.vcf --recode12 --out ${data}admixture/noPruned_MAF/global_MAF_noPrune

for K in 1 2 3 4 5
    do
    admixture --cv ${data}admixture/MAF/global_MAF.ped $K | tee log${K}.out
    admixture --cv ${data}admixture/noPruned_MAF/global_MAF_noPrune.ped $K | tee log${K}.out
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
Total variable, private and fixed SNPs as well as to number of invariant sites per population and those shared between the two:
```
bcftools view -S ${data}AU_samples.tsv global_MAF_noPrune.vcf.gz | bcftools view -H -i 'GT!="alt"' | wc -l # Total SNPs in AU in the SNP file
bcftools view -S ${data}TI_samples.tsv global_MAF_noPrune.vcf.gz | bcftools view -H -i 'GT!="alt"' | wc -l # Total SNPs in TI in the SNP file

bcftools view -S ${data}AU_samples.tsv global_MAF_noPrune.vcf.gz | bcftools view -H -i 'N_PASS(GT=="RR")=15' | wc -l # Total fixed SNPs in AU
bcftools view -S ${data}TI_samples.tsv global_MAF_noPrune.vcf.gz | bcftools view -H -i 'N_PASS(GT=="RR")=11' | wc -l # Total fixed SNPs in TI
```
Finding the number of SNPs private to each group.
```
bcftools view --threads 64 -S AU_samples.tsv final_SNPs/global_MAF_noPrune.vcf.gz |\
  bcftools view --threads 64 -i 'GT=="alt"' -O z -o final_SNPs/AU_MAF_noPrune.vcf.gz
bcftools view --threads 64 -S TI_samples.tsv final_SNPs/global_MAF_noPrune.vcf.gz |\
  bcftools view -i 'GT=="alt"' -O z -o final_SNPs/TI_MAF_noPrune.vcf.gz

bcftools index final_SNPs/AU_MAF_noPrune.vcf.gz
bcftools index final_SNPs/TI_MAF_noPrune.vcf.gz

bcftools isec final_SNPs/AU_MAF_noPrune.vcf.gz final_SNPs/TI_MAF_noPrune.vcf.gz -p final_SNPs/

grep -v "#" final_SNPs/0000.vcf | wc -l # Sites private to AU
grep -v "#" final_SNPs/0001.vcf | wc -l # Sites private to TI
grep -v "#" final_SNPs/0002.vcf | wc -l # Shared sites
```
Mean Depth per SNP per Individual
```
while read -r line
    do
    echo "Averaging depth for $line..."
    dp=$(bcftools view -s $line ${data}final_SNPs/global_MAF_noPrune.vcf.gz | bcftools query -f '[%DP]\n' | awk '{sum += $1}; END {print sum/NR}')
    echo "${line},${dp}" >> indiv_dp.csv
done < samples.txt
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

bcftools view -S ${data}AU_samples.tsv ${data}global_MAF.vcf.gz | bcftools view -e 'N_PASS(GT=="RR")=15' | bcftools stats | grep "TSTV"
bcftools view  -S ${data}TI_samples.tsv ${data}global_MAF.vcf.gz | bcftools view -e 'N_PASS(GT=="RR")=11' | bcftools stats | grep "TSTV"
bcftools view -e 'N_PASS(GT=="RR")=26' ${data}global_MAF.vcf.gz | bcftools stats | grep "TSTV"
```
Finally, an overall summary of global statistics was reviewed by:
```
bcftools stats -S samples.tsv global_MAF_LDpruned.vcf.gz > total_pruned_stats.vchk
plot-vcfstats --sample-names total_pruned_stats.vchk -p bcftools_plots/
```