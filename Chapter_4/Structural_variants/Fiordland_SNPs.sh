#!/bin/sh
########################################################################################################
# Here I have used Lara's fiordland SNP set found at:
# /nesi/nobackup/uoo02695/Kakapo/05_GWAS_herit_partitioning/plink_final/hwe/no_maf/fiordland_hw_snps.tsv
# to identify variable SNPs through RH's lineage and how they correlate with SVs. 
########################################################################################################

# First convert file format to chr pos format for bcftools 
awk'{print $2}' fiordland_hw_snps.tsv | sed 's/_/\t/g' | awk '{print $1"\t"$2}' > fiordland_hw_snps_pos.tsv

# Then made sure I was only looking at SNPs that were variable in RH...
bcftools view -R Fiordland_hw_snps_pos.tsv \
    -s Richard_Henry Trained.bcf \
    | grep -v "#" \
    | grep -v "0/0" > \
    RH_fiordland_variable.tsv
bcftools view -O z -o RH_fiordland_snps_variable.vcf.gz \
    -R RH_fiordland_variable.tsv \
    -s Richard_Henry \
    Trained.bcf

# Now to see how many of these SNPs unique to RH were inherited through the F1 and F2 crosses. First 
# created vcf's for RH offspring with sites variable in RH

while read -r line
    do
    echo "Working on counting SNPs for ${line}..."
    echo ${line} >> RH_SNP_counts.tsv
    bcftools view -R RH_variable_fiordland.tsv \
        -s ${line} \
        Trained.bcf \
        | grep -v "#" \
        | grep -v "0/0" \
        | wc -l >> RH_SNP_counts.tsv
done < RH_lineage.txt

# Here I had a look at the distance between all unfiltered variable SVs in RH and SNPs unique to RH lineage.

zcat RH_fiordland_snps_variable.vcf.gz | grep -v "#" | 