#!/bin/sh
#############################################################################################################
#A brief script for getting the size ranges for each SV into a format to plot in R.
#############################################################################################################

cat all_Mendel.vcf | grep -v "^#" | awk '{print $1, $2, $3, $8}' | tr ":" "\t" | awk '{print $1, $2, $3, $100}'| tr ";" "\t" | \
awk '{print$1,$2,$3,$6}' | sed 's/SVLEN=//g' | grep -v MantaBND | grep -v MantaDUP | \
grep -v CIPOS > Mendel_INS_INV_DEL_sizes.tsv

zcat all_noMendel.vcf.gz | grep -v "^#" | awk '{print $1, $2, $3, $8}' | tr ":" "\t" | awk '{print $1, $2, $3, $100}'| \
tr ";" "\t" | awk '{print$1,$2,$3,$6}' | sed 's/SVLEN=//g' | grep -v MantaBND | grep -v MantaDUP | \
grep -v CIPOS > noMendel_INS_INV_DEL_sizes.tsv

zcat kakapo_SVs.vcf.gz | grep -v "^#" | awk '{print $1, $2, $3, $8}' | tr ":" "\t" | awk '{print $1, $2, $3, $100}'| \
tr ";" "\t" | awk '{print$1,$2,$3,$6}' | sed 's/SVLEN=//g' | grep -v MantaBND | grep -v MantaDUP | \
grep -v CIPOS > unfiltered_SV_sizes_INS_INV_DEL.tsv

cat all_Mendel.vcf | grep -v "^#" | awk '{print $1, $2, $3, $8}' | tr ":" "\t" | awk '{print $1, $2, $3, $100}'| tr ";" "\t" | \
awk '{print$1,$2,$3,$6}' | sed 's/SVLEN=//g' > Mendel_INS_INV_DEL_sizes.tsv