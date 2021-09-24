#!/bin/sh
######################################################################################################################
# In this script, I'm trying to see how often SVs and SNPs intersect with gene and/or gene regions.
# First, I took the Smoove VCF and simplified it to:
# CHROM START_POSITION END_POSITION SVTYPE
######################################################################################################################

echo "Creating bed file with SVs..."
zcat kakapo_noW_Mendel.vcf.gz \
	| grep -v "#" \
	| awk '{print $1"\t"$2"\t"$8}' \
	| tr ";" "\t" \
	| awk '{print $1"\t"$2"\t"$5"\t"$3}' | grep -v BND > DEL_DUP_INV.bed

echo "Sorting GFF with bedtools..."
bedtools sort \
	-i GCF_004027225.2_bStrHab1.2.pri_genomic.gff.gz > GCF_004027225.2_bStrHab1.2.pri_genomic.sorted.gff

echo "Sorting SV bed file with bedtools..."
bedtools sort \
	-i DEL_DUP_INV.bed > DEL_DUP_INV.sorted.bed

echo "Running bedtools to find overlapping characters..."
bedtools closest \
	-a GCF_004027225.2_bStrHab1.2.pri_genomic.sorted.gff \
	-b DEL_DUP_INV.sorted.bed -d

echo "Simplifying the output..."

sed 's/collected-by=Tim Raemaekers%2C Daryl Eason%2C Andrew Digby of Department of Conservation%2C New Zealand;collection-date=2013-04-13;country=New Zealand: Anchor Island;gbkey=Src;genome=chromosome;isolate=Jane;lat-lon=45.757406 S 166.504927 E;mol_type=genomic DNA;sex=female;tissue-type=blood in lysis buffer or EtOH//g' trial.bed
