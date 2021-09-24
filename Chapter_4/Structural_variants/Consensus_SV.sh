#!/bin/bash
########################################################################################################
# Here I have use the BioConductor package Structural Variants and followed the example provided on the 
# GRIDSS github page to create a filtered BED file for each batch. To see the overlap with MANTA calls I
# did the following protocol.
########################################################################################################

bedtools merge -i GRIDSS_females_B01.simple.bed -i GRIDSS_females_B02.simple.bed -i GRIDSS_females_B03.simple.bed \
    -i GRIDSS_females_B04.simple.bed -i GRIDSS_females_B05.simple.bed -i GRIDSS_females_B06.simple.bed \
    -i GRIDSS_females_B07.simple.bed -i GRIDSS_males_B01.simple.bed -i GRIDSS_males_B02.simple.bed -i GRIDSS_males_B03.simple.bed \
    -i GRIDSS_males_B04.simple.bed -i GRIDSS_males_B05.simple.bed -i GRIDSS_males_B06.simple.bed -i GRIDSS_males_B07.simple.bed \
    -i GRIDSS_males_B08.simple.bed -i GRIDSS_males_B09.simple.bed > GRIDSS_merged.simple.bed

# Created a consensus call set between GRIDSS and the trio tested MANTA SV's. The output of this file is
# in MANTA format.
bcftools view -O z -o MANTA_kakapo_shared_SV.vcf.gz -R GRIDSS_merged.simple.bed manta/all_MANTA_mendel.vcf.gz

# Checked the counts of remaining SV's with:
zcat MANTA_kakapo_shared_SV.vcf.gz | grep -v "^#" | awk '{print $1, $2, $3}' | tr ":" " " | awk '{print $1, $2, $3}' \
    > SV_consensus_calls_from_manta.tsv

# Then counted SV's per chromosome with:
cat SV_consensus_calls_from_manta.tsv | awk '{print $1, $3}' | sort | uniq -c

# And total consensus SV's with:
cat SV_consensus_calls_from_manta.tsv | awk '{print $3}' | sort | uniq -c

# This indicated that there were 160 SV's shared between GRIDSS and MANTA. However, it's important to note
# that these are locations with identified 'events'. That is to say that there are 160 locations where 
# each programme noted a SV, the called type of SV may not be the same. So, the break down of the events as
# called in MANTA are as follows:
# 3 MantaBND
# 100 MantaDEL
# 21 MantaDUP
# 17 MantaINS
# 19 MantaINV

# To see whether the called sites are the same type in GRIDSS, I also filtered the GRIDSS consensus calls with this\
# merged bed file as per:
bcftools view -O z -o GRIDSS_kakapo_shared_SV.vcf.gz -R GRIDSS_merged.simple.bed GRIDSS_annotated_kakapo.vcf.gz
 