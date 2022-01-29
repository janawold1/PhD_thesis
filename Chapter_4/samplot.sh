#!/bin/sh -e
##########################################################################################################################
# A brief script for running the SAMplot software package on a filtered SV data set. 
##########################################################################################################################
data=/kakapo-data/bwa/

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\n' ${data}smoove/bwa_smoove_MSHQ3_DHFFC0.7_DHBFC1.3.vcf.gz > ${data}smoove/samplot/smoove_samplot.tsv

while read -r line
    do
    echo "Running SAMplot for ${line}"
    chrom=$(echo ${line} | awk '{print $1}')
    start=$(echo ${line} | awk '{print $2}')
    end=$(echo ${line} | awk '{print $3}')
    type=$(echo ${line} | awk '{print $4}')
    time samplot plot \
        -n Bella Bill Blades Richard_Henry Margaret-Maree Rangi Sue \
        -b /kakapo-data/bwa/bwa_female/nodup_bam/batch02/Bella_nodup.bam \
        /kakapo-data/bwa/bwa_male/nodup_bam/batch02/Bill_nodup.bam \
        /kakapo-data/bwa/bwa_male/nodup_bam/batch02/Blades_nodup.bam \
        /kakapo-data/bwa/bwa_male/nodup_bam/batch10/Richard_Henry_nodup.bam \
        /kakapo-data/bwa/bwa_female/nodup_bam/batch06/Margaret-Maree_nodup.bam \
        /kakapo-data/bwa/bwa_male/nodup_bam/batch10/Rangi_nodup.bam \
        /kakapo-data/bwa/bwa_female/nodup_bam/batch09/Sue_nodup.bam \
        -o ${data}smoove/samplot/${chrom}_${start}_${end}_${type}.png \
        -c ${chrom} \
        -s ${start} \
        -e ${end} \
        -t ${type}
done < ${data}smoove/samplot/smoove_samplot.tsv