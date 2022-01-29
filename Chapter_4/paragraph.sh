#!/bin/sh -e
#####################################################################################
#Script for running the MANTA software package. Intended for structural
#variant calling. MANTA is aprogram for the detection of structural variants
#from paired-end sequence data. Because running all individuals together is taxing on 
#MANTA, this program was run in 7 batches for female kakapo. 
#####################################################################################
data=/kakapo-data/bwa/
work=/kakapo-data/bwa/manta/archived/paragraph/
vcf=/kakapo-data/bwa/manta/archived/paragraph/joint_total.vcf.gz

echo "Creating manifests for samples..."
for fbam in ${data}bwa_female/nodup_bam/batch*/*_nodup.bam
    do
    female=$(basename ${fbam} _nodup.bam)
    path=$(ls ${fbam})
    depth=$(samtools depth ${fbam} | awk '{sum+=$3}; END {print sum/NR}')
    printf "$id\t$path\t$depth\t125\n"
    printf "$base\t$path\t$depth\t125\n" >> ${work}female_manifest.tsv
done &
for mbam in ${data}bwa_male/nodup_bam/batch*/*_nodup.bam
    do
    male=$(basename ${mbam} _nodup.bam)
    path=$(ls ${mbam})
    depth=$(samtools depth ${mbam} | awk '{sum+=$3}; END {print sum/NR}')
    printf "$id\t$path\t$depth\t125\n"
    printf "$base\t$path\t$depth\t125\n" >> ${work}male_manifest.tsv
done
wait

