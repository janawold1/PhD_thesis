#!/bin/bash -e
###############################################################################################################################
# Here is the process for aligning sequence data from the kāākpō125+ consortium. Reads were mapped to sex specific reference
# genomes. 
###############################################################################################################################
fref=/kakapo-data/References/kakapo_full_ref.fa #Reference genome for female alignment
mref=/kakapo-data/References/kakapop_no_Wchromosome.fa #Reference genome for male alignment
female=/kakapo-data/raw_reads/females/ #Directory with female fastq data
male=/kakapo-data/raw_reads/males/ #Directory with female fastq data
fsamdir=/kakapo-data/bwa/bwa_female/sam/ #Female SAM output
msamdir=/kakapo-data/bwa/bwa_male/sam/ #Male SAM output
fbamdir=/kakapo-data/bwa/bwa_female/bam/ #Female BAM file output
mbamdir=/kakapo-data/bwa/bwa_male/bam/ #Male BAM file output
fq1=_R1.fastq.gz #Read 1 suffix
fq2=_R2.fastq.gz #Read 2 suffix
platform="Illumina"

echo "Indexing the reference genome $ref"
time bwa index $fref
time bwa index $mref

echo "Now beginning alignments...."
for fsamp in ${female}*_R1.fq.gz 
    do
    #Remember to be explicit with file location
    base=$(basename $samp _R1.fq.gz)
    infoline=$(zcat ${samp} | head -n 1)
    instrument=`echo ${infoline} | cut -d ':' -f1`
    instrumentrun=`echo ${infoline} | cut -d ':' -f2`
    flowcell=`echo ${infoline} | cut -d ':' -f3`
    lane=`echo ${infoline} | cut -d ':' -f4`
    index=`echo ${infoline} | cut -d ':' -f10`
    name=$(echo ${base} | sed 's/_L00[1-9]//g')

    #Now to incorporate this information into the alignment
    rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
    rgpl="PL:${platform}"
    rgpu="PU:${flowcell}.${lane}"
    rglb="LB:${base}_${lane}"
    rgsm="SM:${name}"

    #Be explicit with file location for read 2 and the sam file output
    echo "Aligning reads for $base" 
    time bwa mem -M -R @RG'\t'${rgid}'\t'${rgpl}'\t'${rgpu}'\t'${rglb}'\t'${rgsm} -t 64 ${fref} ${fsamp} ${female}${base}${fq2} > ${fsamdir}${base}.sam

    echo "Converting sam file to bam file for $base"
    time samtools view -T $fref -b ${fsamdir}${base}.sam > ${fbamdir}${base}.bam
    time samtools sort -@ 64 -o ${fbamdir}${base}.sorted.bam ${fbamdir}${base}.bam
done &

for msamp in ${male}*_R1.fq.gz
    do
    base=$(basename $msamp _R1.fq.gz)
    infoline=$(zcat ${msamp} | head -n 1)
    instrument=`echo ${infoline} | cut -d ':' -f1`
    instrumentrun=`echo ${infoline} | cut -d ':' -f2`
    flowcell=`echo ${infoline} | cut -d ':' -f3`
    lane=`echo ${infoline} | cut -d ':' -f4`
    index=`echo ${infoline} | cut -d ':' -f10`
    name=$(echo ${base} | sed 's/_L00[1-9]//g')

    #Now to incorporate this information into the alignment
    rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
    rgpl="PL:${platform}"
    rgpu="PU:${flowcell}.${lane}"
    rglb="LB:${base}_${lane}"
    rgsm="SM:${name}"

    #Be explicit with file location for read 2 and the sam file output
    echo "Aligning reads for $base" 
    time bwa mem -M -R @RG'\t'${rgid}'\t'${rgpl}'\t'${rgpu}'\t'${rglb}'\t'${rgsm} -t 64 ${mref} ${msamp} ${male}${base}${fq2} > ${msamdir}${base}.sam

    echo "Converting sam file to bam file for $base"
    time samtools view -T $mref -b ${msamdir}${base}.sam > ${mbamdir}${base}.bam
    time samtools sort -@ 64 -o ${mbamdir}${base}.sorted.bam ${mbamdir}${base}.bam
done
wait