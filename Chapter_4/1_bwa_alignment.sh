#!/bin/bash -e
###############################################################################################################################
# Still need to fix.
###############################################################################################################################
ref=/kakapo-data/References/kakapo_full_ref.fa #Reference genome for alignment
datadir=/kakapo-data/subsampled_fasta/fasta20X/female/ #Directory with fastq data
samdir=/kakapo-data/subsampled_alignments/fasta20X/female/sam/ #Sam file output
bamdir=/kakapo-data/subsampled_alignments/fasta20X/female/bam/ #Bam file output
fq1=_R1.fastq.gz #Read 1 suffix
fq2=_R2.fastq.gz #Read 2 suffix
platform="Illumina"

#First index the reference genome
#time bwa index $ref

#Now, retrieving read group and instrument information.
for samp in ${datadir}*_R1.fq.gz 
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
    time bwa mem -M -R @RG'\t'${rgid}'\t'${rgpl}'\t'${rgpu}'\t'${rglb}'\t'${rgsm} -t 64 ${ref} ${samp} ${datadir}${base}${fq2} > ${samdir}${base}.sam

    echo "Converting sam file to bam file for $base"
    time samtools view -T $ref -b ${samdir}${base}.sam > ${bamdir}${base}.bam
    time samtools sort -@ 64 -o ${bamdir}${base}.sorted.bam ${bamdir}${base}.bam
done