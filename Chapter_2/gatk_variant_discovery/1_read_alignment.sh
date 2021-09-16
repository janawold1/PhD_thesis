#!/bin/bash -e
ref=/data/tara_iti_shortreads/reference/bSteHir1.pri.cur.20190820.fasta.gz #Reference genome for alignment
datadir=/data/tara_iti_shortreads/trimmed_reads/lib2/ #Directory with fastq data
samdir=/data/tara_iti_shortreads/alignments/sam/ #Sam file output
bamdir=/data/tara_iti_shortreads/alignments/bam/ #Bam file output
fq1=_R1.fq.gz #Read 1 suffix
fq2=_R2.fq.gz #Read 2 suffix
platform="Illumina"

#First index the reference genome
#time bwa index $ref

#Now, retrieving read group and instrument information.
for samp in ${datadir}*_R1.fq.gz #Remember to be explicit with file location
do
    base=$(basename ${samp} _R1.fq.gz)
    infoline=$(zcat ${samp} | head -n 1)
    instrument=`echo ${infoline} | cut -d ':' -f1`
    instrumentrun=`echo $infoline | cut -d ':' -f2`
    flowcell=`echo $infoline | cut -d ':' -f3`
    lane=`echo $infoline | cut -d ':' -f4`
    index=`echo $infoline | cut -d ':' -f10`
    name=$(basename ${samp} _lib2_R1.fq.gz)

    #Now to incorporate this information into the alignment
    rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
    rgpl="PL:${platform}"
    rgpu="PU:${flowcell}.${lane}"
    rglb="LB:${base}_library1"
    rgsm="SM:${name}"

    echo "Aligning reads for $base" #Be explicit with file location for read 2 and the sam file output
    time bwa mem -M -R @RG'\t'$rgid'\t'$rgpl'\t'$rgpu'\t'$rglb'\t'$rgsm -t 64 $ref $samp ${datadir}${base}${fq2} > ${samdir}${base}.sam

    echo "Converting sam file to bam file for $base"
    time samtools view -@ 32 -T $ref -b ${samdir}${base}.sam | samtools sort - -@ 32 -o ${bamdir}${base}.sorted.bam
    time samtools index -@ 32 ${bamdir}${base}.sorted.bam
done
