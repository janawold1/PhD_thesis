#!/bin/bash -e
############################################################################################################################################
# Script used for the alignment of processed reads
############################################################################################################################################
fref=/kakapo-data/References/kakapo_full_ref #Female reference genome for alignment
mref=/kakapo-data/References/kakapo_no_Wchromosome #Male reference genome for alignment
fdata=/kakapo-data/raw_reads/females/ #Directory with female fastq data
mdata=/kakapo-data/raw_reads/males/ #Directory with male fastq data

mlane1=/kakapo-data/raw_reads/males/lanes1-6/
mlane8=/kakapo-data/raw_reads/males/lanes3-8/
fsam=/kakapo-data/bowtie_female/sam/ #Female sam file output
msam=/kakapo-data/bowtie_male/sam/ #Male sam file output
fbam=/kakapo-data/bowtie_female/bam/ #Female bam file output
mbam=/kakapo-data/bowtie_male/bam/ #Male bam file output
fq1=_R1.fastq.gz #Read 1 suffix
fq2=_R2.fastq.gz #Read 2 suffix
platform="Illumina"

#First index the reference genome
#bowtie2-build /kakapo-data/References/kakapo_no_W_chromosome.fa /kakapo-data/References/kakapo_no_Wchromosome
#bowtie2-build /kakapo-data/References/kakapo_full_ref.fa /kakapo-data/References/kakapo_full_ref

#Now, retrieving read group and instrument information.
for female in ${data}females/
do
    for fnolane in ${female}no_lane/*_R1.fastq.gz
    do
        base=$(basename $fnolane _R1.fastq.gz)
        infoline=$(zcat ${fnolane} | head -n 1)
        instrument=`echo ${infoline} | cut -d ':' -f1`
        instrumentrun=`echo ${infoline} | cut -d ':' -f2`
        flowcell=`echo ${infoline} | cut -d ':' -f3`
        lane=`echo ${infoline} | cut -d ':' -f4`
        index=`echo ${infoline} | cut -d ':' -f10`

        #Now to incorporate this information into the alignment
        rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
        rgpl="PL:${platform}"
        rgpu="PU:${flowcell}.${lane}"
        rglb="LB:${base}_${lane}"
        rgsm="SM:${base}"

        echo "Aligning reads for $base" #Be explicit with file location for read 2 and the sam file output
        time bowtie2 -x ${fref} \
            -p 24 \
            --rg-id $rgid --rg $rgpl --rg $rgpu --rg $rglb --rg $rgsm \
            -1 ${fnolane} -2 ${fdatadir}${base}${fq2} \
            -S ${fsamdir}no_lane/${base}.sam

        echo "Converting sam file to bam file for $base"
        time samtools view -T ${fref}.fa -b ${fsamdir}no_lane/${base}.sam | samtools sort - -@ 64 -o ${fbamdir}no_lane/${base}.bam
        samtools index -@ 64 ${fbamdir}no_lane/${base}.bam
    done &
    for lane1 in ${female}lanes1-6/*_L001_R1.fastq.gz
    do
        base=$(basename ${lane1} _R1.fastq.gz)
        infoline=$(zcat ${lane1} | head -n 1)
        instrument=`echo ${infoline} | cut -d ':' -f1`
        instrumentrun=`echo ${infoline} | cut -d ':' -f2`
        flowcell=`echo ${infoline} | cut -d ':' -f3`
        lane=`echo ${infoline} | cut -d ':' -f4`
        index=`echo ${infoline} | cut -d ':' -f10`
        name=$(echo ${base} | sed "s/_L00[1-6]//g")

        #Now to incorporate this information into the alignment
        rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
        rgpl="PL:${platform}"
        rgpu="PU:${flowcell}.${lane}"
        rglb="LB:${base}_${lane}"
        rgsm="SM:${name}"

        echo "Aligning reads for $base" #Be explicit with file location for read 2 and the sam file output
        time bowtie2 -x ${fref} \
            -p 24 \
            --rg-id $rgid --rg $rgpl --rg $rgpu --rg $rglb --rg $rgsm \
            -1 ${lane1} -2 ${fdatadir}lanes_1-6/${base}${fq2} \
            -S ${fsamdir}lanes1-6/${base}.sam

        echo "Converting sam file to bam file for $base"
        time samtools view -T ${fref}.fa -b ${fsamdir}lanes_1-6/${base}.sam | samtools sort - -@ 64 -o ${fbamdir}lanes_1-6/${base}.bam
        samtools index -@ 64 ${fbamdir}lanes_1-6/${base}.bam
    done &
    for lane8 in ${female}lanes3-8/*_L001_R1.fastq.gz
    do
        base=$(basename ${lane8} _R1.fastq.gz)
        infoline=$(zcat ${lane8} | head -n 1)
        instrument=`echo ${infoline} | cut -d ':' -f1`
        instrumentrun=`echo ${infoline} | cut -d ':' -f2`
        flowcell=`echo ${infoline} | cut -d ':' -f3`
        lane=`echo ${infoline} | cut -d ':' -f4`
        index=`echo ${infoline} | cut -d ':' -f10`
        name=$(echo ${base} | sed "s/_L00[1-6]//g")

        #Now to incorporate this information into the alignment
        rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
        rgpl="PL:${platform}"
        rgpu="PU:${flowcell}.${lane}"
        rglb="LB:${base}_${lane}"
        rgsm="SM:${name}"

        echo "Aligning reads for $base" #Be explicit with file location for read 2 and the sam file output
        time bowtie2 -x ${fref} \
            -p 24 \
            --rg-id $rgid --rg $rgpl --rg $rgpu --rg $rglb --rg $rgsm \
            -1 ${lane8} -2 ${fdatadir}lanes_3-8/${base}${fq2} \
            -S ${fsamdir}lanes3-8/${base}.sam

        echo "Converting sam file to bam file for $base"
        time samtools view -T ${fref}.fa -b ${fsamdir}lanes_3-8/${base}.sam | samtools sort - -@ 64 -o ${fbamdir}lanes_3-8/${base}.bam
        samtools index -@ 64 ${fbamdir}lanes_3-8/${base}.bam
    done
 wait
done &

for male in ${data}males/
do
    for mnolane in ${male}no_lane/*_R1.fastq.gz
    do
        base=$(basename $mnolane _R1.fastq.gz)
        infoline=$(zcat ${mnolane} | head -n 1)
        instrument=`echo ${infoline} | cut -d ':' -f1`
        instrumentrun=`echo ${infoline} | cut -d ':' -f2`
        flowcell=`echo ${infoline} | cut -d ':' -f3`
        lane=`echo ${infoline} | cut -d ':' -f4`
        index=`echo ${infoline} | cut -d ':' -f10`

        #Now to incorporate this information into the alignment
        rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
        rgpl="PL:${platform}"
        rgpu="PU:${flowcell}.${lane}"
        rglb="LB:${base}_${lane}"
        rgsm="SM:${base}"

        echo "Aligning reads for $base" #Be explicit with file location for read 2 and the sam file output
        time bowtie2 -x ${mref} \
            -p 24 \
            --rg-id $rgid --rg $rgpl --rg $rgpu --rg $rglb --rg $rgsm \
            -1 ${mnolane} -2 ${mdatadir}${base}${fq2} \
            -S ${msamdir}no_lane/${base}.sam

        echo "Converting sam file to bam file for $base"
        time samtools view -T ${fref}.fa -b ${fsamdir}no_lane/${base}.sam | samtools sort - -@ 64 -o ${fbamdir}no_lane/${base}.bam
        samtools index -@ 64 ${fbamdir}no_lane/${base}.bam
    done &
    for lane1 in ${male}lanes1-6/*_L001_R1.fastq.gz
    do
        base=$(basename ${lane1} _R1.fastq.gz)
        infoline=$(zcat ${lane1} | head -n 1)
        instrument=`echo ${infoline} | cut -d ':' -f1`
        instrumentrun=`echo ${infoline} | cut -d ':' -f2`
        flowcell=`echo ${infoline} | cut -d ':' -f3`
        lane=`echo ${infoline} | cut -d ':' -f4`
        index=`echo ${infoline} | cut -d ':' -f10`
        name=$(echo ${base} | sed "s/_L00[1-6]//g")

        #Now to incorporate this information into the alignment
        rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
        rgpl="PL:${platform}"
        rgpu="PU:${flowcell}.${lane}"
        rglb="LB:${base}_${lane}"
        rgsm="SM:${name}"

        echo "Aligning reads for $base" #Be explicit with file location for read 2 and the sam file output
        time bowtie2 -x ${fref} \
            -p 24 \
            --rg-id $rgid --rg $rgpl --rg $rgpu --rg $rglb --rg $rgsm \
            -1 ${lane1} -2 ${fdatadir}lanes_1-6/${base}${fq2} \
            -S ${fsamdir}lanes1-6/${base}.sam

        echo "Converting sam file to bam file for $base"
        time samtools view -T ${fref}.fa -b ${fsamdir}lanes_1-6/${base}.sam | samtools sort - -@ 64 -o ${fbamdir}lanes_1-6/${base}.bam
        samtools index -@ 64 ${fbamdir}lanes_1-6/${base}.bam
    done &
    for lane8 in ${male}lanes3-8/*_L001_R1.fastq.gz
    do
        base=$(basename ${lane8} _R1.fastq.gz)
        infoline=$(zcat ${lane8} | head -n 1)
        instrument=`echo ${infoline} | cut -d ':' -f1`
        instrumentrun=`echo ${infoline} | cut -d ':' -f2`
        flowcell=`echo ${infoline} | cut -d ':' -f3`
        lane=`echo ${infoline} | cut -d ':' -f4`
        index=`echo ${infoline} | cut -d ':' -f10`
        name=$(echo ${base} | sed "s/_L00[1-6]//g")

        #Now to incorporate this information into the alignment
        rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
        rgpl="PL:${platform}"
        rgpu="PU:${flowcell}.${lane}"
        rglb="LB:${base}_${lane}"
        rgsm="SM:${name}"

        echo "Aligning reads for $base" #Be explicit with file location for read 2 and the sam file output
        time bowtie2 -x ${fref} \
            -p 24 \
            --rg-id $rgid --rg $rgpl --rg $rgpu --rg $rglb --rg $rgsm \
            -1 ${lane8} -2 ${fdatadir}lanes_3-8/${base}${fq2} \
            -S ${fsamdir}lanes3-8/${base}.sam

        echo "Converting sam file to bam file for $base"
        time samtools view -T ${fref}.fa -b ${fsamdir}lanes_3-8/${base}.sam | samtools sort - -@ 64 -o ${fbamdir}lanes_3-8/${base}.bam
        samtools index -@ 64 ${fbamdir}lanes_3-8/${base}.bam
    done
 wait
done