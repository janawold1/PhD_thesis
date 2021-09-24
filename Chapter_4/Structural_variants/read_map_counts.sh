#!/bin/sh
#########################################################################################
# Estimating relative proportion of reads mapped in each individual.
#########################################################################################

for bam in alignments_female/*.sorted.bam
    do
    base=$(basename ${bam} .sorted.bam)
    echo "Calculating read counts for ${base}..."
    paired=$(samtools view -@ 64 -c -f 1 ${bam}) #Count all reads
    primary=$(samtools view -@ 64 -c -F 256 ${bam}) #Count all primary reads
    not_prime=$(samtools view -@ 64 -c -f 256 ${bam}) #Count all non primary reads
    mapped=$(samtools view -@ 64 -c -f 3 -F 256 ${bam}) #Count all primary reads with both ends mapped
    unmap_unmap=$(samtools view -@ 64 -c -f 13 -F 256 ${bam}) #Count all primary reads with both ends not mapped
    unmap_map=$(samtools view -@ 64 -c -f 5 -F 264 ${bam}) #Count all primary reads with the first end not mapped and the second end mapped
    map_unmap=$(samtools view -@ 64 -c -f 9 -F 260 ${bam}) #Count all promary reads with the first end mapped and the second end not mapped
    echo "${base},${paired},${primary},${not_prime},${mapped},${unmap_unmap},${unmap_map},${map_unmap}" >> female_read_maps.txt
done

for bam in alignments_male/*.sorted.bam
    do
    base=$(basename ${bam} .sorted.bam)
    echo "Calculating read counts for ${base}..."
    paired=$(samtools view -@ 64 -c -f 1 ${bam}) #Count all reads
    primary=$(samtools view -@ 64 -c -F 256 ${bam}) #Count all primary reads
    not_prime=$(samtools view -@ 64 -c -f 256 ${bam}) #Count all non primary reads
    mapped=$(samtools view -@ 64 -c -f 3 -F 256 ${bam}) #Count all primary reads with both ends mapped
    unmap_unmap=$(samtools view -@ 64 -c -f 13 -F 256 ${bam}) #Count all primary reads with both ends not mapped
    unmap_map=$(samtools view -@ 64 -c -f 5 -F 264 ${bam}) #Count all primary reads with the first end not mapped and the second end mapped
    map_unmap=$(samtools view -@ 64 -c -f 9 -F 260 ${bam}) #Count all promary reads with the first end mapped and the second end not mapped
    echo "${base},${paired},${primary},${not_prime},${mapped},${unmap_unmap},${unmap_map},${map_unmap}" >> male_read_maps.txt
done

#########################################################################################
# Estimating relative proportion of reads mapped to each chromosome for each individual.
#########################################################################################

for bam in alignments_female/*.sorted.bam
    do
    base=$(basename ${bam} .sorted.bam)
    while read -r line
        do
        echo "Calculating chromosomal read counts for ${base} at chromosome ${line}..."
        paired=$(samtools view -@ 64 -c -f 1 ${bam} ${line}) #Count all reads
        primary=$(samtools view -@ 64 -c -F 256 ${bam} ${line}) #Count all primary reads
        not_prime=$(samtools view -@ 64 -c -f 256 ${bam} ${line}) #Count all non primary reads
        mapped=$(samtools view -@ 64 -c -f 3 -F 256 ${bam} ${line}) #Count all primary reads with both ends mapped
        unmap_map=$(samtools view -@ 64 -c -f 5 -F 264 ${bam} ${line}) #Count all primary reads with the first end not mapped and the second end mapped
        map_unmap=$(samtools view -@ 64 -c -f 9 -F 260 ${bam} ${line}) #Count all promary reads with the first end mapped and the second end not mapped
        echo "${base},${line},${primary},${not_prime},${mapped},${unmap_map},${map_unmap}" >> female_chr_read_maps.txt
    done < chromosomes.txt
done

for bam in alignments_male/*.sorted.bam
    do
    base=$(basename ${bam} .sorted.bam)
    while read -r line
        do
        echo "Calculating chromosomal read counts for ${base} at chromosome ${line}..."
        paired=$(samtools view -@ 64 -c -f 1 ${bam} ${line}) #Count all reads
        primary=$(samtools view -@ 64 -c -F 256 ${bam} ${line}) #Count all primary reads
        not_prime=$(samtools view -@ 64 -c -f 256 ${bam} ${line}) #Count all non primary reads
        mapped=$(samtools view -@ 64 -c -f 3 -F 256 ${bam} ${line}) #Count all primary reads with both ends mapped
        unmap_map=$(samtools view -@ 64 -c -f 5 -F 264 ${bam} ${line}) #Count all primary reads with the first end not mapped and the second end mapped
        map_unmap=$(samtools view -@ 64 -c -f 9 -F 260 ${bam} ${line}) #Count all promary reads with the first end mapped and the second end not mapped
        echo "${base},${line},${primary},${not_prime},${mapped},${unmap_map},${map_unmap}" >> male_chr_read_maps.txt
    done < chromosomes.txt
done

#########################################################################################
# Estimating average mapping score per bp in 1kb windows.
#########################################################################################

for bam in alignments_female/*.sorted.bam
    do
    base=$(basename ${bam} .sorted.bam)
    while read -r line
        do
        region=$(echo ${line} | awk '{print $1":"$2"-"$3}')
        chr=$(echo ${line} | awk '{print $1}')
        begin=$(echo ${line} | awk '{print $2}')
        end=$(echo ${line} | awk '{print $3}')
        echo "Calculating average bp mapping quality for ${base} at ${region}..."
        qual=$(samtools view ${bam} ${region} | awk '{sum+=$5} END { print sum/1000}')
        echo "${base},${chr},${begin},${end},${qual}" >> read_qual/female_chr_read_qual_1kb.csv
    done < read_qual/chrom_1kb_windows.bed
done

for bam in alignments_male/*.sorted.bam
    do
    base=$(basename ${bam} .sorted.bam)
    while read -r line
        do
        region=$(echo ${line} | awk '{print $1":"$2"-"$3}')
        chr=$(echo ${line} | awk '{print $1}')
        begin=$(echo ${line} | awk '{print $2}')
        end=$(echo ${line} | awk '{print $3}')
        echo "Calculating average bp mapping quality for ${base} at ${region}..."
        qual=$(samtools view ${bam} ${region} | awk '{sum+=$5} END { print sum/1000}')
        echo "${base},${chr},${begin},${end},${qual}" >> read_qual/male_chr_read_qual_1kb.csv
    done < read_qual/chrom_1kb_windows.bed
done