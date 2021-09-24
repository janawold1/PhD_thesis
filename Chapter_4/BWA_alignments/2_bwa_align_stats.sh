#!/bin/bash -e
#########################################################################################
# Script augmented from Joseph Guhlin to estimate alignment stats. The samtools, mosdepth
# and qualimap software packages are installed under the align_stat conda environment.
#########################################################################################
mref=/kakapo-data/References/kakapo_no_Wchromosome.fa
fref=/kakapo-data/References/kakapo_full_ref.fa
mdata=/kakapo-data/bwa_male/
fdata=/kakapo-data/bwa_female/
mstat=/kakapo-data/bwa_male/stat/
fstat=/kakapo-data/bwa_female/stat/
mout=/kakapo-data/bwa_male/markdup/
fout=/kakapo-data/bwa_female/markdup/

for i in {01..11}
        do
        for female in ${fdata}batch${i}/*.bam
                do
                base=$(basename ${female} .bam)
                echo "Running calculating stats for ${base}..."
                samtools stats -@64 ${female} > ${fstat}${base}.stats &
                mosdepth -n --fasta ${fref} ${fstat}${base} ${female} &
                qualimap bamqc -bam ${female} -nw 10000 -nt 32 -c -outdir ${fstat}${base}.graphmap --java-mem-size=8G
                wait
                echo "Now preparing to mark duplicates for ${base}..."
                samtools sort -@ 64 -n -o ${fout}${base}.nsorted.bam ${female}
                samtools fixmate -@ 64 -r -m -c ${fout}${base}.nsorted.bam ${fout}${base}.fixmate.bam
                samtools sort -@ 64 -o ${fout}${base}.fixmate.sorted.bam ${fout}${base}.fixmate.bam
                samtools markdup -@ 64 ${fout}${base}.fixmate.sorted.bam ${fout}${base}_markdup.bam
        done &
done &

for j in {01..14}
        do
        for male in ${mdata}batch${j}/*.bam
                do
                mbase=$(basename ${male} .bam)
                echo "Running calculating stats for ${mbase}..."
                samtools stats -@ 64 ${male} > ${mstat}${mbase}.stats &
                mosdepth -n --fasta ${mref} ${mstat}${mbase} ${male} &
                qualimap bamqc -bam ${male} -nw 10000 0-nt 32 -c -outdir ${mstat}${mbase}.graphmap --java-mem-size=8G
                wait
                echo "Now preparing to mark duplicates for ${base}..."
                samtools sort -@ 64 -n -o ${mout}${mbase}.nsorted.bam ${male}
                samtools fixmate -@ 64 -r -m -c ${mout}${mbase}.nsorted.bam ${mout}${mbase}.fixmate.bam
                samtools sort -@ 64 -o ${mout}${mbase}.fixmate.sorted.bam ${mout}${mbase}.fixmate.bam
                samtools markdup -@ 64 ${mout}${mbase}.fixmate.sorted.bam ${mout}${mbase}_markdup.bam
        done &
done
