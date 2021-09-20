#########################################################################################
# Script augmented from Joseph Guhlin to estimate alignment stats. The samtools, mosdepth
# and qualimap software packages are installed under the align_stat conda environment.
#########################################################################################
mref=/kakapo-data/References/kakapo_no_Wchromosome.fa
fref=/kakapo-data/References/kakapo_full_ref.fa
mdata=/kakapo-data/bowtie_male/
fdata=/kakapo-data/bowtie_femadata
markdup/
for female in ${fdata}bam/*.bam
    do
        base=$(basename ${female} .bam)
        echo "Running calculating stats for ${base}..."
        samtools stats ${female} > ${fdata}stat/${base}.stats &
        mosdepth -n --fasta ${fref} ${fdata}stat/${base} ${female} &
        qualimap bamqc -bam ${female} -nw 10000 -nt 16 -c -outdir ${fdata}stat/${base}.graphmap --java-mem-size=8G
        wait
        echo "Now preparing to mark duplicates for ${base}..."
        samtools sort -@ 32 -n -o ${fdata}markdup/${base}.nsorted.bam ${female}
        samtools fixmate -@ 32 -r -m -c ${fdata}markdup/${base}.nsorted.bam ${fdata}markdup/${base}.fixmate.bam
        samtools sort -@ 32 -o ${fdata}markdup/${base}.fixmate.sorted.bam ${fdata}markdup/${base}.fixmate.bam
        samtools markdup -@ 32 ${fdata}markdup/${base}.fixmate.sorted.bam ${fdata}markdup/${base}_markdup.bam
done &

for male in ${mdata}bam/*.bam
    do
        mbase=$(basename ${male} .bam)
        echo "Running calculating stats for ${mbase}..."
        samtools stats ${male} > ${mdata}stat/${mbase}.stats &
        mosdepth -n --fasta ${mref} ${mdata}stat/${mbase} ${male} &
        qualimap bamqc -bam ${male} -nw 10000 -nt 16 -c -outdir ${mdata}stat/${mbase}.graphmap --java-mem-size=8G
        wait
        echo "Now preparing to mark duplicates for ${base}..."
        samtools sort -@ 32 -n -o ${mdata}markdup/${mbase}.nsorted.bam ${male}
        samtools fixmate -@ 32 -r -m -c ${mdata}markdup/${mbase}.nsorted.bam ${mdata}markdup/${mbase}.fixmate.bam
        samtools sort -@ 32 -o ${mdata}markdup/${mbase}.fixmate.sorted.bam ${mdata}markdup/${mbase}.fixmate.bam
        samtools markdup -@ 32 ${mdata}markdup/${mbase}.fixmate.sorted.bam ${mdata}markdup/${mbase}_markdup.bam
done
