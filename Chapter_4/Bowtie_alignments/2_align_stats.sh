#########################################################################################
# Script augmented from Joseph Guhlin to estimate alignment stats. The samtools, mosdepth
# and qualimap software packages are installed under the align_stat conda environment.
#########################################################################################
mref=/kakapo-data/References/kakapo_no_Wchromosome.fa
fref=/kakapo-data/References/kakapo_full_ref.fa
mdata=/kakapo-data/bowtie_male/bam/
fdata=/kakapo-data/bowtie_female/bam/
mstat=/kakapo-data/bowtie_male/stat/
fstat=/kakapo-data/bowtie_female/stat/
mout=/kakapo-data/bowtie_male/markdup/
fout=/kakapo-data/bowtie_female/markdup/

for female in ${fdata}*.bam
    do
        base=$(basename ${female} .bam)
        echo "Running calculating stats for ${base}..."
        samtools stats ${female} > ${fstat}${base}.stats &
        mosdepth -n --fasta ${fref} ${fstat}${base} ${female} &
        qualimap bamqc -bam ${female} -nw 10000 -nt 16 -c -outdir ${fstat}${base}.graphmap --java-mem-size=8G
        wait
        echo "Now preparing to mark duplicates for ${base}..."
        samtools sort -@ 32 -n -o ${fout}${base}.nsorted.bam ${female}
        samtools fixmate -@ 32 -r -m -c ${fout}${base}.nsorted.bam ${fout}${base}.fixmate.bam
        samtools sort -@ 32 -o ${fout}${base}.fixmate.sorted.bam ${fout}${base}.fixmate.bam
        samtools markdup -@ 32 ${fout}${base}.fixmate.sorted.bam ${fout}${base}_markdup.bam
done &

for male in ${mdata}*.bam
    do
        mbase=$(basename ${male} .bam)
        echo "Running calculating stats for ${mbase}..."
        samtools stats ${male} > ${mstat}${mbase}.stats &
        mosdepth -n --fasta ${mref} ${mstat}${mbase} ${male} &
        qualimap bamqc -bam ${male} -nw 10000 -nt 16 -c -outdir ${mstat}${mbase}.graphmap --java-mem-size=8G
        wait
        echo "Now preparing to mark duplicates for ${base}..."
        samtools sort -@ 32 -n -o ${mout}${mbase}.nsorted.bam ${male}
        samtools fixmate -@ 32 -r -m -c ${mout}${mbase}.nsorted.bam ${mout}${mbase}.fixmate.bam
        samtools sort -@ 32 -o ${mout}${mbase}.fixmate.sorted.bam ${mout}${mbase}.fixmate.bam
        samtools markdup -@ 32 ${mout}${mbase}.fixmate.sorted.bam ${mout}${mbase}_markdup.bam
done