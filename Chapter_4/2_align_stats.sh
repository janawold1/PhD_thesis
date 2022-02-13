#########################################################################################
# Script augmented from Joseph Guhlin to estimate alignment stats. The samtools, mosdepth
# and qualimap software packages are installed under the align_stat conda environment.
#########################################################################################
data=/data/common_tern/alignments/
ref=/data/reference/bSteHir1.pri.cur.20190820.fasta

for bam in ${data}merged_bam/*.bam
do
    base=$(basename ${bam} .bam)
    echo "Now preparing to mark duplicates for ${base}..."
    samtools sort -@ 8 -n -o ${data}nodup_bam/${base}.nsorted.bam ${bam}
    samtools fixmate -@ 8 -r -m -c ${data}nodup_bam/${base}.nsorted.bam \
        ${data}nodup_bam/${base}.fixmate.bam
    samtools sort -@ 8 -o ${data}nodup_bam/${base}.fixmate.sorted.bam \
        ${data}nodup_bam/${base}.fixmate.bam
    samtools nodup -@ 8 ${data}nodup_bam/${base}.fixmate.sorted.bam \
        ${data}nodup_bam/${base}_nodup.bam &
done
wait
for bam in ${data}nodup_bam/*_nodup.bam
    do
    base=$(basename ${bam} _nodup.bam)
    echo "Running Qualimap for ${base}..."
    qualimap bamqc \
        -bam ${bam} \
        -nw 10000 \
        -nt 16 -c \
        -outdir ${data}nodup_bam_stats/${base}.graphmap \
        --java-mem-size=8G
    echo "Running calculating stats for ${base}..."
    samtools stats ${bam} > ${data}nodup_bam_stats/${base}.stats
    mosdepth -n --fasta ${ref} ${data}nodup_bam_stats/${base} ${bam}
done