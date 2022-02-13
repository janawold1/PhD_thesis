#########################################################################################
# Script augmented from Joseph Guhlin to estimate alignment stats. The samtools, mosdepth
# and qualimap software packages are installed under the align_stat conda environment.
#########################################################################################
data=/data/common_tern/alignments/
ref=/data/reference/bSteHir1.pri.cur.20190820.fasta

for bam in ${data}merged_bam/*.bam
do
    base=$(basename ${bam} .bam)
    echo "Running calculating stats for ${base}..."
    samtools stats ${bam} > ${data}stats/${base}.stats
    mosdepth -n --fasta ${ref} ${data}stats/${base} ${bam}
    echo "Now preparing to mark duplicates for ${base}..."
    samtools sort -@ 8 -n -o ${data}markdup_bam/${base}.nsorted.bam ${bam}
    samtools fixmate -@ 8 -r -m -c ${data}markdup_bam/${base}.nsorted.bam \
        ${data}markdup_bam/${base}.fixmate.bam
    samtools sort -@ 8 -o ${data}markdup_bam/${base}.fixmate.sorted.bam \
        ${data}markdup_bam/${base}.fixmate.bam
    samtools markdup -@ 8 ${data}markdup_bam/${base}.fixmate.sorted.bam \
        ${data}markdup_bam/${base}_markdup.bam &
done
wait
for align in ${data}*.bam
    do
    base=$(basename ${align} .sorted.bam)
    echo "Running Qualimap for ${base}..."
    qualimap bamqc \
        -bam ${align} \
        -nw 10000 \
        -nt 16 -c \
        -outdir ${data}stats/${base}.graphmap \
        --java-mem-size=8G
done