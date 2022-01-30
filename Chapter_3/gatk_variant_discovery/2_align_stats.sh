data=/data/tara_iti_shortreads/alignments/bam/
out=/data/tara_iti_shortreads/alignments/stats/
ref=/data/tara_iti_shortreads/reference/bSteHir1.pri.cur.20190820.fasta

for bam in ${data}*.sorted.bam
do
	base=$(basename ${bam} .sorted.bam)
	echo "Running calculating stats for ${base}..."
	samtools stats ${bam} > ${out}${base}.stats
	mosdepth -n --fasta ${ref} ${out}${base} ${bam}
	echo "Now preparing to mark duplicates for ${base}..."
	samtools sort -@ 32 -n -o ${data}${base}.nsorted.bam ${bam}
	samtools fixmate -@ 32 -r -m -c ${data}${base}.nsorted.bam \
		${data}${base}.fixmate.bam
	samtools sort -@ 32 -o ${data}${base}.fixmate.sorted.bam \
		${data}${base}.fixmate.bam
	samtools markdup -@ 32 ${data}${base}.fixmate.sorted.bam \
		${data}${base}_markdup.bam
done
