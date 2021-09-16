#########################################################################################
# Script augmented from Joseph Guhlin to estimate alignment stats. The samtools, mosdepth
# and qualimap software packages are installed under the align_stat conda environment.
#########################################################################################
data=/data/tara_iti_shortreads/alignments/bam/
merged=/data/tara_iti_shortreads/alignments/merged_bam/
stat=/data/tara_iti_shortreads/alignments/stats/
ref=/data/tara_iti_shortreads/reference/bSteHir1.pri.cur.20190820.fasta

for bam in ${merged}*{AU30,AU33,SND04}.bam
do
	base=$(basename ${bam} .bam)
	echo "Running calculating stats for ${base}..."
	samtools stats ${bam} > ${stat}${base}.stats
	mosdepth -n --fasta ${ref} ${stat}${base} ${bam}
	echo "Now preparing to mark duplicates for ${base}..."
	samtools sort -@ 32 -n -o ${merged}${base}.nsorted.bam ${bam}
	samtools fixmate -@ 32 -r -m -c ${merged}${base}.nsorted.bam \
		${merged}${base}.fixmate.bam
	samtools sort -@ 32 -o ${merged}${base}.fixmate.sorted.bam \
		${merged}${base}.fixmate.bam
	samtools markdup -@ 32 ${merged}${base}.fixmate.sorted.bam \
		${merged}${base}_markdup.bam &
done &

for align in ${data}*.sorted.bam
do
base=$(basename ${align} .sorted.bam)
echo "Running Qualimap for ${base}..."
qualimap bamqc \
	-bam ${align} \
	-nw 10000 \
	-nt 16 -c \
	-outdir /data/tara_iti_shortreads/qualimap/${base}.graphmap \
	--java-mem-size=8G &
done
wait
