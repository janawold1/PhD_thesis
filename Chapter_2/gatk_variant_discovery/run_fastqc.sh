data=/data/tara_iti_shortreads/trimmed_reads/lib2/
out=/data/tara_iti_shortreads/trim_reports/lib2/

for fq in ${data}*_R1.fq.gz
do
	base=$(basename ${fq} _R1.fq.gz)
	echo "Running fastqc for ${base} trimmed reads..."
	fastqc -o ${out} ${fq} &
	fastqc -o ${out} ${data}${base}_R2.fq.gz &
	fastqc -o ${out} ${data}${base}_R1_unpaired.fq.gz &
	fastqc -o ${out} ${data}${base}_R2_unpaired.fq.gz
	wait
done

