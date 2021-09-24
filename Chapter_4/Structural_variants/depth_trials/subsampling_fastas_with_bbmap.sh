####################################################################
# Subsampling processed fasta files, outputs are based on the target
# sequencing depth of 30X coverage. 
####################################################################

workdir=/kakapo-data/raw_reads
outdir=/kakapo-data/subsampled_fasta

for fq1 in ${workdir}/*_R1.fastq.gz
do
    fq2=$(basename ${fq1} _R1.fastq.gz)
    echo "Subsampling ${fq2} to 20X coverage..."
    reformat.sh in1=${fq1} \
        in2=${workdir}/${fq2}_R2.fastq.gz \
        out1=${outdir}/fasta20X/${fq2}_R1.fq.gz \
        out2=${outdir}/fasta20X/${fq2}_R2.fq.gz \
        samplerate=0.67 \
        -Xmx30g
    echo "Subsampling ${fq2} to 15X coverage..."
    reformat.sh in1=${fq1} \
        in2=${workdir}/${fq2}_R2.fastq.gz \
        out1=${outdir}/fasta15X/${fq2}_R1.fq.gz \
        out2=${outdir}/fasta15X/${fq2}_R2.fq.gz \
        samplerate=0.5 \
        -Xmx30g
    echo "Subsampling ${fq2} to 10X coverage..."
    reformat.sh in1=${fq1} \
        in2=${workdir}/${fq2}_R2.fastq.gz \
        out1=${outdir}/fasta10X/${fq2}_R1.fq.gz \
        out2=${outdir}/fasta10X/${fq2}_R2.fq.gz \
        samplerate=0.33 \
        -Xmx30g
    echo "Subsampling ${fq2} to 5X coverage..."
    reformat.sh in1=${fq1} \
        in2=${workdir}/${fq2}_R2.fastq.gz \
        out1=${outdir}/fasta5X/${fq2}_R1.fq.gz \
        out2=${outdir}/fasta5X/${fq2}_R2.fq.gz \
        samplerate=0.17 \
        -Xmx30g
done