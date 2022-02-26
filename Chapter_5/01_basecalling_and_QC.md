# ONT basecalling, read quality and trimming
Raw MinION fast5 files were basecalled using guppy v6.0.1 under the super 
```
config=dna_r9.4.1_450bps_sup.cfg
guppy_basecaller -i ${data} \
    -s ${output} \
    --config ${config} \
    --verbose_logs \
    --compress_fastq \
    --device cuda:all:100%\
    --recursive \
    --calib_detect \
    --detect_mid_strand_adapter
```
## Read trimming
```
for fq in ${data}fastq/*.fastq.gz
    do
    base=$(basename ${fq} .fastq.gz)
    echo "Trimming reads for ${base}..."
    porechop -i ${fq} -o ${data}porechop/${base}_porechop.fastq.gz --discard_middle
    for length in {3000,4000,5000,10000}
        do
        gunzip -c ${data}porechop/${base}_porechop.fastq.gz | \
        NanoLyse -r ${lambda} | \
        NanoFilt -q 10 -l ${length} > ${data}trimmed/${base}_q10_${length}trim.fastq
    done &
done
conda deactivate nanofilt
```