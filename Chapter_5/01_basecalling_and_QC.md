# ONT basecalling, read quality and trimming

```

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