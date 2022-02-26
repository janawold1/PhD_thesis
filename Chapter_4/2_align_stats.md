# Alignment Statistics
The scripts provided here were augmented from Joseph Guhlin to estimate alignment stats. [Samtools v1.11](https://github.com/samtools/samtools), [mosdepth v0.3.3](https://github.com/brentp/mosdepth) and [qualimap v2.2.2](http://qualimap.conesalab.org/) software packages were used.

Aligned bam files were sorted, mates fixed and PCR duplicates were removed with SAMtools.
```
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
        ${data}nodup_bam/${base}_nodup.bam
    samtools stats ${bam} > ${data}nodup_bam_stats/${base}.stats
done
```

```
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
    mosdepth --threads 24 --fast-mode --by 50 ${data}nodup_bam_stats/${base} ${bam}
done
```
For ease of comparisons, mosdepth outputs were also plotted with:
```
python ~/anaconda3/envs/mosdepth/scripts/plot-dist.py ${data}nodup_bam_stats/*.global.dist.txt
```