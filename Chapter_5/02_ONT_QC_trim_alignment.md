# ONT long-read alignment with minimap2

```
data=/kakapo-data/ONT/
mref=/kakapo-data/References/kakapo_no_Wchromosome.fa
fref=/kakapo-data/References/kakapo_full_ref.fa
lambda=${data}DNA_CS.fasta
aves=aves_odb10
source ~/anaconda3/etc/profile.d/conda.sh
```

## ONT read alignment with minimap2

```
for male in Bill Blades Gulliver Rangi Sass Smoko
    do
    echo "Running minmap2 alignment for ${male} ONT trimmed reads..."
    minimap2 -ax map-ont ${mref} ${data}trimmed/${male}_q10_5kbtrim.fastq.gz > ${data}minimap/alignment/sam/${male}.sam &
done
wait

for female in Bella C1 C2 Kuia Margaret-Maree Sue
    do
    echo "Running minimap2 alignment for ${female} ONT trimmed reads..."
    minimap2 -ax map-ont ${fref} ${data}trimmed/${female}_q10_5kbtrim.fastq.gz > ${data}minimap/alignment/sam/${female}.sam &
done
wait

for male in Bill Blades Gulliver Rangi Sass Smoko
    do
    echo "Converting SAM to BAM for ${male}..."
    samtools view -@ 64 -T ${mref} -b ${data}minimap/alignment/sam/${male}.sam | \
        samtools sort -@ 64 -o ${data}minimap/alignment/bam/${male}.bam
done

for female in Bella C1 C2 Kuia Margaret-Maree Sue
    do
    echo "Converting SAM to BAM for ${female}..."
    samtools view -@ 64 -T ${fref} -b ${data}minimap/alignment/sam/${female}.sam | \
        samtools sort -@ 64 -o ${data}minimap/alignment/bam/${female}.bam
done
```
Chromosome coverage was then visualised as per:
```
for bam in ${data}minimap/alignment/bam/*.bam
    do
    base=$(basename $bam .bam)
    echo "Estimating coverage for $base..."
    samtools coverage $bam -o ${base}_coverage.tsv &
done
```
## ONT read alignment with winnowmap

```
cd ${data}winnowmap

meryl count k=15 ${fref} output kakapo_full.meryl
meryl count k=15 ${mref} output kakapo_noW.meryl

meryl print greater-than distinct=0.9998 kakapo_full.meryl > kakapo_full_repetitive_k15.txt
meryl pring greater-than distinct=0.9998 kakapo_noW.meryl > kakapo_noW_repetitive_k15.txt

for male in Bill Blades Gulliver Rangi Sass Smoko
    do
    winnowmap -W kakapo_noW_repetitive_k15.txt -ax map-ont ${mref} ${data}trimmed/${male}_q10_5kbtrim.fastq.gz > sam/${male}.sam
    samtools view -@64 -T ${mref} -b sam/${male}.sam | samtools sort -@ 64 -o bam/${male}.bam
done
```

## Calling SVs with Sniffles
Aiming to use the Jasmine pipeline, as such sniffles was run using sensitive parameters an all supporting reads were reported as recommended below:

```
for male in Bill Blades Gulliver Rangi Sass Smoko
    do
    echo "Running sniffles for ${male}..."
    sniffles --input ${data}winnowmap/alignment/bam/${male}.bam \
        --reference ${mref} \
        --vcf ${data}winnowmap/sniffles/${male}.vcf \
        --threads 64 
done

for female in Bella C1 C2 Kuia Margaret-Maree Sue
    do
    echo "Running sniffles for ${female}..."
    sniffles --input ${data}winnowmap/alignment/bam/${female}.bam \
        --reference ${fref} \
        --vcf ${data}winnowmap/sniffles/${female}.vcf \
        --threads 64
done
```

## Running initial cuteSV calls

```
for bam in ${data}winnowmap/alignment/bam/*.bam
    do
    echo "Running cuteSV for ${base}..."
    if [[ $base == @(Bill|Blades|Gulliver|Rangi|Sass|Smoko) ]]
        then
        cuteSV ${data}winnowmap/alignment/bam/${base}.bam ${mref} \
            ${data}winnowmap/cuteSV/${base}.vcf cuteSV/ \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads 24
        else
        cuteSV ${data}winnowmap/alignment/bam/${base}.bam ${fref} \
            ${data}winnowmap/cuteSV/${base}.vcf cuteSV/ \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads 24
    fi &
done
```