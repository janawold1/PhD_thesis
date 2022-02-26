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
## Read Mapping Quality
Depths were estimated using [mosdepth v0.3.3](https://github.com/brentp/mosdepth). Plots were generated using the provided python scripts as outlined below. For future inclusion, samples had to have a coverage of >=4x with fragments >=1kb. 
```
for bam in ${data}winnowmap/alignment/bam/*.bam
    do
    indiv=$(basename $bam .bam)
    echo "Running mosdepth for $indiv"
    mosdepth --threads 24 --fast-mode --by 1000 ${data}winnowmap/alignment/mosdepth/$indiv $bam
done
cd ${data}winnowmap/alignment/mosdepth
python ~/anaconda3/envs/mosdepth/scripts/plot-dist.py *global.dist.txt
```
After review the plots, Bella and Smoko fell below the minimum depth thresholds and were removed from future analyses. 

## Initial SV calls with cuteSV and Sniffles 
Aiming to use the Jasmine pipeline, as such sniffles was run using sensitive parameters an all supporting reads were reported as recommended below:

```
for male in Bill Blades Gulliver Rangi Sass
    do
    echo "Running sniffles for ${male}..."
    sniffles --input ${data}winnowmap/alignment/bam/${male}.bam \
        --reference ${mref} \
        --vcf ${data}winnowmap/sniffles/${male}.vcf \
        --minsvlen 50 \
        --minsupport 2 \
        --threads 64 
done

for female in C1 C2 Kuia Margaret-Maree Sue
    do
    echo "Running sniffles for ${female}..."
    sniffles --input ${data}winnowmap/alignment/bam/${female}.bam \
        --reference ${fref} \
        --vcf ${data}winnowmap/sniffles/${female}.vcf \
        --minsvlen 50 \
        --minsupport 2 \
        --threads 64
done
```

## Running initial cuteSV calls

```
for bam in ${data}winnowmap/alignment/bam/*.bam
    do
    echo "Running cuteSV for ${base}..."
    if [[ $base == @(Bill|Blades|Gulliver|Rangi|Sass) ]]
        then
        cuteSV ${bam} ${mref} \
            ${data}winnowmap/cuteSV/${base}.vcf cuteSV/ \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads 24
        else
        cuteSV ${bam} ${fref} \
            ${data}winnowmap/cuteSV/${base}.vcf cuteSV/ \
            --max_cluster_bias_INS 100 \
            cuteSV ${bam} ${fref} \
            ${data}winnowmap/cuteSV/${base}.vcf cuteSV/ \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads 24 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads 24
    fi &
done
```
## Running Jasmine
After SV discovery calling, breakpoint refinement for duplication and insertions, and SV normalisation was completed.
```
ref=/kakapo-data/References/kakapo_full_ref.fa

echo "Refining sniffles calls"
cd /kakapo-data/ONT/winnowmap/sniffles/
for vcf in *.vcf
    do
    indiv=$(basename $vcf .vcf)
    echo "Running Jasmine and IRIS for $indiv SVs discovered by sniffles..."
    jasmine --dup_to_ins --preprocess_only file_list=${vcf} --comma_filelist genome_file=${ref}
    iris threads=24 vcf_in=output/${indiv}_dupToIns.vcf genome_in=${ref} reads_in=${data}winnowmap/alignment/bam/${indiv}.bam vcf_out=output/${indiv}_dupToIns_irisRefined.vcf
    echo "Normalising SV types for $indiv..."
    jasmine --preprocess_only --pre_normalize file_list=output/${indiv}_dupToIns_irisRefined.vcf --comma_filelist
done
```
After the normalisation of the files, high-confidence calls were marked. Unfortunately, the long-read data assessed here varied widely from 4x - 12x coverage. A minimum depth of 15x would have been ideal. 

The low sequence depth of our samples makes is somewhat challenging to curate high quality calls with the jasmine pipeline. This is due to the recommended minimum number of reads a variant needs to pass (```spec_reads```) is about 25% of the average coverage, which would place some samples < 2x coverage to call a high-quality variant. In light of this, a minimum ```spec_reads``` of 2 was used for all 8 samples with < 10x coverage and 3x was used for the remaining 2 samples with >=10x coverage.

```
for indiv in C1 C2 Gulliver Kuia Margaret-Maree Rangi Sass Sue
    do
    echo "Calling high-quality calls for ${indiv}..."
    jasmine file_list=output/${indiv}_dupToIns_irisRefined_normalizeTypes.vcf --comma_filelist --preprocess_only --mark_specific --spec_reads=2 spec_len=50
done

for indiv in Bill Blades
    do
    echo "Calling high-quality calls for ${indiv}..."
    jasmine file_list=output/${indiv}_dupToIns_irisRefined_normalizeTypes.vcf --comma_filelist --preprocess_only --mark_specific --spec_reads=3 spec_len=50
done

for vcf in output/*_dupToIns_irisRefined_normalizeTypes_markedSpec.vcf
    do
    indiv=$(basename $vcf _dupToIns_irisRefined_normalizeTypes_markedSpec.vcf )
    echo "Removing duplicate calls for ${indiv}..."
    jasmine file_list=$vcf --comma_filelist max_dist=200 --allow_intrasample out_file=${indiv}_jasmine.vcf --nonlinear_dist
done
```
SV calls were merged across all samples, insertions were converted back to duplications.
```
ls *_jasmine.vcf > vcf_list.txt

jasmine file_list=vcf_list.txt out_file=sniffles_mergedSVs.vcf

jasmine --dup_to_ins --postprocess_only out_file=sniffles_mergedSVs.vcf
```
The output of the conversion step of insertions back to duplications for sniffles was: ```Number of insertions converted back to duplications: 54 out of 23457 total variants```

while the conversion step output for cuteSV was: ```Number of insertions converted back to duplications: 97 out of 5263 total variants```

Finally, low-confidence / imprecise call were removed.
```
cat sniffles_mergedSVs.vcf | grep -v 'IMPRECISE;' > sniffles_mergedSVs_precise.vcf
cat sniffles_mergedSVs_precise.vcf | grep -v 'IS_SPECIFIC=0' > sniffles_mergedSVs_filtered.vcf
```
These steps were then repeated for the cuteSV calls for comparison. 

After reviewing the number of SVs per chromosome, as well as the relative mapping rates to each chromosome, autosome 1 was chosen for comparisons of a long-read alignment approach and the *de novo* assembly approach.