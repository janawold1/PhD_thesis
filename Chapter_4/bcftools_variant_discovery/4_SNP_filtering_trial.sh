#!/bin/sh
##########################################################################################################
# SNP filtering trial for tara iti and Australian fairy tern. Here we tested a number of outcomes for SNP
# filtering of tara iti specific, Australian fairy tern specific, and global SNPs. 
##########################################################################################################
input=/data/tara_iti_shortreads/
out=${input}filter_trial/

printf "\nCREATING BCFs FOR AU AND TI SAMPLES...\n"
bcftools view --threads 8 -S ${input}metadata/AU_samples.tsv -O b -o ${input}AU_VariantCalls.bcf ${input}global_VariantCalls.bcf &
bcftools view --threads 8 -S ${input}metadata/TI_samples.tsv -O b -o ${input}TI_VariantCalls.bcf ${input}global_VariantCalls.bcf
wait

for bcf in ${input}*_VariantCalls.bcf
    do
    base=$(basename ${bcf} _VariantCalls.bcf)
    for dp in {4..5}
        do
        for miss in {0.8,0.9,1}
            do
            for gq in {0,10,20}
                do
                echo "Filtering SNPs for ${base}...." 
                vcftools --bcf ${bcf} \
                    --out ${out}bcfs/${base}_${dp}x_${miss}missing_${gq}minGQ_0LD \
                    --minDP ${dp} \
                    --maxDP 100  \
                    --max-missing ${miss} \
                    --maf 0.05 \
                    --minQ 20 \
                    --minGQ ${gq} \
                    --remove-indels \
                    --remove-filtered-all \
                    --recode-bcf \
                    --recode-INFO-all
            done &
        done &
    done &
done

rename 's/.recode//g' ${out}bcfs/*.recode.bcf

printf "\nFILTERING FOR LINKAGE PARAMETERS...\n"
for bcf in ${out}bcfs/*_0LD.bcf
do
    for ld in {0.4,0.6,0.8}
    do
        base=$(basename ${bcf} _0LD.bcf)
        echo "Running light LD pruning at ${ld} for ${base}...."
        bcftools +prune \
            -l ${ld} \
            -w 1000 \
            -O b \
            -o ${out}bcfs/${base}_${ld}LD.bcf \
            ${bcf}
    done &
done

printf "\nNOW BEGINNING TO CALCULATE STATS...\n"
for bcf in ${out}bcfs/*.bcf
do
    base=$(basename ${bcf} .bcf)
    echo "Calculating depth for ${base}..."
    vcftools --bcf ${bcf} \
        --out ${work}stats/${base} \
        --site-depth &
    vcftools --vcf ${bcf} \
        --out ${work}stats/${base} \
        --depth &
    echo "Calculating missingness for ${base}..."
    vcftools --bcf ${bcf} \
        --out ${work}stats/${base} \
        --missing-site &
    vcftools --vcf ${bcf} \
        --out ${work}stats/${base} \
        --missing-indv &
    echo "Calculating individual heterozygosity for ${base}..."
    vcftools --bcf ${bcf} \
        --out ${work}stats/${base} \
        --het
    wait
done

printf "\nSUMMARISING FILTERING RESULTS...\n"
echo "sample,total_reads,mapped_reads,unmapped_reads,mismatches,average_insert_length,insert_SD,avg_qual" > ${out}mapping_stats.csv
for file in data/tara_iti_shortreads/stats/*stats
    do
    base=$(basename ${file} .stats)
    echo "appending ${base}..."
    total=$(grep "raw total sequences" ${file} | awk 'BEGIN {FS = "\t"}; {print $3}')
    mapped=$(grep "reads mapped" ${file} | grep -v "reads mapped and paired" | awk 'BEGIN {FS = "\t"}; {print $3}')
    unmapped=$(grep "reads unmapped" ${file} | awk 'BEGIN {FS = "\t"}; {print $3}')
    mismatches=$(grep "mismatches:" ${file} | awk 'BEGIN {FS = "\t"}; {print $3}')
    insert=$(grep "insert size average:" ${file} | awk 'BEGIN {FS = "\t"}; {print $3}')
    insert_sd=$(grep "insert size standard deviation:" ${file} | awk 'BEGIN {FS = "\t"}; {print $3}')
    qual=$(grep "average quality" ${file} | awk 'BEGIN {FS = "\t"}; {print $3}')
    echo "${base},${total},${mapped},${unmapped},${mismatches},${insert},${insert_sd},${qual}" >> ${out}mapping_stats.csv
done

printf "\nESTIMATING DEPTH OF BAMS...\n"
echo "sample,depth,median" > /data/tara_iti_shortreads/tara_iti_WGSdepth.csv
for bam in ${input}alignments/merged_bam/*.bam
do
    base=$(basename ${bam} .bam)
    echo "Estimating averave and median depth for ${base}..."
    depth=$(samtools depth ${bam} | awk '{sum +=$3} END {print sum/NR}') &
    median=$(samtools depth ${bam} | awk '{print $3}' | sort -n | head -n 135793321 | tail -n 1)
    wait
    echo "${base},${depth},${median}" >> ${out}tara_iti_WGSdepth.csv
done