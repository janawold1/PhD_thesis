# SNP data filtering trial
First set up global variables:
```
data=/data/common_tern/SNP_filtering_trial/
```
Then created initial VCF for filtering trials by removing indels, and any sites with more than 2 alleles:
```
printf "\nCreating initial VCF for filtering trial...\n"
bcftools view --threads 24 -m 2 -M 2 -v snps -O z -o ${data}global.vcf.gz /data/common_tern/bcftools_variantCalls/Fairy_tern_VariantCalls.sorted.vcf.gz

```
# Filtering parameters
Biallelic SNPs for theglobal data were filtered using a range of thresholds. Specifically, the minimum depth threshold ranged from 4-6, percent missing data from 0 - 20%, and a minimum GQ score of either 0, 10, 20 or 30. A minor allele frequency of 0.05, maximum depth of 100x and minimum Q score of 20 were constant across all data sets.
```
for pop in AU TI global
    do
    vcftools --gzvcf ${data}${POP}.vcf.gz --freq2 --out ${data}stats/${POP} --max-alleles 2 &
    vcftools --gzvcf ${data}${POP}.vcf.gz --depth --out ${data}stats/${POP} &
    vcftools --gzvcf ${data}${POP}.vcf.gz --site-mean-depth --out ${data}stats/${POP} &
    vcftools --gzvcf ${data}${POP}.vcf.gz --site-quality --out ${data}stats/${POP} &
    vcftools --gzvcf ${data}${POP}.vcf.gz --missing-indv --out ${data}stats/${POP} &
    vcftools --gzvcf ${data}${POP}.vcf.gz --missing-site --out ${data}stats/${POP} &
    vcftools --gzvcf ${data}${POP}.vcf.gz --het --out ${data}stats/${POP}
done
wait



for POP in AU TI global
    do
    if (( "$POP" == global ))
        then
        for dp in {4..6}
            do
            for miss in {0.8,0.9,1}
                do
                for gq in {0,10,20,30}
                    do
                    echo "Filtering SNPs for ${POP}...."
                    vcftools --gzvcf ${data}${POP}_biallelic.vcf.gz \
                        --out ${data}vcfs/${POP}_${dp}x_${miss}missing_${gq}minGQ_0LD \
                        --minDP ${dp} \
                        --maxDP 260  \
                        --max-missing ${miss} \
                        --maf 0.05 \
                        --minQ 20 \
                        --minGQ ${gq} \
                        --remove-indels \
                        --remove-filtered-all \
                        --recode \
                        --recode-INFO-all
                done &
            done &
        done &
    else
        for dp in {4..6}
            do
            for miss in {0.8,0.9,1}
                do
                for gq in {0,10,20,30}
                    do
                    echo "Filtering SNPs for ${POP}...."
                    vcftools --gzvcf ${data}${POP}_biallelic.vcf.gz \
                        --out ${data}vcfs/${POP}_${dp}x_${miss}missing_${gq}minGQ_0LD \
                        --minDP ${dp} \
                        --maxDP 150  \
                        --max-missing ${miss} \
                        --maf 0.05 \
                        --minQ 20 \
                        --minGQ ${gq} \
                        --remove-indels \
                        --remove-filtered-all \
                        --recode \
                        --recode-INFO-all
                done &
            done &
        done &
    fi &
done
wait
```
Once these initial filtering thresholds were complete, the default ```.recode.vcf``` from each file was removed with ```rename 's/.recode//g' ${data}vcfs/*.recode.vcf```. Linkage disequilibrium was filtered using BCFtools for *r <sup>2</sup>* 0.4, 0.6 and 0.8 in 1kb windows:

```
printf "\nFILTERING FOR LINKAGE PARAMETERS...\n"
for vcf in ${data}vcfs/*_0LD.vcf
do
    for ld in {0.4,0.6,0.8}
    do
        base=$(basename ${vcf} _0LD.vcf)
        echo "Running LD pruning at ${ld} for ${base}...."
        bcftools +prune -m ${ld} -w 1000 ${vcf} -O v -o ${data}vcfs/${base}_${ld}LD.vcf
    done &
done
```

## Estimating preliminary stats

```
printf "\nNOW BEGINNING TO CALCULATE STATS...\n"
for vcf in ${data}vcfs/*.vcf
do
    base=$(basename ${vcf} .vcf)
    echo "Calculating stats for ${base}..."
    vcftools --vcf ${vcf} --out ${work}stats/${base} --site-depth &
    vcftools --vcf ${vcf} --out ${work}stats/${base} --depth &
    vcftools --vcf ${vcf} --out ${work}stats/${base} --missing-site &
    vcftools --vcf ${vcf} --out ${work}stats/${base} --missing-indv &
    vcftools --vcf ${vcf} --out ${work}stats/${base} --het
done
```

# Summarising results
Results were then summarised as below for visualisation in R.
```
printf "\nSUMMARISING FILTERING RESULTS...\n"
echo "pop,minDepth,missingness,minGQ,LD,SNP_count,mean_indiv_depth,mean_Ho,mean_He,mean_indiv_miss,mean_site_depth,mean_site_missingness" > ${data}fairy_tern_SNPsummary.csv
for stat in ${data}stats/*.log
    do
    base=$(basename $stat .log)
    pop=$(echo ${base} | tr "_" " " | awk '{print $1}')
    minDP=$(echo ${base} | tr "_" " " | awk '{print $2}' | sed 's/x//g')
    miss=$(echo ${base} | tr "_" " " | awk '{print $3}' | sed 's/missing//g')
    minGQ=$(echo ${base} | tr "_" " " | awk '{print $4}' | sed 's/minGQ//g')
    LD=$(echo ${base} | tr "_" " " | awk '{print $5}' | sed 's/LD//g')
    SNPs=$(tail -n +2 stats/${base}.ldepth | wc -l)
    indivDP=$(cat ${data}stats/${base}.idepth | awk '{sum += $3}; END {print sum/NR}')
    Ho=$(cat ${data}stats/${base}.het | tail -n +2 | awk '{print $2/$4}' | awk '{sum += $1}; END {print sum/NR}')
    He=$(cat ${data}stats/${base}.het | tail -n +2 | awk '{print $3/$4}' | awk '{sum += $1}; END {print sum/NR}')
    indivmiss=$(cat ${data}stats/${base}.imiss | awk '{sum += $5}; END {print sum/NR}')
    meansiteDP=$(cat ${data}stats/${base}.ldepth | awk '{sum += $3}; END {print sum/NR}')
    meansiteMiss=$(cat ${data}stats/${base}.lmiss | awk '{sum += $6}; END {print sum/NR}')
    echo "$pop,$minDP,$miss,$minGQ,$LD,$SNPs,$indivDP,$Ho,$He,$indivmiss,$meansiteDP,$meansiteMiss"
    echo "$pop,$minDP,$miss,$minGQ,$LD,$SNPs,$indivDP,$Ho,$He,$indivmiss,$meansiteDP,$meansiteMiss" >> ${data}fairy_tern_SNPsummary.csv
done
```

```
echo "sample,total_reads,mapped_reads,unmapped_reads,mismatches,average_insert_length,insert_SD,avg_qual" > ${data}mapping_stats.csv
for file in /data/common_tern/alignments/nodup_bam_stats/*stats
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
    echo "${base},${total},${mapped},${unmapped},${mismatches},${insert},${insert_sd},${qual}" >> ${data}mapping_stats.csv
done
```

```
printf "\nESTIMATING DEPTH OF BAMS...\n"
echo "sample,depth,median" > /data/tara_iti_shortreads/tara_iti_WGSdepth.csv
for bam in /data/common_tern/alignments/nodup_bam/*.bam
do
    base=$(basename ${bam} .bam)
    echo "Estimating average and median depth for ${base}..."
    depth=$(samtools depth ${bam} | awk '{sum +=$3} END {print sum/NR}') &
    median=$(samtools depth ${bam} | awk '{print $3}' | sort -n | head -n 135793321 | tail -n 1)
    wait
    echo "${base},${depth},${median}" >> ${data}tara_iti_WGSdepth.csv
done
```