# Running the SMOOVE software package

First set up global variables and a ulimit size of 5000. Smoove initially uses many small files in early SV discovery and can crash as the number of temporary files created can excede the system default limits. 
```
ref=/kakapo-data/References/kakapo_full_ref.fa
data=/kakapo-data/bwa/
out=/kakapo-data/bwa/smoove/
remove=/kakapo-data/bwa/smoove/unplaced_scaffolds.bed
annotate=/kakapo-data/metadata/annotation/GCF_004027225.2_bStrHab1.2.pri_genomic.gff.gz
TMPDIR=/kakapo-data/bwa/smoove/temp/
ulimit -Sn 5000
```
## Sex specific SV discovery
```
for i in {01..11}
    do
    for fbam in ${data}bwa_female/bam/batch${i}/*.bam
        do
        fbase=$(basename ${fbam} .sorted.bam)
        echo "Running SMOOVE call for ${fbase}..."
        smoove call --name ${fbase} --fasta ${ref} --outdir ${out}SV_calls_female \
            --exclude ${remove} -p 1 --genotype ${fbam}
    done &
done
for j in {01..14}
    do
    for mbam in ${data}bwa_male/bam/batch${j}/*.bam
        do
        mbase=$(basename ${mbam} .sorted.bam)
        echo "Running SMOOVE call for ${mbase}..."
        smoove call --name ${mbase} --fasta ${ref} --outdir ${out}SV_calls_male \
            --exclude ${remove} -p 1 --genotype ${mbam}
    done &
done
wait
```

```
echo "Merging all called variants..."
smoove merge --name bwa_smoove -f ${ref} --outdir ${out} ${out}SV_calls_male/*.genotyped.vcf.gz ${out}SV_calls_female/*.genotyped.vcf.gz
```
## Genotyping

```
for i in {01..11}
    do
    for fbam in ${data}bwa_female/bam/batch${i}/*.bam
        do
        fbase=$(basename ${fbam} .sorted.bam)
        echo "Creating individual genotyped VCF for ${fbase}...."
        smoove genotype -d -x -p 1 --name ${fbase} --fasta ${ref} --outdir ${out}genotypes \
            --duphold --vcf ${out}bwa_smoove.sites.vcf.gz ${fbam}
    done &
done
for j in {01..14}
    do
    for mbam in ${data}bwa_male/bam/batch${j}/*.bam
        do
        mbase=$(basename ${mbam} .sorted.bam)
        echo "Creating individual genotyped VCF for ${mbase}..."
        smoove genotype -d -x -p 1 --name ${mbase} --fasta ${ref} --outdir ${out}genotypes \
            --duphold --vcf ${out}bwa_smoove.sites.vcf.gz ${mbam}
    done &
done
wait
```
## Initial results

```
echo "Creating total raw VCF..."
smoove paste --name ${out}bwa_smoove.genos ${data}genotypes/*.vcf.gz
echo "Annotating raw VCF..."
smoove annotate --gff ${annotate} ${out}bwa_smoove.genos.smoove.square.vcf.gz | bgzip -c > ${out}bwa_smoove.annotated.vcf.gz
```


## Filtering

```
Used bwa_smoove_annotated.vcf.gz for unfiltered SV summary stats. Created filtered file as per:
bcftools view -t ^NC_044302.2 -O v -o ${out}01_smoove_unfiltered.vcf ${out}bwa_smoove.annotated.vcf.gz
bcftools view -i '(SVTYPE = "DEL" & FMT/DHFFC[0-168] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0-168] > 1.3) | (SVTYPE = "INV")' \
    -O v -o ${out}02_smoove_SVfiltered.vcf ${out}01_smoove_unfiltered.vcf
bcftools view -i '(MSHQ>=3)' -O v -o ${out}03_smoove_genofiltered.vcf ${out}02_smoove_SVfiltered.vcf
```

### Initial summaries
```
bcftools query -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tsmoove_unfiltered\n' ${out}01_smoove_unfiltered.vcf > ${out}smoove_summary.tsv
bcftools query -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tsmoove_SVfiltered\n' ${out}02_smoove_SVfiltered.vcf >> ${out}smoove_summary.tsv
bcftools query -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tsmoove_genofiltered\n' ${out}03_smoove_genofiltered.vcf >> ${out}smoove_summary.tsv
```

## Mendelian inheritance tests

```
bcftools +mendelian -m a -T ${trio} -O v -o ${out}04_smoove_SVfiltered_trio.vcf \
    ${out}02_smoove_SVfiltered.vcf
bcftools +mendelian -m a -T ${trio} -O v -o ${out}05_smoove_genofiltered_trio.vcf \
    ${out}03_smoove_genofiltered.vcf

bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0_fail_smoove_SVfilter\n' ${out}04_smoove_SVfiltered_trio.vcf >> ${out}smoove_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.05' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t<=0.05_fail_smoove_SVfilter\n' \
    ${out}04_smoove_SVfiltered_trio.vcf >> ${out}smoove_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.1' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t<=0.1_fail_smoove_SVfilter\n' \
    ${out}04_smoove_SVfiltered_trio.vcf >> ${out}smoove_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.2' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t<=0.2_fail_smoove_SVfilter\n' \
    ${out}04_smoove_SVfiltered_trio.vcf >> ${out}smoove_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0_fail_smoove_genofilter\n' \
    ${out}05_delly_genofiltered_trio.vcf >> ${out}smoove_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.05' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t<=0.05_fail_smoove_genofilter\n' \
    ${out}05_delly_genofiltered_trio.vcf >> ${out}smoove_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.1' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t<=0.1_fail_smoove_genofilter\n' \
    ${out}05_delly_genofiltered_trio.vcf >> ${out}smoove_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.2' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t<=0.2_fail_smoove_genofilter\n' \
    ${out}05_delly_genofiltered_trio.vcf >> ${out}smoove_mendel.tsv
```

## Lineage Comparisons

```
mkdir -p ${out}lineage_comparisons

bcftools view -s Richard_Henry ${out}05_smoove_genofiltered_trio.vcf | bcftools view -i 'GT!="RR" & GT!="mis"' -O z -o ${out}lineage_comparisons/RH_variants.vcf.gz
bcftools view -s ^Richard_Henry,Kuia,Gulliver,Sinbad,Adelaide,Henry,Marian,Gertrude ${out}05_delly_genofiltered_trio.vcf | bcftools view -i 'GT!="RR" & GT!="mis"' -O z -o ${out}lineage_comparisons/SI_variants.vcf.gz

tabix ${out}lineage_comparisons/RH_variants.vcf.gz
tabix ${out}lineage_comparisons/SI_variants.vcf.gz

bcftools isec ${out}lineage_comparisons/RH_variants.vcf.gz \
    ${out}lineage_comparisons/SI_variants.vcf.gz \
    -p ${out}lineage_comparisons/
```
Summarising numbers of SVs per individual
```
bcftools view -h $${out}05_delly_genofiltered_trio.vcf | grep CHROM | tr "\t" "\n" | tail -n 169 > ${out}samples.txt
bcftools query -f '%CHROM\t%POS\n' ${out}lineage_comparisons/0000.vcf > ${out}lineage_comparisons/Fiordland_sites.txt
bcftools query -f '%CHROM\t%POS\n' ${out}lineage_comparisons/0001.vcf > ${out}lineage_comparisons/SI_sites.txt
bcftools query -f '%CHROM\t%POS\n' ${out}lineage_comparisons/0002.vcf > ${out}lineage_comparisons/shared_sites.txt
```
while read -r line
    do
    echo "Counting SVs for ${line}..."
    si_del=$(bcftools view -T ${out}lineage_comparisons/SI_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="DEL" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    si_dup=$(bcftools view -T ${out}lineage_comparisons/SI_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="DUP" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    si_ins=$(bcftools view -T ${out}lineage_comparisons/SI_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="INS" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    si_inv=$(bcftools view -T ${out}lineage_comparisons/SI_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="INV" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    sh_del=$(bcftools view -T ${out}lineage_comparisons/shared_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="DEL" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    sh_dup=$(bcftools view -T ${out}lineage_comparisons/shared_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="DUP" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    sh_ins=$(bcftools view -T ${out}lineage_comparisons/shared_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="INS" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    sh_inv=$(bcftools view -T ${out}lineage_comparisons/shared_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="INV" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    rh_del=$(bcftools view -T ${out}lineage_comparisons/Fiordland_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="DEL" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    rh_dup=$(bcftools view -T ${out}lineage_comparisons/Fiordland_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="DUP" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    rh_ins=$(bcftools view -T ${out}lineage_comparisons/Fiordland_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="INS" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    rh_inv=$(bcftools view -T ${out}lineage_comparisons/Fiordland_sites.txt -s ${line} ${out}bwa_smoove_filtered_trios.vcf | bcftools query -i 'SVTYPE=="INV" & GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    printf "${line}\tDeletions\t${si_del}\tsmoove_SI\n${line}\tDuplications\t${si_dup}\tsmoove_SI\n${line}\tInsertions\t${si_ins}\tsmoove_SI\n${line}\tInversions\t${si_inv}\tsmoove_SI\n" >> ${out}lineage_comparisons/smoove_indiv_SVtypecounts.tsv
    printf "${line}\tDeletions\t${sh_del}\tsmoove_shared\n${line}\tDuplications\t${sh_dup}\tsmoove_shared\n${line}\tInsertions\t${sh_ins}\tsmoove_shared\n${line}\tInversions\t${sh_inv}\tsmoove_shared\n" >> ${out}lineage_comparisons/smoove_indiv_SVtypecounts.tsv
    printf "${line}\tDeletions\t${rh_del}\tsmoove_RH\n${line}\tDuplications\t${rh_dup}\tsmoove_RH\n${line}\tInsertions\t${rh_ins}\tsmoove_RH\n${line}\tInversions\t${rh_inv}\tsmoove_RH\n" >> ${out}lineage_comparisons/smoove_indiv_SVtypecounts.tsv
done < ${out}samples.txt

while read -r line
    do
    echo "Counting SVs for ${line}..."
    si=$(bcftools view -T ${out}lineage_comparisons/SI_sites.txt -s ${line} ${out}05_delly_genofiltered_trio.vcf  | bcftools query -i 'GT!="RR" & GT!="mis"' -f '%SVTYPE\t%SVLEN\n' | wc -l)
    sh=$(bcftools view -T ${out}lineage_comparisons/shared_sites.txt -s ${line} ${out}05_delly_genofiltered_trio.vcf  | bcftools query -i 'GT!="RR" & GT!="mis"' -f '%SVTYPE\t%SVLEN\n' | wc -l)
    rh=$(bcftools view -T ${out}lineage_comparisons/Fiordland_sites.txt -s ${line} ${out}05_delly_genofiltered_trio.vcf  | bcftools query -i 'GT!="RR" & GT!="mis"' -f '%SVTYPE\t%SVLEN\n' | wc -l)
    printf "${line}\t${si}\tsmoove_SI\n" >> ${out}lineage_comparisons/smoove_indiv_SVcounts.tsv
    printf "${line}\t${sh}\tsmoove_shared\n" >> ${out}lineage_comparisons/smoove_indiv_SVcounts.tsv
    printf "${line}\t${rh}\tsmoove_RH\n" >> ${out}lineage_comparisons/smoove_indiv_SVcounts.tsv
    bcftools view -T ${out}lineage_comparisons/SI_sites.txt -s ${line} ${out}05_delly_genofiltered_trio.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tsmoove_SI\n' >> detailed_smoove_summary.tsv
    bcftools view -T ${out}lineage_comparisons/shared_sites.txt -s ${line} ${out}05_delly_genofiltered_trio.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tsmoove_shared\n' >> detailed_smoove_summary.tsv
    bcftools view -T ${out}lineage_comparisons/Fiordland_sites.txt -s ${line} ${out}05_delly_genofiltered_trio.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tsmoove_RH\n' >> detailed_smoove_summary.tsv
done < ${out}samples.txt

while read -r line
    do
    indiv=$(echo $line | awk '{print $1}')
    gen=$(echo $line | awk '{print $2}')
    echo "Counting SVs for $indiv..."
    bcftools view -s ${indiv} ${out}05_smoove_genofiltered_trio.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\n' >> ${out}smoove_indiv_counts.tsv
    echo "Assigning $gen generation to $indiv"
    grep "^$indiv" ${out}smoove_indiv_counts.tsv | awk -v var="$gen" '{print $0"\t"var}' >> ${out}smoove_generations.tsv
done < /kakapo-data/metadata/generations.tsv