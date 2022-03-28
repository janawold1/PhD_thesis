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


bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0_fail_smoove_genofilter\n' \
    ${out}05_smoove_genofiltered_trio.vcf >> ${out}smoove_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.05' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t<=0.05_fail_smoove_genofilter\n' \
    ${out}05_smoove_genofiltered_trio.vcf >> ${out}smoove_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.1' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t<=0.1_fail_smoove_genofilter\n' \
    ${out}05_smoove_genofiltered_trio.vcf >> ${out}smoove_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.2' -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t<=0.2_fail_smoove_genofilter\n' \
    ${out}05_smoove_genofiltered_trio.vcf >> ${out}smoove_mendel.tsv
```

## Lineage Comparisons

```
mkdir -p ${out}lineage_comparisons/{SVfiltered,unfiltered}

bcftools view -s M ${out}01_smoove_unfiltered_trio.vcf | bcftools view -i 'GT=="alt"' -O z -o ${out}lineage_comparisons/unfiltered/RH_unfiltered_variants.vcf.gz
bcftools view -s ^M,G,F,R,N,P,O,S ${out}01_smoove_unfiltered.vcf | bcftools view -i 'GT!="RR" & GT!="mis"' -O z -o ${out}lineage_comparisons/SI_variants.vcf.gz

bcftools index ${out}lineage_comparisons/unfiltered/RH_unfiltered_variants.vcf.gz
bcftools index ${out}lineage_comparisons/unfiltered/SI_unfiltered_variants.vcf.gz

bcftools isec ${out}lineage_comparisons/unfiltered/RH_unfiltered_variants.vcf.gz \
    ${out}lineage_comparisons/unfiltered/SI_unfiltered_variants.vcf.gz \
    -p ${out}lineage_comparisons/unfiltered/
```
Summarising numbers of SVs per individual
```
bcftools view -h $${out}05_delly_genofiltered_trio.vcf | grep CHROM | tr "\t" "\n" | tail -n 169 > ${out}samples.txt
bcftools query -f '%CHROM\t%POS\n' ${out}lineage_comparisons/0000.vcf > ${out}lineage_comparisons/Fiordland_unfiltered_private_sites.txt
bcftools query -f '%CHROM\t%POS\n' ${out}lineage_comparisons/0001.vcf > ${out}lineage_comparisons/SI_unfiltered_private_sites.txt
bcftools query -f '%CHROM\t%POS\n' ${out}lineage_comparisons/0002.vcf > ${out}lineage_comparisons/shared_unfiltered_sites.txt

while read -r line
    do
    indiv=$(echo $line | awk '{print $1}')
    gen=$(echo $line | awk '{print $2}')
    echo "Counting SVs for $indiv..."
    bcftools view -s ${indiv} -R lineage_comparisons/unfiltered/Fiordland_unfiltered_private_sites.txt 01_smoove_unfiltered.vcf.gz | \
        bcftools query -i 'GT=="alt"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVTYPE\tFiordland_unfiltered_lineage\n' >> smoove_lineage_counts.tsv
    bcftools view -s ${indiv} -R lineage_comparisons/unfiltered/Rakiura_unfiltered_private_sites.txt 01_smoove_unfiltered.vcf.gz | \
        bcftools query -i 'GT=="alt"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVTYPE\tRakiura_unfiltered_lineage\n' >> smoove_lineage_counts.tsv
    bcftools view -s ${indiv} -R lineage_comparisons/unfiltered/Shared_unfiltered_sites.txt 01_smoove_unfiltered.vcf.gz | \
        bcftools query -i 'GT=="alt"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVTYPE\tShared_unfiltered_lineage\n' >> smoove_lineage_counts.tsv
done < /kakapo-data/metadata/generations.tsv

while read -r line
    do
    indiv=$(echo $line | awk '{print $1}')
    gen=$(echo $line | awk '{print $2}')
    echo "Counting SVs for $indiv..."
    bcftools view -s ${indiv} ${out}05_smoove_genofiltered_trio.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\n' >> ${out}smoove_indiv_counts.tsv
    echo "Assigning $gen generation to $indiv"
    grep "^$indiv" ${out}smoove_indiv_counts.tsv | awk -v var="$gen" '{print $0"\t"var}' >> ${out}smoove_generations.tsv
done < /kakapo-data/metadata/generations.tsv
```