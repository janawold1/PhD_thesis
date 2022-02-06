# Running Smoove to identify fixed deletions in tara iti and Australian fairy tern populations
Defining global variables:
```
ref=/kakapo-data/References/kakapo_full_ref.fa
data=/kakapo-data/bwa/
out=/kakapo-data/bwa/smoove/
remove=/kakapo-data/bwa/smoove/unplaced_scaffolds.bed
annotate=/kakapo-data/metadata/annotation/GCF_004027225.2_bStrHab1.2.pri_genomic.gff.gz
TMPDIR=/kakapo-data/bwa/smoove/temp/
```
# SV discovery
```
ulimit -Sn 5000
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
## Merging variant calls, genotyping and annotating output
```
echo "Merging all called variants..."
smoove merge --name bwa_smoove -f ${ref} --outdir ${out} ${out}SV_calls_male/*.genotyped.vcf.gz ${out}SV_calls_female/*$for i in {01..11}
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
echo "Creating total raw VCF..."
smoove paste --name ${out}bwa_smoove.genos ${data}genotypes/*.vcf.gz
echo "Annotating raw VCF..."
smoove annotate --gff ${annotate} ${out}bwa_smoove.genos.smoove.square.vcf.gz | bgzip -c > ${out}bwa_smoove.annotated.vcf.gz
```
# Used bwa_smoove_annotated.vcf.gz for unfiltered SV summary stats. Created filtered file as per:
bcftools view -t ^NC_044302.2 -O v -o ${out}01_smoove_unfiltered.vcf ${out}bwa_smoove.annotated.vcf.gz
bcftools view -i '(SVTYPE = "DEL" & FMT/DHFFC[0-168] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0-168] > 1.3) | (SVTYPE = "INV")' \
    -O v -o ${out}02_smoove_SVfiltered.vcf ${out}01_smoove_unfiltered.vcf
bcftools view -i '(MSHQ>=3)' -O v -o ${out}03_smoove_genofiltered.vcf ${out}02_smoove_SVfiltered.vcf

bcftools query -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tsmoove_unfiltered\n' ${out}01_smoove_unfiltered.vcf > ${out}smoove_summary.tsv
bcftools query -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tsmoove_SVfiltered\n' ${out}02_smoove_SVfiltered.vcf >> ${out}smoove_summary.tsv
bcftools query -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tsmoove_genofiltered\n' ${out}03_smoove_genofiltered.vcf >> ${out}smoove_summary.tsv


# Lineage Comparisons

mkdir -p ${out}lineage_comparisons

bcftools view -s Richard_Henry ${out}05_smoove_genofiltered_trio.vcf | bcftools view -i 'GT!="RR" & GT!="mis"' -O z -o ${out}lineage_comparisons/RH_variants.vcf.gz
bcftools view -s ^Richard_Henry,Kuia,Gulliver,Sinbad,Adelaide,Henry,Marian,Gertrude ${out}05_delly_genofiltered_trio.vcf | bcftools view -i 'GT!="RR" & GT!="mis"' -O z -o ${out}lineage_comparisons/SI_variants.vcf.gz

tabix ${out}lineage_comparisons/RH_variants.vcf.gz
tabix ${out}lineage_comparisons/SI_variants.vcf.gz

bcftools isec ${out}lineage_comparisons/RH_variants.vcf.gz \
    ${out}lineage_comparisons/SI_variants.vcf.gz \
    -p ${out}lineage_comparisons/