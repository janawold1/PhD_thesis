# Running Smoove to identify fixed deletions in tara iti and Australian fairy tern populations
Defining global variables:
```
ref=/data/reference/bSteHir1.pri.cur.20190820.fasta
data=/data/common_tern/alignments/nodup_bam/
out=/kakapo-data/bwa/smoove/
remove=/data/common_tern/metadata/
TMPDIR=/data/common_tern/smoove/temp/
```
## SV discovery
Similar to SV discovery in Chapter 3, Smoove was run as below:
```
ulimit -Sn 5000
for i in {01..11}
    do
    for fbam in ${data}bwa_female/bam/batch${i}/*.bam
        do
        fbase=$(basename ${fbam} .sorted.bam)
        echo "Running SMOOVE call for ${fbase}..."
        smoove call --name ${fbase} --fasta ${ref} --outdir ${out}SV_calls_female \
            -p 1 --genotype ${fbam}
    done &
done
wait
```
## Merging variant calls and genotyping output
Called variants were then merged and genotyped. 
```
echo "Merging all called variants..."
smoove merge --name fairy_tern -f ${ref} --outdir ${out} ${out}SV_calls/QC_pass/*.vcf.gz
while read -r line
    do
    id=$(echo $line | tr "/" " " | awk '{print $5}' | sed 's/_nodup.bam//g')
    echo "Creating individual genotyped VCF for ${id}...."
    smoove genotype -d -x -p 1 --name ${id} --fasta ${ref} --outdir ${out}genotypes \
        --duphold --vcf ${out}01_fairy_tern.sites.vcf.gz ${line}
done < ${out}samples.bamlist
wait
echo "Creating total raw VCF..."
smoove paste --name ${out}02_fairy_tern.genos ${data}genotypes/*.vcf.gz
```

## SV Filtering
Finally, SVs were sorted into population specific data sets, filtered to include variable SVs only and then the number of private and shared variants were identified. 
```
bcftools view -t chroms -i '(SVTYPE = "DEL" & FMT/DHFFC[0-25] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0-25] > 1.3) | (SVTYPE = "INV")' \
    -O z -o ${out}03_fairy_tern_SVfiltered.vcf.gz ${out}02_fairy_tern.genos.smoove.square.vcf.gz

bcftools view -i 'N_PASS(GT="mis")=0' -O ${out}04_fairy_tern_filtered_missing.vcf.gz ${out}03_fairy_tern_SVfiltered.vcf.gz

bcftools view -S ${out}AU_samples.tsv ${out}04_fairy_tern_filtered_missing.vcf.gz \
    bcftools view -i 'GT="alt"' \
    -O z -o ${out}pops/AU_variable.vcf.gz

bcftools view -S ${out}TI_samples.tsv ${out}04_fairy_tern_filtered_missing.vcf.gz \
    bcftools view -i 'GT="alt")' \
    -O z -o ${out}pops/TI_variable.vcf.gz

bcftools index ${out}pops/AU_variable.vcf.gz
bcftools index ${out}pops/TI_variable.vcf.gz
```
The proportion of these variants fixed in each population were found with:
```
mv ${out}pops/0000.vcf ${out}pops/AU_private_variable.vcf
mv ${out}pops/0001.vcf ${out}pops/TI_private_variable.vcf
mv ${out}pops/0003.vcf ${out}pops/TI_shared_variable.vcf
mv ${out}pops/0002.vcf ${out}pops/AU_shared_variable.vcf

bgzip ${out}pops/AU_shared_variable.vcf
bcftools index ${out}pops/AU_shared_variable.vcf.gz
bgzip ${out}pops/TI_shared_variable.vcf
bcftools index ${out}pops/TI_shared_variable.vcf.gz

bcftools merge -m id --threads 24 -O z -o shared_merged.vcf.gz AU_shared_variable.vcf.gz TI_shared_variable.vcf.gz

bcftools query -i 'N_PASS(GT=="AA")=15' -f '%SVTYPE\n' ${out}pops/AU_private_variable.vcf | sort | uniq -c
bcftools query -i 'N_PASS(GT=="AA")=11' -f '%SVTYPE\n' ${out}pops/TI_private_variable.vcf | sort | uniq -c
bcftools query -i 'N_PASS(GT=="AA")=26' -f '%SVTYPE\n' ${out}pops/shared_merged.vcf.gz | sort | uniq -c
```
## SV Summary

```

```