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
Finally, SVs were sorted into population specific data sets, filtered to include fixed deletions only and then the private variants identified. 

```
bcftools view -S ${out}AU_samples.tsv -i 'SVTYPE=="DEL" & (FMT/DHFFC[0:15] < 0.7)' ${out}02_fairy_tern.genos.smoove.square.vcf.gz | \
    bcftools view -i 'N_PASS(FMT/GT =="AA") = 15' -O z -o ${out}AU_filtered.vcf.gz

bcftools view -S ${out}TI_samples.tsv -i 'SVTYPE=="DEL" & (FMT/DHFFC[15:25] < 0.7)' ${out}02_fairy_tern.genos.smoove.square.vcf.gz | \
    bcftools view -i 'N_PASS(FMT/GT =="AA") = 11' -O z -o ${out}TI_filtered.vcf.gz

bcftools view -S TI_samples.tsv -i 'SVTYPE == "DEL" & (FMT/DHFFC[15-25] < 0.7)' 02_fairy_tern.genos.smoove.square.vcf.gz | bcftools view -i 'FMT/GT[0:10] == "AA"' -O z -o TI_filtered.vcf.gz

tabix ${out}AU_filtered.vcf.gz
tabix ${out}TI_filtered.vcf.gz

bcftools isec ${out}AU_filtered.vcf.gz ${out}TI_filtered.vcf.gz -p ${out}
```