# Script for SV discovery with Delly.
Here is a brief walkthrough of how I ran Delly for SV discovery in Australian fairy tern (*Sterna nereis nereis*) and tara iti (*Sterna nereis davisae*) populations. It is notable that we are looking at reads aligned to a conspecific reference genome (Common tern, *Sterna hirundo*). In light of this, we aim to identify Deletions with a focus on those that are fixed and private in one population or the other. 

To start, we defined global variables as with ANGSD script. 
```
work=/data/alignments/nodup_bam/
out=/data/delly/
ref=/data/reference/bSteHir1.pri.cur.20190820.fasta
bcf=${out}fairy_tern_SVsites.bcf
target=/data/metadata/common_tern_autosomes.bed
TI=/data/metadata/TI_samples.tsv
AU=/data/metadata/AU_samples.tsv
```
## SV Discovery
Putative SV sites were identified as below.
```
printf "\nRunning Delly...\n"
time for i in {01..09}
	do
	for bam in ${work}batch${i}/*_nodup.bam
		do
		base=$(basename ${bam} _nodup.bam)
		printf "\nRunning Delly for ${base}...\n"
		delly call -g ${ref} -o ${out}SVcalls/${base}.bcf  ${bam}
	done &
done
wait
delly merge -o ${out}01_fairy_tern_minSize300bp_sites.bcf -m 300 SVcalls/QC_pass/*.bcf
```
# SV genotyping
Sites were then merged, with a minimum size threshold of 300bp and individuals genotyped. 
```
printf "\nRunning Delly genotyping...\n"
while read -r line
    do
    id=$(echo $line | tr "/" " " | awk '{print $5}' | sed 's/_nodup.bam//')
    printf "\nRunning Delly genotyping for ${mbase}..."
    delly call -g ${ref} -v ${out}01_fairy_tern_minSize300bp_sites.bcf -o ${out}genotypes/${mbase}.geno.minsize.bcf ${line}
done
```
## SV filtering
Finally, all genotyped individuals were merged and calls were filtered with Delly's germline setting (somatic filtering works with tumour/normal pairs). Deletions were filtered for those that passed all record level filters in each population. Finally, the interesection between the two populations was used to identify private deletions. 
```
bcftools merge -m id -O b -o ${out}02_fairy_tern_genotypes.bcf ${out}genotypes/*.bcf
bcftools index ${out}02_fairy_tern_genotypes.bcf

delly filter -f germline -o ${out}03_fairy_tern_germline_filtered.bcf \
	${out}02_fairy_tern_genotypes.bcf

bcftools view -S ../smoove/AU_samples.tsv -i 'FILTER=="PASS" & SVTYPE == "DEL"' 03_fairy_tern_germline_filtered.bcf | bcftools view -i 'N_PASS(GT == "AA") = 15' -O z -o AU_filtered.vcf.gz

bcftools view -S ../smoove/TI_samples.tsv -i 'FILTER=="PASS" & SVTYPE == "DEL"' 03_fairy_tern_germline_filtered.bcf | bcftools view -i 'N_PASS(GT == "AA") = 11' -O z -o TI_filtered.vcf.gz

tabix ${out}${pop}_fixed_DEL.vcf.gz

bcftools isec ${out}TI_fixed_DEL.vcf.gz ${out}AU_fixed_DEL.vcf.gz -p ${out}
```
Now we are left with the fixed deletions that are private to each population and shared, spoiler alert there are none.