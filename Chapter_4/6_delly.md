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
delly merge -o ${out}01_fairy_tern_sites.bcf ${out}SVcalls/QC_pass/*.bcf
```
# SV genotyping
Sites were then merged, with a minimum size threshold of 300bp and individuals genotyped. 
```
printf "\nRunning Delly genotyping...\n"
for bcf in ${out}SVcalls/QC_pass/*.bcf
    do
    id=$(basename ${bcf} .bcf)
    printf "\nRunning Delly genotyping for ${id}..."
	delly call -g ${ref} -v ${out}01_fairy_tern_sites.bcf -o ${out}genotypes/$(id).geno.bcf /data/common_tern/alignments/nodup_bam/${id}_nodup.bam &
	done
```
## SV filtering
Finally, all genotyped individuals were merged and calls were filtered with Delly's germline setting (somatic filtering works with tumour/normal pairs). 
```
bcftools merge -m id -O b -o ${out}02_fairy_tern_genotypes.bcf ${out}genotypes/*.geno.bcf
bcftools index ${out}02_fairy_tern_genotypes.bcf

delly filter -f germline --pass -o ${out}03_fairy_tern_germline_filtered.bcf \
	${out}02_fairy_tern_genotypes.bcf

delly filter -f germline --pass -m 50 -o ${out}04_fairy_tern_minSize50bp_filtered.bcf \
	${out}02_fairy_tern_genotypes.bcf

delly filter -f germline --pass -m 300 -o ${out}05_fairy_tern_minSize300bp_filtered.bcf \
	${out}02_fairy_tern_genotypes.bcf

bcfools view -t $target -i '(SVTYPE = "INS")' \
	-O b -o ${out}INS.bcf \
	${out}03_fairy_tern_germline_filtered.bcf
bcftools view -t $target -i '(SVTYPE = "DEL")' \
    -O b -o ${out}DEL.bcf \
	${out}04_fairy_tern_minSize50bp_filtered.bcf
bcftools view -t $target -i '(SVTYPE = "INV" | SVTYPE = "DUP")' \
    -O b -o ${out}INV_DUP.bcf \
	${out}05_fairy_tern_minSize300bp_filtered.bcf

bcftools index ${out}INS.bcf
bcftools index ${out}DEL.bcf
bcftools index ${out}INV_DUP.bcf

bcftools concat -a -O v -o 06_delly_SVfilter.vcf ${out}INS.bcf ${out}DEL.bcf ${out}INV_DUP.bcf

bcftools view -i '(N_PASS(GT!="mis") = 26)' -O v -o 07_delly_genofilter.vcf 06_delly_SVfilter.vcf
```
## SV summary
Once the final SV data was filtered for, number of variable private and shared as well as fixed private and shared SVs were found as per:
```
bcftools view -S ${out}AU_samples.tsv ${out}07_delly_genofilter.vcf | \
    bcftools view -i 'GT="alt"' \
    -O z -o ${out}pops/AU_variable.vcf.gz

bcftools view -S ${out}TI_samples.tsv ${out}07_delly_genofilter.vcf | \
    bcftools view -i 'GT="alt"' \
    -O z -o ${out}pops/TI_variable.vcf.gz

bcftools index ${out}pops/AU_variable.vcf.gz
bcftools index ${out}pops/TI_variable.vcf.gz
```
The proportion of these variants fixed in each population were found with:
```
bcftools isec ${out}pops/AU_variable.vcf.gz ${out}pops/TI_variable.vcf.gz -p ${out}pops/

mv ${out}pops/0000.vcf ${out}pops/AU_private_variable.vcf
mv ${out}pops/0001.vcf ${out}pops/TI_private_variable.vcf
mv ${out}pops/0002.vcf ${out}pops/AU_shared_variable.vcf
mv ${out}pops/0003.vcf ${out}pops/TI_shared_variable.vcf

bgzip ${out}pops/AU_shared_variable.vcf
bgzip ${out}pops/TI_shared_variable.vcf
bcftools index ${out}pops/AU_shared_variable.vcf.gz
bcftools index ${out}pops/TI_shared_variable.vcf.gz

bcftools merge -m id --threads 24 -O z -o ${out}pops/shared_merged.vcf.gz ${out}pops/AU_shared_variable.vcf.gz ${out}pops/TI_shared_variable.vcf.gz

bcftools query -i 'N_PASS(GT=="AA")=15' -f '%SVTYPE\n' ${out}pops/AU_private_variable.vcf | sort | uniq -c
bcftools query -i 'N_PASS(GT=="AA")=11' -f '%SVTYPE\n' ${out}pops/TI_private_variable.vcf | sort | uniq -c
bcftools query -i 'N_PASS(GT=="AA")=26' -f '%SVTYPE\n' ${out}pops/shared_merged.vcf.gz | sort | uniq -c
```
