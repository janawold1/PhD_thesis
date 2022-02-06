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
```
Sites were then merged, with a minimum size threshold of 300bp. And individuals genotyped. 
```
delly merge -o ${out}fairy_tern_delly_sites_minsize.bcf --minsize 300 ${out}SVcalls/*.bcf
printf "\nRunning Delly genotyping...\n"
time for i in {01..09}
	do
	for bam in ${work}batch${i}/*_markdup.bam
		do
		base=$(basename ${bam} .bam)
		printf "\nRunning Delly genotyping for ${base}..."
		delly call -g ${ref} -v ${out}fairy_tern_delly_sites_minsize.bcf -o ${out}genotypes/${base}.minsize.geno.bcf ${bam}
	done &
done
```
Finally, all genotyped individuals were merged and calls were filtered with Delly's germline setting (somatic filtering works with tumour/normal pairs). Deletions were filtered for those that passed all record level filters in each population. Finally, the interesection between the two populations was used to identify private deletions. 
```
bcftools merge -m id -O b -o ${out}fairy_tern_delly_minsize_genotypes.bcf \
	${out}genotypes/*.minsize.geno.bcf

delly filter -f germline -o ${out}fairy_tern_delly_minsize_genotypes_filtered.bcf \
	${out}delly_minsize_genotypes.bcf

for group in /data/metadata/*_samples.tsv
	do
	pop=$(basename $group _samples.tsv)
	echo "Identifying fixed deletions for $pop..."
	if [[ $pop == TI]]
		then
		bcftools view -T $target -S $group -i 'N-\
			-O z -o ${out}${pop}_fixed_DEL.vcf.gz \
			${out}delly_minsize_genotypes.bcf
		tabix ${out}${pop}_fixed_DEL.vcf.gz
		else
		bcftools view -T $target -S $group -i 'N_PASS(GT=="AA")>=19 & FILTER == "PASS"' \
			-O z -o ${out}${pop}_fixed_DEL.vcf.gz \
			${out}delly_minsize_genotypes.bcf
		tabix ${out}${pop}_fixed_DEL.vcf.gz
	fi
done

bcftools isec ${out}TI_fixed_DEL.vcf.gz ${out}AU_fixed_DEL.vcf.gz -p ${out}
```
Now we are left with the fixed deletions that are private to each population and shared, spoiler alert there are none.