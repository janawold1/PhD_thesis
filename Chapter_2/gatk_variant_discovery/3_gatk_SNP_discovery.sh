#!/bin/sh -e
###########################################################################################
# Script for SNP discovery with GATK.
###########################################################################################
ref=/data/tara_iti_shortreads/reference/bSteHir1.pri.cur.20190820.fasta
data=/data/tara_iti_shortreads/alignments/markdup_bam/
out=/data/tara_iti_shortreads/gatk_variants/

#gatk CreateSequenceDictionary -R ${ref}

printf "\nRunning GATK HaplotypeCaller...\n"
time for i in {01..09}
	do
	for bam in ${data}batch${i}/*.bam
	do
		base=$(basename $bam _markdup.bam)
		printf "\nRunning GATK for ${base}...\n\n"
		gatk --java-options "-Xmx16g" HaplotypeCaller \
		    -I ${bam} \
		    -R ${ref} \
		    -O ${out}${base}.gvcf.gz \
		    -ERC GVCF
	done &
done

#printf "\nConverting GVCF to BCF..."
#bcftools convert \
#    --gvcf2vcf \
#    -f ${ref} \
#    -O bcf \
#    -o ${bcf} \
#    ${gvcf}
#bcftools index ${bcf}
#bcftools stats --threads 16 -F ${ref} ${bcf} > ${out}fairy_tern_identity_unfiltered.stats
