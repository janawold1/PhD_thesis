#!/bin/bash/ -e
#####################################################################################
#Script for running the Delly software package. Intended for structural variant
#calling. Delly is a program for the detection of structural variants from paired-end
#sequence data.
#####################################################################################

##Setting fixed variables
workbam=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/alignments_male/bam/sorted/
ref_noW=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/References/kakapo_no_Wchromosome.fa
control=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/alignments_male/jane/jane_sim_noW.sorted.bam
output=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/delly/germline_calling/male/

for samp in ${workbam}*.sorted.bam
        do
        base=$(basename ${samp} .sorted.bam)
        echo "Running Delly for ${base}"
        time delly call -g ${ref_noW} -o ${output}${base}.bcf ${workbam}${base}.sorted.bam
done

###Merging SV sites into a unified site list

samps=$(ls -lh ${output}*.bcf | awk '{print $9}' | xargs)
export LC_ALL=C; unset LANGUAGE #Known issue in Delly where you have to reset the locale prior to running 'merge'
time delly merge -o ${output}male_sites.bcf ${samps}

###Gentoyping the merged SV site list across all samples. Can be run in parallel for each sample.
for samp in ${workbam}*.sorted.bam
        do
        base=$(basename ${samp} .sorted.bam)
        echo "Genotyping SVs for ${base}"
        time delly call -g ${ref_noW} -v ${output}male_sites.bcf -o ${base}.geno.bcf ${samp}
done

###Merge all the genotyped samples into single BCF
samps=$(ls -lh ${output}*.geno.bcf | awk '{print $9}' | xargs)
time bcftools merge -m id -O b -o ${output}merged_germline_male.bcf ${samps}

###Apply the Delly germline SV filter

time delly filter -f germline -o ${output}merged_germline_male_SVfiltered.bcf ${output}merged_germline_male.bcf