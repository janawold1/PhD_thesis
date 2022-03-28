#!/bin/bash/
#####################################################################################
#Script for running the MANTA software package. Intended for structural
#variant calling. MANTA is aprogram for the detection of structural variants
#from paired-end sequence data. Because running all individuals together is taxing on 
#MANTA, this program was run in 7 batches for female kakapo. 
#####################################################################################
##Setting fixed variables
ref=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/References/kakapo_full_ref.fa
work=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/alignments_female/bam/sorted/
out=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/manta/female/
MANTA_INSTALL_PATH=/scale_wlg_nobackup/filesets/nobackup/uoo02695/software/manta-1.6.0.centos6_x86_64/
#To create inputs
#ls -lh alignments_female/bam/sorted/ | awk '{print $9}' | grep -v .bai | grep bam | sort | sed -e 's%^%--bam /scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/alignments_female/bam/sorted/%' | awk '{print $0, "\"}'

##First need to create the configuration file
${MANTA_INSTALL_PATH}/bin/configManta.py \
--bam /scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/alignments_female/bam/sorted/N.sorted.bam \
--referenceFasta ${ref} \
--runDir ${out}

echo "Running MANTA workflow..."
/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/manta/female/runWorkflow.py