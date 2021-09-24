#!/bin/bash/
#####################################################################################
#Script for running the MANTA software package. Intended for structural
#variant calling. MANTA is aprogram for the detection of structural variants
#from paired-end sequence data. Because running all individuals together is taxing on 
#MANTA, this program was run in 9 batches for male kakapo. 
#####################################################################################
manta=/scale_wlg_nobackup/filesets/nobackup/uoo02695/software/manta-1.6.0.centos6_x86_64/
input=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/alignments_male/bam/sorted/
jane=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/alignments_female/jane/jane_sim_noW.sorted.bam
output=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/manta/male/batch01/
ref=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/References/kakapo_no_Wchromosome.fa

${manta}/bin/configManta.py \
--bam ${jane} \
--bam ${input}Al.sorted.bam \
--bam ${input}Arab.sorted.bam \
--bam ${input}Ariki.sorted.bam \
--bam ${input}Attenborough.sorted.bam \
--bam ${input}Awhero.sorted.bam \
--bam ${input}Barnard.sorted.bam \
--bam ${input}Basil.sorted.bam \
--bam ${input}Ben.sorted.bam \
--bam ${input}Bill.sorted.bam \
--bam ${input}Blades.sorted.bam \
--referenceFasta ${ref} \
--runDir ${output}

##########################################################################################################################################
#Once the configuration directories and files are created, time to run the workflow pipeline.
##########################################################################################################################################

/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/manta/male/runWorkflow.py