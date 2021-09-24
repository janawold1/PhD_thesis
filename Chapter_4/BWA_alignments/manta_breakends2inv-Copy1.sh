#!/bin/bash

convert=/scale_wlg_nobackup/filesets/nobackup/uoo02695/software/manta-1.6.0.centos6_x86_64/libexec/convertInversion.py
sam=/scale_wlg_nobackup/filesets/nobackup/uoo02695/software/manta-1.6.0.centos6_x86_64/libexec/samtools
fref=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/References/kakapo_full_ref.fa
mref=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/References/kakapo_no_Wchromosome.fa
finput=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/manta/joint_calling/female_diploidSV.vcf
minput=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/manta/joint_calling/male_diploidSV.vcf

${convert} ${sam} ${fref} ${finput} > female_inv_diploidSV.vcf

${convert} ${sam} ${mref} ${minput} > male_inv_diploidSV.vcf