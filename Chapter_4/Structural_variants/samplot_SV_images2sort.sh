#!/bin/sh
########################################################################################################
# From the samplot GitHub: samplot is a command line tool for rapid, multi-sample structural variant
# visualization. samplot takes SV coordinates and bam files and produces high-quality images that
# highlight any alignment and depth signals that substantiate the SV.
########################################################################################################
bin=/scale_wlg_nobackup/filesets/nobackup/uoo02695/software/samplot/bin/
output=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/samplot_img
male=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/backup_bams/alignments_male
female=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/backup_bams/alignments_female


# Plotting SV's unique to RH lineage...
while read -r line
do
chrom=$(echo ${line} | awk '{print $1}')
starts=$(echo ${line} | awk '{print $2}')
end=$(echo ${line} | awk '{print $4}')
type=$(echo ${line} | awk '{print $3}')
echo "Plotting ${type} between ${starts} and ${end} on scaffold ${chrom}..."
time ${bin}samplot plot \
    -n Richard_Henry Flossie Kuia Gulliver Sinbad Awarua Elliot Felix Ian Kuihi Roha Sandra Scratch Suzanne Tutoko Waa\
    -b ${male}/Richard_Henry.sorted.bam \
    ${female}/Flossie.sorted.bam \
    ${female}/Kuia.sorted.bam \
    ${male}/Gulliver.sorted.bam \
    ${male}/Sinbad.sorted.bam \
    ${female}/Awarua.sorted.bam \
    ${male}/Elliot.sorted.bam \
    ${male}/Felix.sorted.bam \
    ${male}/Ian.sorted.bam \
    ${female}/Kuihi.sorted.bam \
    ${female}/Roha.sorted.bam \
    ${female}/Sandra.sorted.bam \
    ${male}/Scratch.sorted.bam \
    ${female}/Suzanne.sorted.bam \
    ${male}/Tutoko.sorted.bam \
    ${female}/Waa.sorted.bam \
    -o ${output}/2sort/${chrom}_${starts}_${end}_${type}.png \
    -c ${chrom} \
    -s ${starts} \
    -e ${end} \
    -t ${type}
done < ${output}/shared_samplot.tsv# uniq2SI_samplot.tsv uniq2RH_samplot.tsv