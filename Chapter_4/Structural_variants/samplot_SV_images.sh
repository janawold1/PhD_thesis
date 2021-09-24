#!/bin/sh
########################################################################################################
# From the samplot GitHub: samplot is a command line tool for rapid, multi-sample structural variant
# visualization. samplot takes SV coordinates and bam files and produces high-quality images that
# highlight any alignment and depth signals that substantiate the SV.
########################################################################################################
bin=/scale_wlg_nobackup/filesets/nobackup/uoo02695/software/samplot/bin/
output=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/samplot_img/joint_calls
male=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/alignments_male/bam/sorted
female=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/alignments_female/bam/sorted


# Plotting SV's unique to RH lineage...
while read -r line
do
chrom=$(echo ${line} | awk '{print $1}')
starts=$(echo ${line} | awk '{print $2}')
end=$(echo ${line} | awk '{print $4}')
type=$(echo ${line} | awk '{print $3}')
echo "Plotting ${type} between ${starts} and ${end} on scaffold ${chrom} for the first generation of RH lineage..."
time ${bin}samplot plot \
    -n Richard_Henry Flossie Kuia Gulliver Sinbad \
    -b ${male}/Richard_Henry.sorted.bam \
    ${female}/Flossie.sorted.bam \
    ${female}/Kuia.sorted.bam \
    ${male}/Gulliver.sorted.bam \
    ${male}/Sinbad.sorted.bam \
    -o ${output}/Richard_H_lineage/F1_${chrom}_${starts}_${end}_${type}.png \
    -c ${chrom} \
    -s ${starts} \
    -e ${end} \
    -t ${type}
done < uniq2RH_samplot.tsv

# Plotting SV's unique to everyone else...
while read -r line
do
chrom=$(echo ${line} | awk '{print $1}')
starts=$(echo ${line} | awk '{print $2}')
end=$(echo ${line} | awk '{print $4}')
type=$(echo ${line} | awk '{print $3}')
echo "Plotting ${type} between ${starts} and ${end} on scaffold ${chrom} for individuals with no RH lineage..."
time ${bin}samplot plot \
    -n Whiskas Flossie Scratch Wiremu \
    -b ${male}/Whiskas.sorted.bam \
    ${female}/Flossie.sorted.bam \
    ${male}/Scratch.sorted.bam \
    ${male}/Wiremu.sorted.bam \
    -o ${output}/noSI_lineage/${chrom}_${starts}_${end}_${type}.png \
    -c ${chrom} \
    -s ${starts} \
    -e ${end} \
    -t ${type}
done < uniq2noSI_samplot.tsv

#Plotting SV's shared between groups...
while read -r line
do
chrom=$(echo ${line} | awk '{print $1}')
starts=$(echo ${line} | awk '{print $2}')
end=$(echo ${line} | awk '{print $4}')
type=$(echo ${line} | awk '{print $3}')
echo "Plotting ${type} between ${starts} and ${end} on scaffold ${chrom} shared between groups..."
time ${bin}samplot plot \
    -n Richard_Henry Flossie Kuia Gulliver Sinbad Whiskas Wiremu Scratch \
    -b ${male}/Richard_Henry.sorted.bam \
    ${female}/Flossie.sorted.bam \
    ${female}/Kuia.sorted.bam \
    ${male}/Gulliver.sorted.bam \
    ${male}/Sinbad.sorted.bam \
    ${male}/Whiskas.sorted.bam \
    ${male}/Wiremu.sorted.bam \
    ${male}/Scratch.sorted.bam \
    -o ${output}/shared/F1_${chrom}_${starts}_${end}_${type}.png \
    -c ${chrom} \
    -s ${starts} \
    -e ${end} \
    -t ${type}
done < shared_samplot.tsv