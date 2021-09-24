#!/bin/sh
###################################################################################
#Here are some merging strategies I tested to identify structural variants shared
#across batches in MANTA.
###################################################################################

#Creating fixed variables for directories.
Minput=/manta/male
Moutput=/manta/male
Finput=/manta/female
Foutput=/manta/female

#In this attempt, I tried to create VCF's for files that contained only sites shared across
#batches. 

for i in {1..9}
    do
    echo "Creating intersection of sites for batch ${i}..."
    bcftools isec -O v -o ${Moutput}/males_isec${i}.vcf -n +9 -w ${i} \
    ${Minput}/batch01_manta_INV_conversion_male.vcf.gz \
    ${Minput}/batch02_manta_INV_conversion_male.vcf.gz \
    ${Minput}/batch03_manta_INV_conversion_male.vcf.gz \
    ${Minput}/batch04_manta_INV_conversion_male.vcf.gz \
    ${Minput}/batch05_manta_INV_conversion_male.vcf.gz \
    ${Minput}/batch06_manta_INV_conversion_male.vcf.gz \
    ${Minput}/batch07_manta_INV_conversion_male.vcf.gz \
    ${Minput}/batch08_manta_INV_conversion_male.vcf.gz \
    ${Minput}/batch09_manta_INV_conversion_male.vcf.gz
done

#As you can see, this doens't exactly take care of merging all the vcfs into one
#file. But at least only the sites shared in all batches are called for each output
#file. To address the issue of having multiple batch files, I attempted to merge
#with...

bcftools merge -O v -o ${Moutput}/manta_males.vcf \
${Minput}/males_isec1.vcf.gz \
${Minput}/males_isec2.vcf.gz \
${Minput}/males_isec3.vcf.gz \
${Minput}/males_isec4.vcf.gz \
${Minput}/males_isec5.vcf.gz \
${Minput}/males_isec6.vcf.gz \
${Minput}/males_isec7.vcf.gz \
${Minput}/males_isec8.vcf.gz \
${Minput}/males_isec9.vcf.gz

#And repeated for the females. There were less batches of females, so less files to
#find the intersection of for merging.

for i in {1..7}
    do
    echo "Creating intersection of sites for batch ${i}..."
    bcftools isec -O v -o ${Foutput}/females_isec${i}.vcf -n +7 -w ${i} \
    ${Finput}/batch01_manta_INV_conversion_female.vcf.gz \
    ${Finput}/batch02_manta_INV_conversion_female.vcf.gz \
    ${Finput}/batch03_manta_INV_conversion_female.vcf.gz \
    ${Finput}/batch04_manta_INV_conversion_female.vcf.gz \
    ${Finput}/batch05_manta_INV_conversion_female.vcf.gz \
    ${Finput}/batch06_manta_INV_conversion_female.vcf.gz \
    ${Finput}/batch07_manta_INV_conversion_female.vcf.gz
done

bcftools merge -O v -o ${Foutput}/females.vcf \
${Finput}/females_isec1.vcf.gz \
${Finput}/females_isec2.vcf.gz \
${Finput}/females_isec3.vcf.gz \
${Finput}/females_isec4.vcf.gz \
${Finput}/females_isec5.vcf.gz \
${Finput}/females_isec6.vcf.gz \
${Finput}/females_isec7.vcf.gz \
${Finput}/females_isec8.vcf.gz \
${Finput}/females_isec9.vcf.gz

