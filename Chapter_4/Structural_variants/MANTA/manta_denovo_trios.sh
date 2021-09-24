#!/bin/bash
###################################################################################
# Here I am attempting to conduct the trio binning for the called SV's. The 
# sample_trios file was modified from the Pedigree_19_July_2020.filtered.txt file. 
# In all, there are 124 trios for comparison.
###################################################################################

# First began by creating fixed variables
manta=/home/rccuser/anaconda3/envs/manta/share/manta-1.6.0-0/libexec/denovo_scoring.py
work=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/manta/joint_calling/
trio=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/manta/joint_calling/trios/sample_trios
output=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/manta/joint_calling/trios
cd ${output}
module load BCFtools

# After converting breakends to inversion symbolic call, we need to merge all 
# discovered SVs for genotyping.Have removed the W to make genotyping easier.

bcftools merge -O z -o joint_kakapo_inv_diploidSV.vcf.gz \
    -m all -t ^NC_044301.2 \
    ${work}male_inv_diploidSV.vcf.gz \
    ${work}female_inv_diploidSV.vcf.gz

# The manta script automatically writes out to a new file. The name this file is 
# written to cannot be changed as far as I know. To address this, I'm creating vcf's
# for each individual trio comparison and running the de novo pipeline as per below.
while read -r line
    do
    fam=$(echo ${line} | sed 's/ /,/g')
    indiv=$(echo ${line} | awk '{print $1}')
    mantafam=$(echo ${line})
    echo "Creating trio VCF for ${line}"
    bcftools view -O v -o ${indiv}_trio.vcf -s ${fam} \
    ${work}joint_kakapo_inv_diploidSV.vcf
    echo "Testing trio for ${indiv}..."
    ${manta} ${indiv}_trio.vcf ${mantafam}
done < ${trio}

# Finally, I'm merging the stats generated from MANTA into one file.
for stat in *_trio.de_novo.stats.txt
    do
    base=$(basename ${stat} _trio.de_novo.stats.txt)
    echo ${base} >> all_stats.txt
    cat ${stat} >> all_stats.txt
done

# Now to see how many of these non-mendelian sites are shared in the population
for vcf in *_trio.de_novo.vcf
    do
    bgzip ${vcf}
    tabix ${vcf}.gz
done

# To get all the sites not adhering to mendelian inheritance, I first created
# a file with these individuals and sites.

for vcf in *_trio.de_novo.vcf.gz
    do
    indiv=$(basename ${vcf} _trio.de_novo.vcf.gz)
    echo "Creating VCF of failed sites for ${indiv}..."
    bcftools view -O z -o ${indiv}_noMendel.vcf.gz -i '(DQ)==60' -s ${indiv} ${vcf}
    echo "Now indexing VCF of failed sites for ${indiv}..."
    tabix ${indiv}_noMendel.vcf.gz
done

ls -lh *_noMendel.vcf.gz | awk '{print $9}' | xargs > 2merge

while read -r line
    do
    echo "Merging all failed sites..."
    bcftools merge -O v -o NoMendel.vcf -m all ${line}
done < 2merge

# Now to remove all failing sites from the main data set.
awk '{print $1, $2}' NoMendel.vcf | tr " " "\t" > sites2rmv.tsv

bcftools isec -O v -o Mendel.vcf -C ${work}joint_kakapo_inv_diploidSV.vcf \
    ${output}/NoMendel.vcf