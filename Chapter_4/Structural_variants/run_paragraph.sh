##############################################################################################################
# There is a bit of data wrangling before you can start with genotyping with Paragraph. To start, you need
# to generate a sample manifest. This is a .tsv that has the following required fields: UniqueID, path/to/bam,
# average depth across the genome, and average read length across the genome. 
##############################################################################################################

# First step was to calculate the average read depth across the genome in males where the W-chr was exluded...
for mbam in alignments_male/*.sorted.bam
    do
    base=$(basename ${mbam} .sorted.bam)
    echo "Average depth for ${base}..."
    samtools depth -a ${mbam} | awk '{sum+=$3} END { print "Average = ",sum/1111171819}'
done >> paragraph/male_depth.txt

# And in females... The size of the female genome is different since it includes the W-chr.
for fbam in alignments_female/*.sorted.bam
    do
    basef=$(basename ${fbam} .sorted.bam)
    echo "Average depth for ${basef}..."
    samtools depth -a ${fbam} | awk '{sum+=$3} END { print "Average = ",sum/1148597545}'
done > paragraph/female_depth.txt 

# Once the average depth across the genome is calculated, paragraph recommends creating a unique individual 
# sample manifest for every sample you want to genotype in a population. I did this as per below...
cat paragraph/female_depth.txt \
    | tr "\n" "," \
    | sed 's/Average = //g' \
    | sed 's/Average depth for //g' \
    | sed -e 's/..., /\t/g' \
    | tr "," "\n" \
    | awk '{print $1",/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/alignments_female/"$1".sorted.bam,"$2",125"}' \
    > paragraph/female_manifest.csv 
    
# To create individual files for the sample manifests. Although paragraph says a tsv will work, I could only get it to run with csv...
# Each column had to be labelled in each manifest as "id,path,deph,read length"
while read -r line
    do
    indiv=$(echo ${line} | tr "," "\t" | awk '{print $1}')
    deets=$(echo ${line})
    echo creating manifest for $indiv
    echo "id,path,depth,read length" >> female_manifests/${indiv}.csv
    echo ${deets} >> female_manifests/${indiv}.csv
done < female_manifest.csv

# This was repeated for males. 
# Paragraph cannot genotype insertions or breakpoint/translocations. For this, both INS 
# and BND were removed from the Manta total SV call set as per below. This took the 
# number of SV's for genotyping from 145,705 to 71,637.
bcftools view \
    -e 'SVTYPE="BND" | SVTYPE="INS"' \
    -G -O v -o female_DEL_DUP_INV.vcf \
    female_inv_diploidSV.vcf.gz
    
bcftools filter -e 'SVTYPE="BND"' -o joint_noINS_noBND.vcf.gz -O z joint_noINS.vcf.gz
# 'Bad' sites that reference a BND were identified with...
zcat joint_noINS_noBND.vcf.gz | grep -v "#" | awk '{print $5}' | grep ","
# Locations where INV related to a BND to output to sites2rmv.tsv
bcftools view -O z -o joint_noINS_noBND_nobadsites.vcf.gz -T ^ sites2rmv.tsv joint_noINS_noBND.vcf.gz

# Then, paragraph was then run as below...
module load Anaconda3/2019.03-gimkl-2018b
conda activate paragraph

candidate_svs=joint_noINS_noBND_nobadsites.vcf.gz
ref=/kakapo-data/References/kakapo_full_ref.fa
output=/kakapo-data/paragraph/female_genos

for female in paragraph/female_manifests/*.csv
    do
    base=$(basename ${female} .csv)
    echo "Genotyping ${base}..."
    mkdir ${output}/${base}
    multigrmpy.py -i ${candidate_svs} \
        -m ${female} \
        -r ${ref} \
        -o ${output}/${base}
done