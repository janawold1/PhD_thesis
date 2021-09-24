##############################################################################################################
# There is a bit of data wrangling before you can start with genotyping with Paragraph. To start, you need
# to generate a sample manifest. This is a .tsv that has the following required fields: UniqueID, path/to/bam,
# average depth across the genome, and average read length across the genome. 
##############################################################################################################

# First step was to calculate the average read depth across the genome in males where the W-chr was exluded...
for mbam in alignments_male/*.sorted.bam
    do
    base=$(basename ${mbam} .sorted.bam)
    echo "Calculating depth for ${base}..."
    samtools depth -a ${mbam} | awk '{sum+=$3} END { print "Average = ",sum/1111171819}'
done > paragraph/male_depth.txt

# And in females... The size of the female genome is different since it includes the W-chr.
for fbam in alignments_female/*.sorted.bam
    do
    basef=$(basename ${fbam} .sorted.bam)
    echo "Average depth for ${basef}..."
    samtools depth -a ${fbam} | awk '{sum+=$3} END { print "Average = ",sum/1148597545}'
done > paragraph/female_depth.txt 

# Cannot merge raw  Manta files for males and females before separating by SV type.
# This is because 2 different SV types are merged if they overlap. To avoid this, I split each raw male and female
# file Inversion BND calls have been converted to Inversion
# symbolic alleles.
bcftools merge \
    --threads 32 \
    --merge none \
    --force-samples \
    -o kakapo_joint.vcf \
    -O v \
    female_inv_diploidSV.vcf.gz \
    male_inv_diploidSV.vcf.gz 

# Count number of each type of SV in merged file.
zcat kakapo_joint_inv_diploid.vcf.gz |\
    grep -v "#" |\
    awk '{print $3}' |\
    awk 'BEGIN {FS = ":"}; {print $1}' |\
    sort |\
    uniq -c

# Which tells us there are:
#  72272 MantaBND
#   3911 MantaDEL
#   1439 MantaDUP
#   1849 MantaINS
#  66298 MantaINV

##############################################################################################################
# Because I have had a lot of issues genotyping multiple SV types at once in the past, I'm going to attempt
# genotyping the each SV type on its own. 
##############################################################################################################
bcftools view \
    -i 'SVTYPE == "DEL" && FILTER == "PASS"' \
    -G \
    -O v \
    -o joint_male_DEL.vcf \
    male_inv_diploid.vcf.gz

bcftools view \
     -i 'SVTYPE == "DUP" && FILTER == "PASS"' \
     -G \
     -O z \
     -o joint_male_DUP.vcf.gz \
     male_inv_diploidSV.vcf.gz

bcftools view \
     -i 'SVTYPE == "INS" && FILTER == "PASS"' \
     -G \
     -O z \
     -o joint_male_INS.vcf.gz \
     male_inv_diploidSV.vcf.gz

zcat input_vcfs/joint_male_INS.vcf.gz \
    | grep "#" >> input_vcfs/joint_male_INS_filtered.vcf
zcat input_vcfs/joint_male_INS.vcf.gz \
    | grep -v "#" \
    | grep -v "RIGHT_SVINSSEQ" \
    | grep -v "LEFT_SVINSSEQ" >> input_vcfs/joint_male_INS_filtered.vcf

bcftools view \
     -i 'SVTYPE == "INV" && FILTER == "PASS"' \
     -G \
     -O z \
     -o joint_male_INV.vcf -G \
     male_inv_diploidSV.vcf.gz
grep "#" joint_male_INV.vcf >> joint_male_INV_filtered.vcf
grep -v "#" joint_male_INV.vcf | grep CIPOS | grep CIEND >> joint_male_INV_filtered.vcf

#Filtered INV by size with:
zcat input_vcfs/joint_male_INV_filtered.vcf.gz \
    | grep -v "#" \
    | awk '{print $1"\t"$2"\t"$8}' \
    | tr ";" "\t" \
    | awk '{print $1"\t"$2"\t"$5}' \
    | sed 's/SVLEN=//g' \
    | awk '$3<100' \
    | awk '{print $1"\t"$2}' > INV_2_rmv.tsv

bcftools view \
    -T ^INV_2_rmv.tsv \
    -O v \
    -o input_vcfs/joint_male_INV_filtered_by_size.vcf.gz \
    input_vcfs/joint_male_INV_filtered.vcf.gz

##############################################################################################################
# Repeat for the ladies
##############################################################################################################

bcftools view \
    -i 'SVTYPE == "DEL" && FILTER == "PASS"' \
    -t ^NC_044301.2 \
    -G \
    -O v \
    -o joint_female_DEL.vcf \
    female_inv_diploid.vcf.gz

bcftools view \
     -i 'SVTYPE == "DUP" && FILTER == "PASS"' \
     -t ^NC_044301.2 \
     -G \
     -O z \
     -o input_vcfs/joint_female_DUP.vcf.gz \
     female_inv_diploidSV.vcf.gz

bcftools view \
     -i 'SVTYPE == "INS" && FILTER == "PASS"' \
     -t ^NC_044301.2 \
     -G \
     -O v \
     -o joint_female_INS.vcf \
     female_inv_diploidSV.vcf.gz
     
grep "#" joint_female_INS.vcf >> input_vcfs/joint_female_INS_filtered.vcf
grep -v "#" joint_female_INS.vcf \
    | grep -v "RIGHT_SVINSSEQ" \
    | grep -v "LEFT_SVINSSEQ" >> input_vcfs/joint_female_INS_filtered.vcf

rm joint_female_INS.vcf

# Tried running Paragraph on these insertions, but it will likely fail.
# Identified which INS to remove, edited the output into the correct format,
# then removed these sites with:
# grep "Exception: Missing key SEQ for <INS> at " paragraph_male_ins_female_geno.out \
#    | awk 'BEGIN {FS = " "}; {print $11}' \
#    | tr ":" "\t" \
#    | sed 's/;//g' \
#    | sort > ins2rmv.tsv
# But didn't need to filter insertions in females as done in males. 

bcftools view \
     -i 'SVTYPE == "INV" && FILTER == "PASS"' \
     -t ^NC_044301.2 \
     -G \
     -O v \
     -o joint_female_INV.vcf \
     female_inv_diploidSV.vcf.gz
     
grep "#" joint_female_INV.vcf >> joint_female_INV_filtered.vcf
grep -v "#" joint_female_INV.vcf | grep CIPOS | grep CIEND >> joint_female_INV_filtered.vcf

#Filtered INV by size with:
grep -v "#" joint_female_INV_filtered.vcf \
    | awk '{print $1"\t"$2"\t"$8}' \
    | tr ";" "\t" \
    | awk '{print $1"\t"$2"\t"$5}' \
    | sed 's/SVLEN=//g' \
    | awk '$3<100' \
    | awk '{print $1"\t"$2}' > INV_2_rmv.tsv

bcftools view \
    -T ^INV_2_rmv.tsv \
    -O z \
    -o input_vcfs/joint_female_INV_filtered_by_size.vcf.gz \
    joint_female_INV_filtered.vcf

#Stitch back together into a single vcf
for file in input_vcfs/*.vcf.gz
do
    tabix ${file}
done

bcftools merge \
    -m none \
    -O v \
    -o input_vcfs/kakapo_paragraph_variants.vcf \
    input_vcfs/joint_male_DEL.vcf.gz \
    input_vcfs/joint_male_DUP.vcf.gz \
    input_vcfs/joint_male_INS_filtered.vcf.gz \
    input_vcfs/joint_male_INV_filtered_by_size.vcf.gz \
    input_vcfs/joint_female_DEL.vcf.gz \
    input_vcfs/joint_female_DUP.vcf.gz \
    input_vcfs/joint_female_INS_filtered.vcf.gz \
    input_vcfs/joint_female_INV_filtered_by_size.vcf.gz

work_dir=/kakapo-data/paragraph
ref_dir=/kakapo-data/References
in_VCF=/kakapo-data/paragraph/manta_joint/input_vcfs/kakapo_paragraph_variants.vcf

nohup multigrmpy.py \
    -i ${in_VCF} \
    -m ${work_dir}/female_manifest.csv \
    -o ${work_dir}/manta_joint/genotypes/females \
    -r ${ref_dir}/kakapo_full_ref.fa \
    --threads 32 \
    --scratch-dir /kakapo-data/paragraph/tmp_json/ \
    --verbose >> paragraph_female_geno.out &
disown

nohup multigrmpy.py \
    -i ${in_VCF} \
    -m ${work_dir}/male_manifest.csv \
    -o ${work_dir}/manta_joint/genotypes/males \
    -r ${ref_dir}/kakapo_no_Wchromosome.fa \
    --threads 32 \
    --scratch-dir /kakapo-data/paragraph/tmp_json/ \
    --verbose >> paragraph_male_geno.out &
disown

##############################################################################################################
# Running paragraph on batched SV calls
##############################################################################################################

bcftools merge \
    --merge none \
    --force-samples \
    -O v \
    -o /kakapo-data/paragraph/manta_batch/input_vcfs/female_batch.vcf \
    /kakapo-data/paragraph/manta_batch/manta_vcfs/*female.vcf.gz

bcftools view \
    -i 'SVTYPE == "DEL" && FILTER == "PASS"' \
    -t ^NC_044301.2 \
    -G \
    -O z \
    -o batch_female_DEL.vcf.gz \
    female_batch.vcf 

bcftools view \
    -i 'SVTYPE == "DUP" && FILTER == "PASS"' \
    -t ^NC_044301.2 \
    -G \
    -O z \
    -o batch_female_DUP.vcf.gz \
    female_batch.vcf 

bcftools view \
    -i 'SVTYPE == "INS" && FILTER == "PASS"' \
    -t ^NC_044301.2 \
    -G \
    -O z \
    -o batch_female_INS.vcf.gz \
    female_batch.vcf

grep "#" batch_female_INS.vcf >> input_vcfs/batch_female_INS_filtered.vcf
grep -v "#" batch_female_INS.vcf \
    | grep -v "RIGHT_SVINSSEQ" \
    | grep -v "LEFT_SVINSSEQ" >> input_vcfs/batch_female_INS_filtered.vcf
    
bcftools view \
    -i 'SVTYPE == "INV" && FILTER == "PASS"' \
    -t ^NC_044301.2 \
    -G \
    -O z \
    -o batch_female_INV.vcf.gz \
    female_batch.vcf

grep "#" batch_female_INV.vcf >> batch_female_INV_filtered.vcf
grep -v "#" batch_female_INV.vcf | grep CIPOS | grep CIEND >> batch_female_INV_filtered.vcf

#Filtered INV by size with:
grep -v "#" batch_female_INV_filtered.vcf \
    | awk '{print $1"\t"$2"\t"$8}' \
    | tr ";" "\t" \
    | awk '{print $1"\t"$2"\t"$5}' \
    | sed 's/SVLEN=//g' \
    | awk '$3<100' \
    | awk '{print $1"\t"$2}' > INV_2_rmv.tsv

################################################################################################
bcftools merge \
    --merge none \
    --force-samples \
    -O v \
    -l /kakapo-data/paragraph/manta_batch/input_vcfs/male_batch.vcf \
    /kakapo-data/paragraph/manta_batch/manta_vcfs/*male.vcf.gz

bcftools view \
    -i 'SVTYPE == "DEL" && FILTER == "PASS"' \
    -G \
    -O z \
    -o batch_male_DEL.vcf.gz \
    male_batch.vcf 

bcftools view \
    -i 'SVTYPE == "DUP" && FILTER == "PASS"' \
    -G \
    -O z \
    -o batch_male_DUP.vcf.gz \
    male_batch.vcf 

bcftools view \
    -i 'SVTYPE == "INS" && FILTER == "PASS"' \
    -G \
    -O z \
    -o batch_male_INS.vcf.gz \
    male_batch.vcf

grep "#" batch_female_INS.vcf >> input_vcfs/batch_female_INS_filtered.vcf
grep -v "#" batch_female_INS.vcf \
    | grep -v "RIGHT_SVINSSEQ" \
    | grep -v "LEFT_SVINSSEQ" >> input_vcfs/batch_female_INS_filtered.vcf

bcftools view \
    -i 'SVTYPE == "INV" && FILTER == "PASS"' \
    -t ^NC_044301.2 \
    -G \
    -O z \
    -o batch_female_INV.vcf.gz \
    female_batch.vcf

grep "#" batch_female_INV.vcf >> batch_female_INV_filtered.vcf
grep -v "#" batch_female_INV.vcf | grep CIPOS | grep CIEND >> batch_female_INV_filtered.vcf

#Filtered INV by size with:
grep -v "#" batch_female_INV_filtered.vcf \
    | awk '{print $1"\t"$2"\t"$8}' \
    | tr ";" "\t" \
    | awk '{print $1"\t"$2"\t"$5}' \
    | sed 's/SVLEN=//g' \
    | awk '$3<100' \
    | awk '{print $1"\t"$2}' > INV_2_rmv.tsv


##############################################################################################################
# 
##############################################################################################################

nohup multigrmpy.py -i male_joint_noBND_noGeno_filteredINS.vcf \
    -m ../female_manifest.csv \
    -o ./maleSVs/females \
    -r /kakapo-data/References/kakapo_no_Wchromosome.fa \
    --threads 32 \
    --scratch-dir /kakapo-data/paragraph/tmp_json/ \
    --verbose \
    > paragraph.out &
disown