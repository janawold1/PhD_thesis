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