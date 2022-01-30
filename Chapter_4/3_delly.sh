#!/bin/bash -e
##################################################################################################################################
# Running the Delly SV discovery and genotyping pipeline. Delly is a programme for structural variant discovery with paired-end
# sequence data. Delly requires sorted bam files with marked duplicates (as in the script 2_align_stat.sh) as input. 
##################################################################################################################################
male=/kakapo-data/bwa/bwa_male/markdup/
female=/kakapo-data/bwa/bwa_female/markdup/
out=/kakapo-data/bwa/delly/markdup/
ref=/kakapo-data/References/kakapo_full_ref.fa
exclude=/kakapo-data/metadata/kakapo_SVexcluded_scaffolds.bed

printf "\nRunning Delly...\n"
for i in {01..14}
       do
       for mbam in ${male}batch${i}/*_markdup.bam
               do
               mbase=$(basename ${mbam} _markdup.bam)
               printf "\nRunning Delly for ${mbase}...\n"
               delly call -g ${ref} -x ${exclude} -o ${out}SV_calls_male/${mbase}.bcf ${mbam}
       done &
done
for j in {01..14}
       do
       for fbam in ${female}batch${j}/*_markdup.bam
               do
               fbase=$(basename ${fbam} _markdup.bam)
               printf "\nRunning Delly for ${fbase}...\n"
               delly call -g ${ref} -x ${exclude} -o ${out}SV_calls_female/${fbase}.bcf ${fbam}
       done &
done
wait
printf "\nMerging SV sites...\n"
delly merge -o ${out}kakapo_bwa_delly_sites.bcf ${out}SV_calls_male/*.bcf ${out}SV_calls_female/*.bcf
delly merge -o ${out}kakapo_bwa_delly_minSize300bp_sites.bcf -m 300 ${out}SV_calls_male/*.bcf ${out}SV_calls_female/*.bcf

printf "\nRunning Delly genotyping...\n"
for i in {01..14}
       do
       for mbam in ${male}batch${i}/*_markdup.bam
               do
               mbase=$(basename ${mbam} _markdup.bam)
               printf "\nRunning Delly genotyping for ${mbase}..."
               delly call -g ${ref} -v ${out}bwa_delly_sites.bcf -o ${out}genotype/${mbase}.geno.bcf ${mbam}
               delly call -g ${ref} -v ${out}bwa_delly_minSize300bp_sites.bcf -o ${out}genotype/${mbase}.geno.minsize.bcf ${mbam}
       done &
done
for j in {01..11}
        do
        for fbam in ${female}batch${j}/*_markdup.bam
                do
                fbase=$(basename ${fbam} _markdup.bam)
                printf "\nRunning Delly genotyping for ${fbase}...\n"
                delly call -g ${ref} -v ${out}bwa_delly_minSize300bp_sites.bcf -o ${out}genotype/${fbase}.geno.minsize.bcf ${fbam}
               delly call -g ${ref} -v ${out}bwa_delly_sites.bcf -o ${out}genotype/${fbase}.geno.bcf ${fbam}
        done &
done
wait

bcftools merge -m id -O b -o ${out}01_bwa_delly_genotypes.bcf --threads 24 ${out}genotype/*geno.bcf # Used for unfiltered SV stats 
#bcftools merge -m id -O b -o ${out}bwa_delly_minsize_genotypes.bcf --threads 24 ${out}genotype/*.minsize.bcf
tabix ${out}01_bwa_delly_genotypes.bcf
#tabix ${out}bwa_delly_minsize_genotypes.bcf

delly filter -f germline -p -m 50 -o ${out}02_bwa_delly_germline_minimum50bp.bcf ${out}01_bwa_delly_genotypes.bcf 
delly filter -f germline -p -m 300 -o ${out}03_bwa_delly_germline_minimum300bp.bcf ${out}01_bwa_delly_genotypes.bcf

bcftools view -t ^NC_044302.2 -i '(SVTYPE = "INS" & FILTER == "PASS")' \
    -O b -o ${out}04_bwa_delly_ins.bcf ${out}01_bwa_delly_genotypes.bcf # Used for filtered INS stats
bcftools view -t ^NC_044302.2 -i '(SVTYPE = "DEL")' \
    -O b -o ${out}05_bwa_delly_del.bcf 02_bwa_delly_germline_minimum50bp.bcf
bcftools view -t ^NC_044302.2 -i '(SVTYPE = "INV") | (SVTYPE = "DUP")' \
    -O b -o ${out}06_bwa_delly_inv_dup.bcf ${out}03_bwa_delly_germline_minimum300bp.bcf
tabix ${out}04_bwa_delly_ins.bcf
tabix ${out}05_bwa_delly_del.bcf
tabix ${out}06_bwa_delly_inv_dup.bcf

bcftools concat -a -O v -o 07_delly_SVfilter.vcf ${out}04_bwa_delly_ins.bcf ${out}05_bwa_delly_del.bcf ${out}06_bwa_delly_dup_inv.bcf
bcftools view -i '(N_PASS(FT!="PASS") / N_PASS(GT!="RR")) >=  0.8' -O v -o 08_delly_genofilter.vcf 07_delly_SVfilter.vcf

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\tdelly_unfiltered\n' ${out}01_bwa_delly_genotypes.bcf > ${out}delly_summary.tsv
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\tdelly_SVfiltered\n' ${out}07_delly_SVfilter.vcf >> ${out}delly_summary.tsv
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\tdelly_genofiltered\n' ${out}08_delly_genofilter.vcf >> ${out}delly_summary.tsv
##################################################################################################################################
# Now for trio filters. 
##################################################################################################################################
bcftools +mendelian -m a -T ${trio} -O v -o ${out}09_delly_SVfilter_trio.vcf \
    ${out}07_delly_SVfilter.vcf
bcftools +mendelian -m a -T ${trio} -O v -o ${out}10_delly_genofilter_trio.vcf \
    ${out}08_delly_genofilter.vcf

bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0' -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\t0_fail_delly_SVfilter\n' \
    ${out}09_delly_SVfilter_trio.vcf >> ${out}delly_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.05' -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\t<=0.05_fail_delly_SVfilter\n' \
    ${out}09_delly_SVfilter_trio.vcf >> ${out}delly_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.1' -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\t<=0.1_fail_delly_SVfilter\n' \
    ${out}09_delly_SVfilter_trio.vcf >> ${out}delly_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.2' -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\t<=0.2_fail_delly_SVfilter\n' \
    ${out}09_delly_SVfilter_trio.vcf >> ${out}delly_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0' -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\t0_fail_delly_genofilter\n' \
    ${out}10_delly_genofilter_trio.vcf >> ${out}delly_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.05' -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\t<=0.05_fail_delly_genofilter\n' \
    ${out}10_delly_genofilter_trio.vcf >> ${out}delly_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.1' -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\t<=0.1_fail_delly_genofilter\n' \
    ${out}10_delly_genofilter_trio.vcf >> ${out}delly_mendel.tsv
bcftools query -i '(MERR / N_PASS(GT!="RR")) <=0.2' -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\t<=0.2_fail_delly_genofilter\n' \
    ${out}10_delly_genofilter_trio.vcf >> ${out}delly_mendel.tsv

##################################################################################################################################
# Lineage Comparisons
##################################################################################################################################
mkdir -p ${out}lineage_comparisons

bcftools view -s Richard_Henry ${out}10_delly_genofilter_trio.vcf | bcftools view -i 'GT!="RR" & GT!="mis"' -O v -o ${out}lineage_comparisons/RH_variants.vcf
bcftools view -s ^Richard_Henry,Kuia,Gulliver,Sinbad,Adelaide,Henry,Marian,Gertrude ${out}10_delly_genofilter_trio.vcf | bcftools view -i 'GT!="RR" & GT!="mis"' -O v -o ${out}lineage_comparisons/SI_variants.vcf

tabix ${out}lineage_comparisons/RH_variants.vcf
tabix ${out}lineage_comparisons/SI_variants.vcf

bcftools isec ${out}lineage_comparisons/RH_variants.bcf \
    ${out}lineage_comparisons/SI_variants.bcf \
    -p ${out}lineage_comparisons/

# Summarising numbers of SVs per individual
bcftools view -h ${out}08_bwa_delly_final_trio.vcf | grep CHROM | tr "\t" "\n" | tail -n 169 > ${out}samples.txt
bcftools query -f '%CHROM\t%POS\n' ${out}lineage_comparisons/0000.vcf > ${out}lineage_comparisons/Fiordland_sites.txt
bcftools query -f '%CHROM\t%POS\n' ${out}lineage_comparisons/0001.vcf > ${out}lineage_comparisons/SI_sites.txt
bcftools query -f '%CHROM\t%POS\n' ${out}lineage_comparisons/0002.vcf > ${out}lineage_comparisons/shared_sites.txt

while read -r line
    do
    echo "Counting SVs for ${line}..."
    si=$(bcftools view -T ${out}lineage_comparisons/SI_sites.txt -s ${line} ${out}08_delly_final_trio.vcf | bcftools query -i 'GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    sh=$(bcftools view -T ${out}lineage_comparisons/shared_sites.txt -s ${line} ${out}08_delly_final_trio.vcf | bcftools query -i 'GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    rh=$(bcftools view -T ${out}lineage_comparisons/Fiordland_sites.txt -s ${line} ${out}08_delly_final_trio.vcf | bcftools query -i 'GT!="RR" & GT!="mis"' -f '%SVTYPE\n' | wc -l)
    printf "${line}\t${si}\tdelly_SI\n" >> ${out}lineage_comparisons/delly_indiv_SVcounts.tsv
    printf "${line}\t${sh}\tdelly_shared\n" >> ${out}lineage_comparisons/delly_indiv_SVcounts.tsv
    printf "${line}\t${rh}\tdelly_RH\n" >> ${out}lineage_comparisons/delly_indiv_SVcounts.tsv
    bcftools view -T ${out}lineage_comparisons/SI_sites.txt -s ${line} ${out}08_delly_final_trio.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tdelly_SI\n' >> detailed_delly_summary.tsv
    bcftools view -T ${out}lineage_comparisons/shared_sites.txt -s ${line} ${out}08_delly_final_trio.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tdelly_shared\n' >> detailed_delly_summary.tsv
    bcftools view -T ${out}lineage_comparisons/Fiordland_sites.txt -s ${line} ${out}08_delly_final_trio.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\tdelly_RH\n' >> detailed_delly_summary.tsv
done < ${out}samples.txt


while read -r line
    do
    indiv=$(echo $line | awk '{print $1}')
    gen=$(echo $line | awk '{print $2}')
    echo "Counting SVs for $indiv..."
    bcftools view -s ${indiv} ${out}10_delly_genofilter_trio.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVTYPE\n' >> ${out}delly_indiv_counts.tsv
    echo "Assigning $gen generation to $indiv"
    grep "^$indiv" ${out}delly_indiv_counts.tsv | awk -v var="$gen" '{print $0"\t"var}' >> ${out}delly_generations.tsv
done < /kakapo-data/metadata/generations.tsv

##################################################################################################################################
# Preparing data for R markdown
##################################################################################################################################
# Converting SV types into long-form names and updating NCBI scaffold names to chromosome
while read -r line
    do
    ncbi=$(echo $line | awk '{print $1}')
    chr=$(echo $line | awk '{print $2}')
    echo Converting $ncbi to $chr
    sed -i "s/$ncbi/$chr/g" ${out}delly_summary.tsv
done < ${out}convert_chr.txt 
# Adding chromosome size to summary file
while read -r line
    do
    chr=$(echo $line | awk '{print $1}')
    size=$(echo $line | awk '{print $2}')
    echo "Adding $size for $chr..."
    grep "^${chr}" ${out}delly_summary.tsv | awk -v var=$size '{print $0"\t"var}' >> ${out}delly_size
done < ${out}chrom_sizes 
mv ${out}delly_size ${out}delly_summary.tsv
