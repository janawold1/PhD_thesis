#!/bin/sh
#####################################################################################
# This is a walkthrough of how I completed local assemblies in regions with putative 
# structural variants. 
#####################################################################################

#I tested 2 local assembly size ranges. The first one included overlapping reads 1kb
# up and down stream of the putative SV. The second included overlapping reads 5kb up 
# and down stream of the putative SV.

input=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/alignments_female/bam/sorted/
work=/scale_wlg_nobackup/filesets/nobackup/uoo02695/Kakapo/Sex_Chromosome/local_assemblies/
spades=/scale_wlg_nobackup/filesets/nobackup/uoo02695/sofware/SPAdes-3.14.1-Linux/bin/spades.py

#module purge
#module load SAMtools
#module load BEDTools

# To loop through the regions of interest, I first converted the bed file to this 
# format: "Chr:start-end"
#awk'{print $1, ":"$2, "-"$3}' ${5kb_start_end_mantaSV.bed} | sed 's/ //g' > ${5kb_start_end_mantaSV.loopinput}

# Then Isolated the reads with...
echo "Now beginning to isolate reads..."

# Before running this loop, it is imperative to make you are isolating read pairs
# properly. Spades cannot handle an uneven number of reads in the forward and reverse
# fastq files....
#mkdir -p ${work}/{5kbreads,5kbbams,5kbregions2assemble,5kb_intermediatebams,5kb_intermediatereads}

for bam in ${input}*.bam
    do
    base=$(basename ${bam} .sorted.bam)
    echo "Beginning local assemblies for ${base}..."
    while read -r line
        do
        filename=$(echo ${line} | sed 's/:/_/g')

        echo "Creating bam files for properly mapped reads in ${base} at ${line}..."
        echo samtools view -@ 16 -u -f 1 -F 12 ${bam} ${line} \
            > ${work}5kb_intermediatebams/${base}_${filename}_map_map.bam
        echo "Creating bam file for paired reads that are unmapped_mapped in ${base} at ${line}..."
        echo samtools view -@ 16 -u -f 4 -F 264 ${bam} ${line} \
            > ${work}5kb_intermediatebams/${base}_${filename}_unmap_map.bam
        echo "Creating bam file for paired reads that are mapped_unmapped in ${base} at ${line}..."
        echo samtools view -@ 16 -u -f 8 -F 260 ${bam} ${line} \
            > ${work}5kb_intermediatebams/${base}_${filename}_map_unmap.bam
            
        echo "Merging bam files for reads with unmapped pairs for ${base} in ${line}..."
        echo samtools merge -u ${work}5kb_intermediatebams/${base}_${filename}_unmap.bam \
            ${work}5kb_intermediatebams/${base}_${filename}_unmap_map.bam \
            ${work}5kb_intermediatebams/${base}_${filename}_map_unmap.bam
            
        echo samtools sort -@ 32 -n ${work}5kb_intermediatebams/${base}_${filename}_map_map.bam \
            > ${work}5kbbams/${base}_${filename}_mapped.bam
        echo samtools sort -@ 32 -n ${work}5kb_intermediatebams/${base}_${filename}_unmap.bam \
            > ${work}5kbbams/${base}_${filename}_unmapped.bam
            
        echo "Converting bam files to fastq files for ${base} in ${line}..."
        echo bamToFastq -i ${work}5kbbams/${base}_${filename}_mapped.bam \
            -fq ${work}/5kb_intermediatereads/${base}_${filename}_mapped.r1.fastq \
            -fq2 ${work}/5kb_intermediatereads/${base}_${filename}_mapped.r2.fastq
            
        echo bamToFastq -i ${work}5kbbams/${base}_${filename}_unmapped.sorted.bam \
            -fq ${work}/5kb_intermediatereads/${base}_${filename}_unmapped.r1.fastq \
            -fq2 ${work}/5kb_intermediatereads/${base}_${filename}_unmapped.r2.fastq
            
        echo cat ${work}/5kb_intermediatereads/${base}_${filename}_mapped.r1.fastq \
            ${work}/5kb_intermediatereads/${base}_${filename}_unmapped.r1.fastq \
            > ${work}/5kbreads/${base}_${filename}.r1.fastq
            
        echo cat ${work}/5kb_intermediatereads/${base}_${filename}_mapped.r2.fastq \
            ${work}/5kb_intermediatereads/${base}_${filename}_unmapped.r2.fastq \
            > ${work}/5kbreads/${base}_${filename}.r2.fastq
            
        echo "Now conducting local assembly for ${base} in region ${line}..."
        echo ${spades} \
             --careful \
             -o ${work}assemblies/5kb/${base}_${filename} \
             -1 ${work}/5kbreads/${base}_${filename}.r1.fastq \
             -2 ${work}/5kbreads/${base}_${filename}.r2.fastq \
             -m 64
    done < ${work}5kb_start_end_mantaSV.loopinput
done