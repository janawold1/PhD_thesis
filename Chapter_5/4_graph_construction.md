# Genome Graph Construction
## Identifying autosome 7 contigs
First did approximate alignment of contigs and found those that aligned to autosome 1.
```
data=/kakapo-data/ONT/assembly/
for male in Bill Blades Gulliver Rangi Sass
    do
    echo "Finding alignments to autosome 1 for ${male}..."
    wfmash --threads 24 -m ${mref} ${data}medaka/racon_polished/${male}/consensus.fasta > ${data}wfmash/${male}_raconMedaka_approx.paf
    wfmash --threads 24 -m ${mref} ${data}NextPolish/racon_NextPolish/${male}_polished/genome.nextpolish.fasta > ${data}wfmash/${male}_raconNextpolish_approx.paf
    grep NC_044283.2 ${data}wfmash/${male}_racon_approx.paf | awk '{print $1}' | sort | uniq > ${data}wfmash/${male}_raconMedaka_chr7.txt
    grep NC_044283.2 ${data}wfmash/${male}_raconNextpolish_approx.paf | awk '{print $1}' | sort | uniq > ${data}wfmash/${male}_raconNextpolish_chr7.txt
    ./faSomeRecords ${data}medaka/racon_polished/${male}/consensus.fasta ${data}wfmash/${male}_raconMedaka_chr7.txt ${data}chr7_scaffolds/${male}_raconMedaka_chr7.fa
    ./faSomeRecords ${data}NextPolish/racon_NextPolish/${male}_polished/genome.nextpolish.fasta ${data}wfmash/${male}_raconNextpolish_chr7.txt ${data}chr7_scaffolds/${male}_raconNextpolish_chr7.fa
done &
for female in C1 C2 Kuia Margaret-Maree Sue
do
    echo "Finding alignments to autosome 1 for ${female}..."
    wfmash --threads 24 -m ${fref} ${data}medaka/racon_polished/${female}/consensus.fasta > ${data}wfmash/${female}_raconMedaka_approx.paf
    wfmash --threads 24 -m ${fref} ${data}NextPolish/racon_NextPolish/${female}_polished/genome.nextpolish.fasta > ${data}wfmash/${female}_raconNextpolish_approx.paf
    grep NC_044283.2 ${data}wfmash/${female}_racon_approx.paf | awk '{print $1}' | sort | uniq > ${data}wfmash/${female}_raconMedaka_chr7.txt
    grep NC_044283.2 ${data}wfmash/${female}_raconNextpolish_approx.paf | awk '{print $1}' | sort | uniq > ${data}wfmash/${female}_raconNextpolish_chr7.txt
    ./faSomeRecords ${data}medaka/racon_polished/${female}/consensus.fasta ${data}wfmash/${female}_raconMedaka_chr7.txt ${data}chr7_scaffolds/${female}_raconMedaka_chr7.fa
    ./faSomeRecords ${data}NextPolish/racon_NextPolish/${female}_polished/genome.nextpolish.fasta ${data}wfmash/${female}_raconNextpolish_chr7.txt ${data}chr7_scaffolds/${female}_raconNextpolish_chr7.fa
done
```
## Estimating genome divergence
Distance of each query sequence to the reference was estimated using [mash](https://mash.readthedocs.io/en/latest/index.html) v2.3.
```
mash dist -p 24 ${data}assembly/chr1_racon_scaffolds/Jane_chr1.fa \
    ${data}assembly/chr7_racon_scaffolds/Bill_chr7.fa \
    ${data}assembly/chr7_racon_scaffolds/Blades_chr7.fa \
    ${data}assembly/chr7_racon_scaffolds/C1_chr7.fa \
    ${data}assembly/chr7_racon_scaffolds/C2_chr7.fa \
    ${data}assembly/chr7_racon_scaffolds/Gulliver_chr7.fa \
    ${data}assembly/chr7_racon_scaffolds/Margaret-Maree_racon_chr7.fa \
    ${data}assembly/chr7_racon_scaffolds/Sue_racon_chr7.fa 
```
The number of matching hashes for the racon & Medaka polished assemblies varied from 524 - 971 out of 1000, and the mash-distance varied from 0.0178312 - 0.000705841. Average genome divergence based off the MASH distances was 0.0052838601.

In contrast, the number of matching hashes for the racon & Medaka polished assemblies varied from 524 - 971 out of 1000, and the mash-distance varied from 0.0178312 - 0.000705841. Average genome divergence based off the MASH distances was 0.0052838601.

## Running PGGB
For pangenome graph construction assembled contigs for all samples were named in accordance with the naming convention:
<sample_name>#<contig_ID>

Files were then concatenated into a single fasta and a graph constructed with [pggb](https://github.com/pangenome/pggb) v0.2.0.
```
for indiv in Bill Blades C1 C2 Gulliver Jane Margaret-Maree Sue
    do
    echo ${indiv}
    cat chr_7/${indiv}_racon_chr7.fa >> chr_7/chr7_filtered.fa
done

for fa in chr7_filtered chr7_asJane
    do
    if (( "$fa" == chr7_asJane ))
    then
        for segment in 3000 10000 50000 100000
        do
            pggb -i ${fa}.fa -o ${fa}_defaultk_${segment}bpSegment -t 24 -p 98 -s $segment -n 8 -T 24 -U -v -L -V 'Jane:#' -m
            pggb -i ${fa}.fa -o ${fa}_k100_${segment}bpSegment -t 24 -p 98 -s $segment -k 100 -n 8 -T 24 -U -v -L -V 'Jane:#' -m
            pggb -i ${fa}.fa -o ${fa}_k311_${segment}bpSegment -t 24 -p 98 -s $segment -k 311 -n 8 -T 24 -U -v -L -V 'Jane:#' -m
        done
    else
        for seg in 10000 50000 100000
        do
            pggb -i ${fa}.fa -o ${fa}_defaultk_${segment}bpSegment -t 24 -p 98 -s $segment -n 8 -T 24 -U -v -L -V 'Jane:#' -m
            pggb -i ${fa}.fa -o ${fa}_k100_${segment}bpSegment -t 24 -p 98 -s $segment -k 100 -n 8 -T 24 -U -v -L -V 'Jane:#' -m
            pggb -i ${fa}.fa -o ${fa}_k311_${segment}bpSegment -t 24 -p 98 -s $segment -k 311 -n 8 -T 24 -U -v -L -V 'Jane:#' -m
        done
    fi
done
```
Graph visualisation and manipulation were conducted using [odgi](https://github.com/pangenome/odgi) v0.6.3. 

After examining overall topology, number of edges, number of nodes, the 

Finding regions of high complexity by estimating the depth over the pangenome.
```
graph_dir=/kakapo-data/ONT/assembly/pggb_raconNextpolish/chr7_asJane_k311_100000bpSegment/
graph_og=/kakapo-data/ONT/assembly/pggb_raconNextpolish/chr7_asJane_k311_100000bpSegment/chr7.fa.9d4992c.4030258.a7f493c.smooth.og

odgi depth -i ${graph_og} -r Jane#chromosome_7 | bedtools makewindows -b /dev/stdin -w 5000 > ${graph_dir}chr7.w5kb.bed
odgi depth -i ${graph_og} -b ${graph_dir}chr7.w5kb.bed --threads 2 | bedtools sort > ${graph_dir}chr7.depth.w5kb.bed
```
## Examining TLRs on autosome 7
Variability in toll-like receptors 1, 2, 2 type 2, 3, and 6 were then visualised in 1D graphs.

```
for dir in chr7_*/
    do
    echo "Running ODGI for $dir..."
    odgi viz -i ${dir}chr7*.smooth.og -o TLR1.png -x 500 -bz -r Jane#chromosome_7:24009134-25000000 &
    odgi viz -i ${dir}chr7*.smooth.og -o TLR2.png -x 500 -bz -r Jane#chromosome_7:44500000-46500000 &
    odgi viz -i ${dir}chr7*.smooth.og -o TLR3.png -x 500 -bz -r Jane#chromosome_7:32400000-34400000 &
    odgi viz -i ${dir}chr7*.smooth.og -o TLR6.png -x 500 -bz -r Jane#chromosome_7:22900000-24009134
    wait
done
```
Extracting sub-graphs for each TLR.
```
odgi extract --threads 24 --progress -i ${graph_og} -o TLR1_subgraph.og -r Jane#chromosome_7:24009134-25000000 -E
odgi extract --threads 24 --progress -i ${graph_og} -o TLR2_subgraph.og -r Jane#chromosome_7:44500000-46500000 -E
odgi extract --threads 24 --progress -i ${graph_og} -o TLR3_subgraph.og -r Jane#chromosome_7:32400000-34400000 -E
odgi extract --threads 24 --progress -i ${graph_og} -o TLR6_subgraph.og -r Jane#chromosome_7:22900000-24009134 -E

odgi stats -i TLR1_subgraph.og -S
odgi stats -i TLR2_subgraph.og -S
odgi stats -i TLR3_subgraph.og -S
odgi stats -i TLR6_subgraph.og -S

odgi sort -i TLR1_subgraph.og -o - -O | odgi viz -i - -o TLR1_subgraph.png -s '#' -P
odgi sort -i TLR2_subgraph.og -o - -O | odgi viz -i - -o TLR2_subgraph.png -s '#' -P
odgi sort -i TLR3_subgraph.og -o - -O | odgi viz -i - -o TLR3_subgraph.png -s '#' -P
odgi sort -i TLR6_subgraph.og -o - -O | odgi viz -i - -o TLR6_subgraph.png -s '#' -P

odgi view -i TLR1_subgraph.og -g > TLR1_subgraph.gfa
odgi view -i TLR2_subgraph.og -g > TLR2_subgraph.gfa
odgi view -i TLR3_subgraph.og -g > TLR3_subgraph.gfa
odgi view -i TLR6_subgraph.og -g > TLR6_subgraph.gfa
```
Estimating depth along each of the TLRs
```
odgi depth -i TLR1_subgraph.og -r Jane#chromosome_7:24008847-25000092 | bedtools makewindows -b /dev/stdin -w 1000 > TLR1.w1kb.bed
odgi depth -i TLR2_subgraph.og -r Jane#chromosome_7:44499742-46500132 | bedtools makewindows -b /dev/stdin -w 1000 > TLR2.w1kb.bed
odgi depth -i TLR3_subgraph.og -r Jane#chromosome_7:32399822-34400344 | bedtools makewindows -b /dev/stdin -w 1000 > TLR3.w1kb.bed
odgi depth -i TLR6_subgraph.og -r Jane#chromosome_7:22899873-24009152 | bedtools makewindows -b /dev/stdin -w 1000 > TLR6.w1kb.bed

odgi depth -i TLR1_subgraph.og -b ${graph_dir}TLR1.w1kb.bed --threads 2 | bedtools sort > ${graph_dir}TLR1.depth.w1kb.bed
odgi depth -i TLR2_subgraph.og -b ${graph_dir}TLR2.w1kb.bed --threads 2 | bedtools sort > ${graph_dir}TLR2.depth.w1kb.bed
odgi depth -i TLR3_subgraph.og -b ${graph_dir}TLR3.w1kb.bed --threads 2 | bedtools sort > ${graph_dir}TLR3.depth.w1kb.bed
odgi depth -i TLR6_subgraph.og -b ${graph_dir}TLR6.w1kb.bed --threads 2 | bedtools sort > ${graph_dir}TLR6.depth.w1kb.bed
```