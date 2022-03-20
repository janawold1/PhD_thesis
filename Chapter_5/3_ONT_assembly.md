# Genome assembly with MinION long-reads
Here are the initial stages of curating reference assemblies for 10 kākāpō. These individuals are inclusive of 5 males and 5 females. All males sequenced here and 3 of the females are founders. The remaining 2 females are relatively young birds hatched....

Global variables were defined as:
```
data=/kakapo-data/ONT/
aves=aves_odb10
source ~/anaconda3/etc/profile.d/conda.sh
```
As per work done by Annabel Whibley at Auckland University, we trialled different minimum read lengths such as 3, 4, 5, and 10 kb in length.
```
conda activate flye
for i in {3,4,5,10}
    do
    for fasta in ${data}trimmed/*_q10_${i}kbtrim.fastq.gz
        do
        base=$(basename ${fasta} .fastq.gz)
        printf "\nRunning flye for $base...\n"
        flye --nano-raw ${fasta} --out-dir ${data}flye_out${base} --genome-size 1.2g \
            --threads 24 --debug
    done &
done
wait
conda deactivate flye
```
[BUSCO](https://busco.ezlab.org/)v5.3.0 scores and [quast](https://github.com/ablab/quast) were run for each assembly were used to identify which assemblies returned the highest scores. Consistent with Dr. Whibley's results, I found that a minimum length of 5kb was a great cutoff across the board.
```
conda activate busco
mref=/kakapo-data/References/kakapo_no_Wchromosome.fa
fref=/kakapo-data/References/kakapo_full_ref.fa
echo "Running first round of BUSCO and quast checks..."
cd ${data}assembly/busco
for male in Bill Blades Gulliver Rangi Sass Smoko
       do
       base=$(basename $assembly)
       quast -t 24 -o ${data}assembly/quast/initial_assemblies/${male} -r $mref ${data}assembly/initial_assemblies/${male}.fasta
       busco --in ${data}assembly/initial_assemblies/${male}.fasta --out ${male}_initial --mode genome --lineage_dataset $aves
done
for female in Bella C1 C2 Kuia Margaret-Maree Sue
    echo "Running BUSCO and quast for ${female}..."
    quast -t 24 -o ${data}assembly/quast/initial_assemblies/${female} -r $fref ${data}assembly/initial_assemblies/${female}.fasta
    busco --in ${data}assembly/initial_assemblies/${female}.fasta --out ${female}_initial --mode genome --lineage_dataset $aves
done
```
## Initial polishing
Polishing for each assembly was conducted using [NextPolish](https://github.com/Nextomics/NextPolish) v1.4.0 and [racon]() v1.4.20. 
Both short-read and long read data were used in polishing with NextPolish for all samples excluding C1 and C2. Only long-read data was used with racon.

Config files for each individual were generated in accordance with these instructions. Upon completion of the initial round of polishing, assembly quality was assessed once more with BUSCO and quast.


Polishing was also conducted with [NextPolish](https://github.com/Nextomics/NextPolish) Config files included short- and long-reads. 
```
mkdir -p ${data}assembly/NextPolish/cfgs

for indiv in Bill Blades Gulliver Kuia Margaret-Maree Rangi Sass Sue
    do
    ls ${data}trimmed/${indiv}_q10_5kbtrim.fastq.gz > ${data}assembly/NextPolish/cfg/${indiv}_lgs.fofn
    ls ${data}assembly/short-reads/${indiv}* > ${data}assembly/NextPolish/cfg/${indiv}_sgs.fofn
done
```

Finally, polishing with Racon was run as per:
```
mkdir -p ${data}assembly/racon{1,2}/{alignments,polishing}

for samp in Bill Blades C1 C2 Gulliver Kuia Margaret-Maree Rangi Sass Sue
    do
    echo "Aligning reads with minimap for ${samp}..."
    minimap2 -ax map-ont ${data}assembly/initial_assemblies/${samp}.fasta ${data}trimmed/${samp}_q10_5kbtrim.fastq.gz > ${data}assembly/racon/racon1/alignments/${samp}.sam
    echo "Running RACON for $samp..."
    racon -t 48 ${data}trimmed/${samp}_q10_5kbtrim.fastq.gz ${data}assembly/racon/racon1/alignments/${samp}.sam ${data}assembly/initial_assemblies/${samp}.fasta > ${data}assembly/racon/racon1/polished_assemblies/${samp}_racon1.fasta
    echo "Running second round of polishing for ${samp}..."
    minimap2 -ax map-ont ${data}aassembly/racon/racon1/polished_assemblies/${samp}_racon1.fasta ${data}trimmed/${samp}_q10_5kbtrim.fastq.gz > ${data}assembly/racon/racon2/alignments/${samp}.sam
    racon -t 48 ${data}trimmed/${samp}_q10_5kbtrim.fastq.gz ${data}assembly/racon/racon2/alignments/${samp}.sam ${data}assembly/racon/racon1/polished_assemblies/${samp}_racon1.fasta > ${data}assembly/racon/racon2/polished_assemblies/${samp}_racon2.fasta
done
```
Quality of these polishing steps were assessed with:
```
cd ${data}assembly/
for indiv in Bill Blades C1 C2 Gulliver Kuia Margaret-Maree Rangi Sass Sue
    do
    echo "Running BUSCO v5.3.0 for $indiv..."
    busco --in racon/racon1/polished_assemblies/${indiv}_racon1.fasta --out busco/busco_racon1/${indiv} --mode genome --lineage_dataset $aves
    busco --in racon/racon2/polished_assemblies/${indiv}_racon2.fasta --out busco/busco_racon2/${indiv} --mode genome --lineage_dataset $aves
    busco --in NextPolish/${indiv}_polished/genome.nextpolish.fasta --out busco/busco_NextPolish1/${indiv} --mode genome --lineage_dataset $aves
done &
for female in C1 C2 Kuia Margaret-Maree Sue
    do
    echo "Running QUAST for ${female}..."
    quast -t 24 -o quast/racon1/${female} -r $fref racon/racon1/${female}.fasta
    quast -t 24 -o quast/racon2/${female} -r $fref racon/racon2/${female}.fasta
    quast -t 24 -o quast/NextPolish1/${female} -r $fref NextPolish/${female}_polished/genome.nextpolish.fasta
done &
for male in Bill Blades Gulliver Rangi Sass
    do
    quast -t 24 -o quast/2_racon1/${male} -r $mref racon/racon1/${male}.fasta 
    quast -t 24 -o quast/3_racon2/${male} -r $mref racon/racon2/${male}.fasta
    quast -t 24 -o quast/2_NextPolish1/${male} -r $mref NextPolish/${male}_polished/genome.nextpolish.fasta
done
```
## Final polishing
A final round of polishing was conducted with [medaka](https://github.com/nanoporetech/medaka) v1.5.0.
```
for indiv in Bill Blades C1 C2 Gulliver Kuia Margaret-Maree Rangi Sass Sue
    do
    echo "Running MEDAKA for ${indiv}...."
    medaka_consensus -i ../trimmed/${indiv}_q10_5kbtrim.fastq.gz -d racon/racon2/polished_assemblies/${indiv}_racon2.fasta -o medaka/racon_polished/${indiv} -t 48 -m r941_min_sup_g507
    medaka_consensus -i ../trimmed/T{indiv}_q10_5kbtim.fastq.gz -d NextPolish/${indiv}_polished/genome.nextpolish.fasta -o medaka/NextPolish_polished/${indiv}
done

for male in Bill Blades Gulliver Rangi Sass
    do
    echo "Running QUAST for ${indiv}...."
    quast -t 24 -o quast/4_medaka_racon/${male} -r ${mref} medaka/racon_polished/${male}/consensus.fasta &
    quast -t 24 -o quast/4_medaka_NextPolish/${male} -r ${mref} medaka/NextPolish_polished/${male}/consensus.fasta
    wait
done
for female in C1 C2 Kuia Margaret-Maree Sue
    do
    quast -t 24 -o quast/4_medaka_racon/${female} -r ${fref} medaka/racon_polished/${female}/consensus.fasta &
    quast -t 24 -o quast/4_medaka_NextPolish/${female} -r ${fref} medaka/NextPolish_polished/${female}/consensus.fasta
    wait
done
quast -t 24 -o quast/Jane ${fref}
wait
for indiv in Bill Blades C1 C2 Gulliver Kuia Margaret-Maree Rangi Sass Sue
    do
    busco --in medaka/racon_polished/${indiv}/consensus.fasta --out busco/busco_racon/${indiv} --cpu 24 --mode genome --lineage_dataset ${aves} &
    busco --in medaka/NextPolish_polished/${indiv}/consensus.fasta --out busco/busco_NextPolish/${indiv} --cpu 24 --mode genome --lineage_dataset ${aves}
    wait
done
```
