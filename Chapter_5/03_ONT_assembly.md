# Brief trial script for running flye.


```
data=/kakapo-data/ONT/
aves=aves_odb10
source ~/anaconda3/etc/profile.d/conda.sh
```

```
conda activate flye
for i in {10}
    do
    for fasta in ${data}trimmed/*_q10_10kbtrim.fastq.gz
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

```
conda activate busco
echo "Running first round of BUSCO checks..."
cd ${data}busco
for assembly in ${data}fly_out/Blades*
       do
       base=$(basename $assembly)
       quast -t 24 -o ${data}quast/${base} -r /kakapo-data/References/kakapo_no_Wchromosome.fa ${assembly}/assembly.fasta
       busco --in ${assembly}/assembly.fasta --out ${base} --mode genome --lineage_dataset $aves
done
```