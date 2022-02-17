# *De novo* assembly of kākāpō founders

```
data=/kakapo-data/ONT/fastq/
aves=aves_odb10
source ~/anaconda3/etc/profile.d/conda.sh
```
## Initial Flye assemblies
```
conda activate flye
for fasta in ${data}*.fastq.gz
do
base=$(basename ${fasta} .fastq.gz)
#mkdir ${out}${base}
printf "\nRunning flye for $base...\n"
flye --nano-raw ${fasta} --out-dir ${out}${base} --genome-size 1.2g --threads 16 --debug
done
conda deactivate flye
```
## BUSCO scores

```
conda activate busco
cd ${data}busco

for asm in ${data}flye_out/*_q10_5000trim/
    do
    base=$(basename $asm)
    echo "Starting BUSCO for ${base}..."
    busco --in ${asm}assembly.fasta --out ${base} --mode genome --lineage_dataset ${aves} &
done
wait
```