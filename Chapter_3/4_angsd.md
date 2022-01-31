# Step 4: Running ANGSD for tara iti and Australian fairy tern
To start, I defined some global variables and activated the Conda environment for ANGSD.
```
data=/data/common_tern/angsd/
ref=/data/reference/bSteHir1.pri.cur.20190820.fasta
ngstools=~/Programs/ngsTools/
source ~/anaconda3/etc/profile.d/conda.sh
conda activate angsd
```
Then I began a QC check for each the tara iti, Australian and global sample sets. The 'global' data includes both the tara iti and Australian samples. For each group, only paired-reads that uniquely mapped to the reference with a min MapQ score of 20 and a max read depth of 150x were included. These settings will be used for all analyses going forward. 
```
mkdir ${data}qc
echo "Estimating group QC metrics..."
for list in ${data}*.bamlist
        do
        group=$(basename ${list} .bamlist)
        angsd -b ${list} -out ${data}qc/${group}.qc -ref ${ref} \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doQsDist 1 \
                -doDepth 1 -doCounts 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -maxDepth 150 &
done
wait
```
Once the QC run was completed, individual plots were generated using the NGStools Rscript. Low quality samples were identified and excluded from further analyses.
```
for qc in ${data}qc/*.qc
        do
        base=$(basename ${qc} .qc
        echo "Plotting QC for ${base}..."
        Rscript ${ngsTools}Scripts/plotQC.R ${qc}
done
```