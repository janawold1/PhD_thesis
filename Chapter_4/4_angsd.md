# Step 4: Running ANGSD for tara iti and Australian fairy tern
The steps outlined here for the most part follow a useful [Tutorial](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md) outlining best practice of using [ANGSD](https://github.com/ANGSD/angsd) and [ngsTools](https://github.com/mfumagalli/ngsTools). ANGSD can be a very challenging program to run due to unclear documentation and dependency clashes (even for the Conda installation). The Conda installation worked well for most ANGSD steps, however the realSFS progam would not work. Other users reported similar challenges, with the only [solution](https://github.com/ANGSD/angsd/issues/396#issuecomment-1004385423) being to use a [singularity container](https://cloud.sylabs.io/library/james-s-santangelo/angsd/angsd) provided by James Santangelo for ANGSD version 0.9 was used.

To start with the Conda installed ANGSD on the University of Canterbury Research Compute Cluster (RCC), I defined some global variables and activated the Conda environment for ANGSD.
```
data=/data/common_tern/angsd/
ref=/data/reference/bSteHir1.pri.cur.20190820.fasta
ngstools=~/Programs/ngsTools/
source ~/anaconda3/etc/profile.d/conda.sh
conda activate angsd
```
# QC of aligned reads
Then I began a QC check for each the tara iti (TI), Australian fairy tern (AU) and global sample sets. The 'global' data includes both the tara iti and Australian samples. For each group, only paired-reads that uniquely mapped to the reference with a max read depth of 150x were included. The ```-C 50``` setting corrects for the number of read mismatches. 
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
Once the QC run was completed, individual plots were generated using the NGStools Rscript. 
```
for qc in ${data}qc/*.qc
        do
        base=$(basename ${qc} .qc
        echo "Plotting QC for ${base}..."
        Rscript ${ngsTools}Scripts/plotQC.R ${qc} &
done
wait
```
For inclusion in further analyses, individual read depth distributions had to peak at a coverage >=5x. This left 11 TI and 15 AU samples. Also found that a minimum mapping quality and minimum base quality score of 20 works well for all 3 data sets. 

# Inbreeding Estimates
Given the length and severity of population bottleneck experienced by the tara iti population, the assumption of HWE is likely not valid. In instances like these, it is recommended that the priors for genotyping liklihoods are adjusted with estimates of inbreeding coefficients to ensure inferences about genetic distance and population structure take this into account. 

To start, genotype liklihoods were estimated as per:

```
for POP in AU TI global
        do
        if [[ "$POP" == AU ]]
                then
                echo "Estimating inbreeding for $POP...."
                angsd -P 24 -b ${POP}.bamlist -ref ${ref} -out ${data}inbreeding/${POP} \
                	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                	-minMapQ 20 -minQ 20 -minInd 14 -setMinDepth 5 -setMaxDepth 150 -doCounts 1 \
        	        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        	        -doGlf 3 -SNP_pval 1e-3
        elif [[ $POP == TI ]]
                then
                angsd -P 24 -b ${POP}.bamlist -ref ${ref} -out ${data}inbreeding/${POP} \
                        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                        -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 5 -setMaxDepth 150 -doCounts 1 \
                        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
                        -doGlf 3 -SNP_pval 1e-3
        else
                angsd -P 24 -b ${POP}.bamlist -ref ${ref} -out ${data}inbreeding/${POP} \
        	        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        	        -minMapQ 20 -minQ 20 -minInd 23 -setMinDepth 5 -setMaxDepth 150 -doCounts 1 \
                        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
                        -doGlf 3 -SNP_pval 1e-3
        fi &
done
wait
```
These genotype liklihoods were then used as input for [ngsF](https://github.com/fgvieira/ngsF) from the [ngsTools](https://github.com/mfumagalli/ngsTools) suite.

Although genotype likelihoods were estimated for the 'global' data set, inbreeding coefficients were estimated for the AU and TI populations independently.

```
for POP in AU TI
        do
        nsamp=$(wc -l ${data}${POP}.bamlist | awk '{print $1}')
        nsites=$(zcat ${data}inbreeding/${POP.mafs.gz | tail -n+2 | wc -l)
        echo "Preparing GLF for $POP for $nsamp samples and $nsites sites..."
        zcat ${data}inbreeding/${POP}.glf.gz > ${data}inbreeding/${POP}.glf
        ${ngstools}ngsF/ngsF.sh --n_ind $nsamp --n_sites $nsites --glf ${data}inbreeding/${POP}.glf --out ${data}inbreeding/${POP}.indF &
done
wait
```
These initial inbreeding coefficients were then used to inform the posterior probabilities for genotype and allele frequency as below. Here, [ANGSD]() was used to estimate the site frequency spectrum (SFS), which is the proportion of site allele frequencies. Since reads were aligned to the [Common tern](https://www.ncbi.nlm.nih.gov/genome/?term=common+tern) reference genome, the ancestral state could not be used to estimate an unfolded SFS, limiting inferences about population demography or selective events. 
```
for POP in AU TI global
        do
        if [[ "$POP" == AU ]]
                then
                echo "Estimating inbreeding for $POP...."
                angsd -P 24 -b ${POP}.bamlist -ref ${ref} -anc ${ref} -out ${data}inbreeding/${POP}.inbred \
                	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                	-minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 5 -setMaxDepth 150 -doCounts 1 \
        	        -GL 1 -doMajorMinor 1 -doMaf -1 -skipTriallelic 1 \
        	        -SNP_pval 1e-3 -doGeno 20 -doPost 1 -doSaf 2 -indF ${data}inbreeding/${POP}.indF
        else
                angsd -P 24 -b ${POP}.bamlist -ref ${ref} -anc ${ref} -out ${data}inbreeding/${POP}.inbred \
                        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                        -minMapQ 20 -minQ 20 -minInd 11 -setMinDepth 5 -setMaxDepth 150 -doCounts 1 \
                        -GL 1 -doMajorMinor 1 -doMaf -1 -skipTriallelic 1 \
                        -SNP_pval 1e-3 -doGeno 20 -doPost 1 -doSaf 2 -indF ${data}inbreeding/${POP}.indF
        fi &
done
wait
```

Once the posterior probabilities (informed by inbreeding coefficients) of sample allele frequencies (SAF) were computed as above, the SFS is computed with ```realSFS```.

Initially examined the output with:
```
realSFS print ${data}inbreeding/AU.saf.idx | less -S
realSFS print ${data}inbreeding/TI.saf.idx | less -S
```
And the overall SFS was estimated and plotted as per:
```
for POP in AU TI
        do
        echo "Estimating SFS for $POP..."
        realSFS ${data}inbreeding/${POP}.saf.idx > ${data}SFS/${POP}.sfs 
done
wait

Rscript ${ngstools}Scripts/plotSFS.R ${data}SFS/AU.sfs-${data}SFS/TI.sfs AU-TI 1 ${data}SFS/ALL.sfs.pdf
```
# Summary Statistics 

# Genetic Distance

# Population Structure
To visualise global population strucgture (i.e., Structure between AU and TI populations) and sub-population structure within AU and TI, I performed a Principal Component Analysis (PCA) and Mulidimensional Scaling (MDS) analysis using genetic distance. Genotype probabilities for each site were estimated using per-site allele frequency as a prior with ```ANGSD```, ```ngsF``` and ```realSFS```, as above. 
## Principal Component Analysis (PCA)
To estimate a PCA, the ngsCovar program was used. 
```
gunzip ${data}SFS/*.geno.gz

for POP in AU TI
        do
        nsites=$(zcat ${data}inbreeding/${POP}.mafs.gz | tail -n+2 | wc -l)
        nind=$(cat ${data}${POP}.bamlist | wc -l)
        ${ngstools}ngsPopGen/ngsCovar -probfile ${data}SFS/$POP.geno -outfile ${data}distance/${POP}.covar -nind
done

Rscript -e 'write.table(cbind(seq(1,30),rep(1,30),c(rep("LWK",10),rep("TSI",10),rep("PEL",10))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="Results/ALL.clst", quote=F)'
Rscript $NGSTOOLS/Scripts/plotPCA.R -i Results/ALL.covar -c 1-2 -a Results/ALL.clst -o Results/ALL.pca.pdf

```
## Multi Dimensional Scaling (MDS) plot
```

```
## Population Admixture 
```

```
# Relatedness Estimates for population comparison
```

```