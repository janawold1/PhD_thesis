# Chapter 4: Incorporating structural variants into the conservation management of Aotearoa New Zealand's most threatened breeding bird, tara iti (*Sternula nereis davisae*)

Here are the methods I used for SNP discovery and exploring relative levels of inbreeding and differentiation between Australian fairy tern (*Sternula nereis nereis*) and tara iti (*S. n. davisae*).

## 1 Read trimming and read quality statistics
After assessing the quality of reads with [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.9, the [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) v0.6.7 package was used to trim raw Illumina data prior to alignment to the reference genome. Trimming was conducted under the two-colour option for Nextera adapters. 

## 2 Read alignment with BWA
Following read trimming, reads were aligned to the reference genome using [BWA](https://github.com/lh3/bwa) v0.7.17 and [SAMtools](https://github.com/samtools/samtools) v1.12 was used to sort and process alignments prior to calculating mapping statistics.

## 3 Mapping statistics
The quality of read mapping was assessed using [SAMtools](https://github.com/samtools/samtools) v1.12, [mosdepth](https://github.com/brentp/mosdepth) v0.3.3 and [qualimap](http://qualimap.conesalab.org/) v2.2.2.

## 4 SNP calling pipeline
A pipeline provided by Roger Moraga of Tea Break Bioinformatics was used for discovery of single nucleotide polymorphisms (SNPs). This pipeline speeds up ```bcftools mpileup``` by segregating regions of the genome into batches that may be run in parallel. [BCFtools](https://github.com/samtools/bcftools) v1.12 was used for SNP discovery.

## 5 SNP filtering parameters
[VCFtools](http://vcftools.sourceforge.net/man_latest.html) v0.1.16 was used to explore different filtering thresholds for SNPs. Outputs were graphed in R using the ggplot2 package.

## 6 SNP analyses
Exploratory analyses of SNPs were conducted with [pixy](https://github.com/ksamuk/pixy) v1.2.6, [VCFtools](http://vcftools.sourceforge.net/man_latest.html), and [ADMIXTURE](https://bioconda.github.io/recipes/admixture/README.html?highlight=admixture#package-package%20&#x27;admixture&#x27;) v1.3.0.

## 7 SV discovery and filtering with Delly
Structural variant (SV) discovery was conducted with [Delly](https://github.com/dellytools/delly) v0.8.7 as per [Chapter 3](https://github.com/janawold1/PhD_thesis/tree/main/Chapter_3)

## 8 SV discovery and filtering with SMOOVE
SV discovery was conducted in [SMOOVE](https://github.com/brentp/smoove) v0.2.8 as per [Chapter 3](https://github.com/janawold1/PhD_thesis/tree/main/Chapter_3). SV filtering was modified as an annotation for the reference genome used here is not yet available. 

## 9 Exploring SNP and SV data in R
Outputs from [pixy](https://github.com/ksamuk/pixy) and [VCFtools](http://vcftools.sourceforge.net/man_latest.html) were graphed. Population structure and differentiation were inferred from pairwise Fst and PCAs in the hierfstat and adegenet R packages. 