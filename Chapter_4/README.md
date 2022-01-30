# Chapter 4: Structural variant discovery in the kākāpō125+ data set

Here are the steps I took to trial different strategies for SV discovery using illumina sequence data. The raw fasta files were preprocessed, but averaged around 25x sequence coverage. Once mapped, this coverage ranged from 9x - 38x with a mean of ~18x. 

Three pipelines were used for SV discovery, Delly, Manta and SMOOVE. SV discovery was conducted for each replicate of down sampled alignments and the total alignment.

## 1 Read alignment with BWA
Reads were aligned to the kākāpō reference genome (VGP) using [BWA](http://bio-bwa.sourceforge.net/). The W sex chromosome (present in females) was removed prior to alignments of males as the W is highly repetitive, and may contain regions that would result in misalignment of reads. Other programs used in this step include: [SAMtools](https://github.com/samtools/samtools)

## 2 Estimates of mapping statistics
Programs used [Mosdepth](https://github.com/brentp/mosdepth), [qualimap](http://qualimap.conesalab.org/) and [SAMtools](https://github.com/samtools/samtools)

## 3 Structural Variant discovery using [Delly](https://github.com/dellytools/delly)
Called and genotyped SVs with Delly. Other programs used in this step include [BCFtoools](http://samtools.github.io/bcftools/)

## 4 Structural Variant discovery using [Smoove](https://github.com/brentp/smoove)
Called and genotyped SVs with the Smoove pipeline. Other programs used in this step include [BCFtoools](http://samtools.github.io/bcftools/).

## 5 Structural Variant discovery using [Manta](https://github.com/Illumina/manta)
Called SVs using two strategies 1) Joint and 2) Batched. Other programs used include BCFtools

## 6 Genotyping Manta outputs with [BayesTyper](https://github.com/bioinformatics-centre/BayesTyper)
Genotyped SVs using BayesTyper for both Joint and Batched Manta outputs. Programs used include [BCFtoools](http://samtools.github.io/bcftools/) and [KMC](https://github.com/refresh-bio/KMC)

## 7 Identified consensus SV calls with [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
Identified overlapping SVs with SURVIVOR. Other progams used:  [BCFtoools](http://samtools.github.io/bcftools/)