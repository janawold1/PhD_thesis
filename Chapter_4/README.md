# Chapter 4: Structural variant discovery in the kākāpō125+ data set

Here are the steps I took to trial different strategies for SV discovery using illumina sequence data. The raw fasta files were preprocessed, but averaged around 25x sequence coverage. 

## 1 Read alignment with BWA
Reads were aligned to the kākāpō reference genome (VGP) using bwa. The W sex chromosome (present in females) was removed prior to alignments of males as the W is highly repetitive, and may contain regions that would result in misalignment of reads. 

## 2 Estimates of mapping statistics
Three pipelines were used for SV discovery, Delly, Manta and SMOOVE. SV discovery was conducted for each replicate of down sampled alignments and the total alignment.

    Programs used BWA, SAMtools (link to their githubs)

## 3 Structural Variant discovery using Delly
Called and genotyped SVs with Delly.

    Programs used Delly, BCFtoools

## 4 Structural Variant discovery using Smoove
Called and genotyped SVs with the Smoove pipeline.

    Programs used Smoove, BCFtools.

## 5 Structural Variant discovery using Manta
Called SVs using two strategies 1) Joint and 2) Batched.

    Programs used Manta, BCFtools

## 6 Genotyping Manta outputs with BayesTyper
Genotyped SVs using BayesTyper for both Joint and Batched Manta outputs

    Programs used BayesTyper, KMC, BCFtools

## 7 Identified consensus SV calls with Survivor
Identified overlapping SVs with Survivor

    Progams used Survivor, BCFtools