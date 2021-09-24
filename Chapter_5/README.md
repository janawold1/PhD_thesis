# Chapter 5: Kākāpō long-read sequencing, genome assembly and SV discovery
Here are the steps I took to trial different strategies for SV discovery using illumina sequence data. The raw fasta files were preprocessed, but averaged around 25x sequence coverage. These reads were then down-sampled to approximately 10x sequence coverage for a second trial set (bbmap). Down-sampling was conducted in triplicate to assess variability in reads retained for alignment. 

## Alignment
Reads were aligned to the kākāpō reference genome (VGP) using bwa. The W sex chromosome (present in females) was removed prior to alignments of males as the W is highly repetitive, and may contain regions that would result in misalignment of reads. 

## SV discovery
Three pipelines were used for SV discovery, Delly, Manta and SMOOVE. SV discovery was conducted for each replicate of down sampled alignments and the total alignment.

## SV filtering
