# Challenges and Opportunities for Integrating Structural Variants into Conservation Genomics
Here are the scripts used for analyses in each Chapter of my thesis. Within each directory is a brief over view of chapters, and all scripts used in formal analysis in Chapters 3, 4, and 5. 

## Chapter 1: General Introduction
In Chapter 1 I provide a brief overview of conservation genomics and highlight the opportunities to enhance conservation outcomes by integrating structural variants (SVs) into conservation genomics studies. 

## Chapter 2: Expanding the conervation genomics toolbox: Incorporating structural variants to enhance genomic studies for species of conservation concern
Published Perspective in the *Molecular Ecology* Special Issue: [Whole-genome sequencing in molecular ecology](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16141).

## Chapter 3: The promise and challenges of structural variant discovery in the critically endangered kākāpō (*Strigops habroptilus*)
A nearly complete [short-read resequence data set](https://www.doc.govt.nz/our-work/kakapo-recovery/what-we-do/research-for-the-future/kakapo125-gene-sequencing/) and [high quality reference genome](https://vgp.github.io/genomeark/Strigops_habroptilus/) were used to investigate structual variantion (SVs) in the critically endangered kākāpō (*Strigops habroptilus*). Structural variants were called using the [Delly](https://github.com/dellytools/delly), [Manta](https://github.com/Illumina/manta) and [SMOOVE](https://github.com/brentp/smoove) pipelines.

## Chapter 4: Incorporating structural variants into the conservation management of Aotearoa New Zealand's most threatened breeding bird, tara iti (*Sternula nereis davisae*)
A population genomics study of the threatened Australian fairy tern (*Sternula nereis nereis*) and the Nationally Critical tara iti (*S. n. davisae*). Here, we aligned whole-genome short-read resequence data to the chromosomally assembled reference genome for the [common tern](https://vgp.github.io/genomeark/Sterna_hirundo/) (*Sterna hirundo*) to explore relative levels of population diversity and differentiation for these threatened taxa.

## Chapter 5: The promise of genome graphs as a conservation genomics resource for the critically endangered kākāpō (*Strigops habroptilus*)
In this chapter I generate a genome graph for kākāpō chromosome 7, which has four toll-like receptor genes using a *de novo* assembly approach to graph construction. Long-read sequence data was generated using an Oxford Nanopore MinION for 12 individuals. Genomes were assembled using [FLYE](https://github.com/fenderglass/Flye) v2.8.3 and multiple strategies to assembly polishing trialled. Graphs were constructed using [pggb](https://github.com/pangenome/pggb) v0.2.0, graph statistics estimated with [odgi](https://github.com/pangenome/odgi) v0.6.3 and visualised with [MultiQC](https://github.com/ewels/MultiQC) v1.11. Graphs were also processed with odgi and the final outputs examined with [Bandage](https://rrwick.github.io/Bandage/) v0.9.0.

## Chapter 6: General Discussion
Finally, I close with a summary of results and identify future directions for research.