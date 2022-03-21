# Challenges and Opportunities for Integrating Structural Variants into Conservation Genomics
Presented here are the scripts used for data analyses in each data chapter of my thesis. Within each directory is a brief overview of chapter contents and details about each of the tools used. 

For a brief outline of chapter contents:

## Chapter 3: The promise and challenges of structural variant discovery in the critically endangered kākāpō (*Strigops habroptilus*)
A nearly complete [short-read resequence data set](https://www.doc.govt.nz/our-work/kakapo-recovery/what-we-do/research-for-the-future/kakapo125-gene-sequencing/) and [high quality reference genome](https://vgp.github.io/genomeark/Strigops_habroptilus/) were used to investigate structual variantion (SVs) in the critically endangered kākāpō (*Strigops habroptilus*). Structural variants were called using the [Delly](https://github.com/dellytools/delly), [Manta](https://github.com/Illumina/manta) and [SMOOVE](https://github.com/brentp/smoove) pipelines.

## Chapter 4: Incorporating structural variants into the conservation management of Aotearoa New Zealand's most threatened breeding bird, tara iti (*Sternula nereis davisae*)
Whole-genome short-read resequence data were aligned to the chromosomally assembled reference genome for the [common tern](https://vgp.github.io/genomeark/Sterna_hirundo/) (*Sterna hirundo*) to explore relative levels of population diversity and differentiation for two fairy tern taxa.

## Chapter 5: The promise of genome graphs as a conservation genomics resource for the critically endangered kākāpō (*Strigops habroptilus*)
Long-read Oxford Nanopore sequence data was used for *de novo* genome assembly. Genomes were assembled using [FLYE](https://github.com/fenderglass/Flye) v2.8.3. Graphs were constructed using [pggb](https://github.com/pangenome/pggb) v0.2.0, and [odgi](https://github.com/pangenome/odgi) v0.6.3.

### Publications
Publications associated with this thesis include:
Published Perspective in the *Molecular Ecology* Special Issue: [Whole-genome sequencing in molecular ecology 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16141).
