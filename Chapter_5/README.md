# Chapter 5: The promise of genome graphs as a consdervation genomics resource for the critically endangered kākāpō (*Strigops habroptilus*)
Here are the steps I used for *de novo* genome assembly and graph contstruction in the critically endangered kākāpō. Here, we specifically examined the graph topology of toll-like receptors found on this chromosome.  

## 1 ONT base calling and read quality statistics
Oxford Nanopore long-read (ONT) data were used in this chapter. Raw reads were base called using Guppy v6.0.1. Basecalled reads were assessed using [MinIONQC](https://github.com/roblanf/minion_qc) v1.4.1. Adaptor trimming was conducted with [porechop](https://github.com/rrwick/Porechop) v0.2.4, lambda DNA removed with [NanoLyse](https://github.com/wdecoster/nanolyse) v1.2.0 and reads filtered for minimum quality (Q-score) and length using [NanoFilt](https://github.com/wdecoster/nanofilt) v2.8.0. Overall reads were visualised with [NanoPlot](https://github.com/wdecoster/NanoPlot) v1.39.0 before and after trimming.

## 2 Coverage inferred from mapping to a linear reference
A number of the samples sequenced here were >16 years old. To explore whether ONT data reads were evenly distributed across the genome, or some regions were more represented than others due to degradation of DNA, reads were aligned to the high-quality reference assembly for [Jane's genome](https://www.ncbi.nlm.nih.gov/genome/?term=kakapo). Reads were aligned using [winnowmap](https://github.com/marbl/Winnowmap) v2.0.3 and coverage was estimated using SAMtools and [mosdepth](https://github.com/brentp/mosdepth) v0.3.3.

## 3 *De novo* genome assembly
The long-read assembler [FLYE](https://github.com/fenderglass/Flye) v2.8.3 was used for *de novo* genome assembly. Multiple polishing regimes were trialled. These included [NextPolish](https://github.com/Nextomics/NextPolish) v1.4.0, [Medaka](https://github.com/nanoporetech/medaka) v1.5.0 and [Racon](https://github.com/lbcb-sci/racon) v1.4.20. Genome statistics (e.g., N50, assembly size, L50) were estimated with [quast](https://github.com/ablab/quast) v5.0.2 and genome 'completeness' was assessed using [BUSCO](https://busco.ezlab.org/)v5.3.0 scores were used to identify the best polished outputs. 

## 4 Genome graph construction
Once the assembly was finalised, contigs aligning to autosome 7 in [Jane's genome](https://www.ncbi.nlm.nih.gov/genome/?term=kakapo) were identified using [Mash](https://github.com/marbl/Mash) v2.3 and isolated using the perl script faSomeRecords. [D-genies](http://dgenies.toulouse.inra.fr/) v1.3.0 was used to scaffold reads into a single alignment prior to graph construction with [pggb](https://github.com/pangenome/pggb) v0.2.0. Graph statistics were estimated using [odgi](https://github.com/pangenome/odgi) v0.5.1. Proportions of node and edges were visualised using [MulitQC](https://github.com/ewels/MultiQC) v1.11. Final graphs were visualised using [Bandage](https://github.com/rrwick/Bandage) v0.9.0.

## 5 Plotting of graph depth
Depth estimated by ODGI were plotted for comparison with graph topology. 