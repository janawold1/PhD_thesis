# Long-read SV discovery and genotyping
Here are some brief notes of how I conducted long-read SV discovery with Sniffles and refined calles with Jasmine. 

## Genotyping with Paragraph
Genotyping SVs was challenging due to formatting SV calls to be compatible with Paragraph. First, Transversions and unplaced scaffolds were excluded while the VCF was normalised.
```
bcftools view -e 'SVTYPE="TRAV"' -O v -o sniffles_merged_autosomes_noTRA.vcf.gz sniffles_merged_autosomes.vcf
bcftools norm -T ../metadata/kakapo_chromosome_scaffolds.bed --threads 24 --check-ref s -D -f /kakapo-data/References/kakapo_full_ref.fa -O v -o candidates.vcf sniffles_merged_autosomes_noTRA.vcf.gz
```

Once this was complete, all calls within 125 bp of chromosomal scaffold start and ends were exluded (n = 2).
```
bcftools view -T ^bad_sites.tsv -O v -o candidates_noBadSites.vcf candidates.vcf
```
