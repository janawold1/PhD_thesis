#!/bin/bash
########################################################################
# Script for changing the scaffold names in a given file file.
########################################################################

file=~/gdrive/Data/Kakapo/Genomic\ Analyses/MANTA/R_analyses/supported_variants/inputs/chr_counts.tsv

sed -i 's/NC_044277.2/S01/g' ${file}
sed -i 's/NC_044278.2/S02/g' ${file}
sed -i 's/NC_044302.2/S03/g' ${file}
sed -i 's/NC_044279.2/S04/g' ${file}
sed -i 's/NC_046358.1/S05/g' ${file}
sed -i 's/NC_044281.2/S06/g' ${file}
sed -i 's/NC_044282.2/S07/g' ${file}
sed -i 's/NC_044283.2/S08/g' ${file}
sed -i 's/NC_044284.2/S09/g' ${file}
sed -i 's/NC_044285.2/S10/g' ${file}
sed -i 's/NC_046359.1/S11/g' ${file}
sed -i 's/NC_046360.1/S12/g' ${file}
sed -i 's/NC_044301.2/S13/g' ${file}
sed -i 's/NC_044288.2/S14/g' ${file}
sed -i 's/NC_044289.2/S15/g' ${file}
sed -i 's/NW_022651054.1/S16/g' ${file}
sed -i 's/NC_044290.2/S17/g' ${file}
sed -i 's/NC_044291.2/S18/g' ${file}
sed -i 's/NC_044292.2/S19/g' ${file}
sed -i 's/NC_044293.2/S20/g' ${file}
sed -i 's/NC_044294.2/S21/g' ${file}
sed -i 's/NC_044295.2/S22/g' ${file}
sed -i 's/NC_044296.2/S23/g' ${file}
sed -i 's/NC_044297.2/S24/g' ${file}
sed -i 's/NC_044298.2/S25/g' ${file}
sed -i 's/NC_044299.2/S26/g' ${file}

#To rename files...
rename NC_044277.2 S01 NC_*
rename NC_044278.2 S02 NC_*
rename NC_044302.2 S03 NC_*
rename NC_044279.2 S04 NC_*
rename NC_046358.1 S05 NC_*
rename NC_044281.2 S06 NC_*
rename NC_044282.2 S07 NC_*
rename NC_044283.2 S08 NC_*
rename NC_044284.2 S09 NC_*
rename NC_044285.2 S10 NC_*
rename NC_046359.1 S11 NC_*
rename NC_046360.1 S12 NC_*
rename NC_044301.2 S13 NC_*
rename NC_044288.2 S14 NC_*
rename NC_044289.2 S15 NC_*
rename NW_022651054.1 S16 NW_*
rename NC_044290.2 S17 NC_*
rename NC_044291.2 S18 NC_*
rename NC_044292.2 S19 NC_*
rename NC_044293.2 S20 NC_*
rename NC_044294.2 S21 NC_*
rename NC_044295.2 S22 NC_*
rename NC_044296.2 S23 NC_*
rename NC_044297.2 S24 NC_*
rename NC_044298.2 S25 NC_*
rename NC_044299.2 S26 NC_*
