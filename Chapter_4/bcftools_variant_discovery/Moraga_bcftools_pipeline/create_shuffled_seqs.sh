#!/bin/sh

mkfifo read1
mkfifo read2
gunzip -c H01391-L1_S6_L003_R1_001.fastq.gz > read1 &
gunzip -c H01391-L1_S6_L003_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01391-L1_S6_L003_shuffled.fastq
gunzip -c H01391-L1_S6_L004_R1_001.fastq.gz > read1 &
gunzip -c H01391-L1_S6_L004_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01391-L1_S6_L004_shuffled.fastq
gunzip -c H01392-L1_S7_L005_R1_001.fastq.gz > read1 &
gunzip -c H01392-L1_S7_L005_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01392-L1_S7_L005_shuffled.fastq
gunzip -c H01392-L1_S7_L006_R1_001.fastq.gz > read1 &
gunzip -c H01392-L1_S7_L006_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01392-L1_S7_L006_shuffled.fastq
gunzip -c H01393-L1_S8_L005_R1_001.fastq.gz > read1 &
gunzip -c H01393-L1_S8_L005_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01393-L1_S8_L005_shuffled.fastq
gunzip -c H01393-L1_S8_L006_R1_001.fastq.gz > read1 &
gunzip -c H01393-L1_S8_L006_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01393-L1_S8_L006_shuffled.fastq
gunzip -c H01394-L1_S9_L005_R1_001.fastq.gz > read1 &
gunzip -c H01394-L1_S9_L005_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01394-L1_S9_L005_shuffled.fastq
gunzip -c H01394-L1_S9_L006_R1_001.fastq.gz > read1 &
gunzip -c H01394-L1_S9_L006_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01394-L1_S9_L006_shuffled.fastq
gunzip -c H01395-L1_S10_L007_R1_001.fastq.gz > read1 &
gunzip -c H01395-L1_S10_L007_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01395-L1_S10_L007_shuffled.fastq
gunzip -c H01395-L1_S10_L008_R1_001.fastq.gz > read1 &
gunzip -c H01395-L1_S10_L008_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01395-L1_S10_L008_shuffled.fastq
gunzip -c H01396-L1_S11_L007_R1_001.fastq.gz > read1 &
gunzip -c H01396-L1_S11_L007_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01396-L1_S11_L007_shuffled.fastq
gunzip -c H01396-L1_S11_L008_R1_001.fastq.gz > read1 &
gunzip -c H01396-L1_S11_L008_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01396-L1_S11_L008_shuffled.fastq
gunzip -c H01397-L1_S12_L007_R1_001.fastq.gz > read1 &
gunzip -c H01397-L1_S12_L007_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01397-L1_S12_L007_shuffled.fastq
gunzip -c H01397-L1_S12_L008_R1_001.fastq.gz > read1 &
gunzip -c H01397-L1_S12_L008_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01397-L1_S12_L008_shuffled.fastq
gunzip -c H01398-L1_S1_L001_R1_001.fastq.gz > read1 &
gunzip -c H01398-L1_S1_L001_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01398-L1_S1_L001_shuffled.fastq
gunzip -c H01398-L1_S1_L002_R1_001.fastq.gz > read1 &
gunzip -c H01398-L1_S1_L002_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01398-L1_S1_L002_shuffled.fastq
gunzip -c H01399-L1_S2_L001_R1_001.fastq.gz > read1 &
gunzip -c H01399-L1_S2_L001_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01399-L1_S2_L001_shuffled.fastq
gunzip -c H01399-L1_S2_L002_R1_001.fastq.gz > read1 &
gunzip -c H01399-L1_S2_L002_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01399-L1_S2_L002_shuffled.fastq
gunzip -c H01400-L1_S3_L001_R1_001.fastq.gz > read1 &
gunzip -c H01400-L1_S3_L001_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01400-L1_S3_L001_shuffled.fastq
gunzip -c H01400-L1_S3_L002_R1_001.fastq.gz > read1 &
gunzip -c H01400-L1_S3_L002_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01400-L1_S3_L002_shuffled.fastq
gunzip -c H01401-L1_S4_L003_R1_001.fastq.gz > read1 &
gunzip -c H01401-L1_S4_L003_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01401-L1_S4_L003_shuffled.fastq
gunzip -c H01402-L1_S5_L003_R1_001.fastq.gz > read1 &
gunzip -c H01402-L1_S5_L003_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01402-L1_S5_L003_shuffled.fastq
gunzip -c H01403-L1_S6_L003_R1_001.fastq.gz > read1 &
gunzip -c H01403-L1_S6_L003_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01403-L1_S6_L003_shuffled.fastq
gunzip -c H01404-L1_S7_L005_R1_001.fastq.gz > read1 &
gunzip -c H01404-L1_S7_L005_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01404-L1_S7_L005_shuffled.fastq
gunzip -c H01404-L1_S7_L006_R1_001.fastq.gz > read1 &
gunzip -c H01404-L1_S7_L006_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01404-L1_S7_L006_shuffled.fastq
gunzip -c H01405-L1_S8_L005_R1_001.fastq.gz > read1 &
gunzip -c H01405-L1_S8_L005_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01405-L1_S8_L005_shuffled.fastq
gunzip -c H01405-L1_S8_L006_R1_001.fastq.gz > read1 &
gunzip -c H01405-L1_S8_L006_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01405-L1_S8_L006_shuffled.fastq
gunzip -c H01406-L1_S9_L005_R1_001.fastq.gz > read1 &
gunzip -c H01406-L1_S9_L005_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01406-L1_S9_L005_shuffled.fastq
gunzip -c H01406-L1_S9_L006_R1_001.fastq.gz > read1 &
gunzip -c H01406-L1_S9_L006_R2_001.fastq.gz > read2 &
shuffleSequences_fastq.pl read1 read2 H01406-L1_S9_L006_shuffled.fastq
