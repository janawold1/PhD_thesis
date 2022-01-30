# SubSampler_SNPcaller
A series of scripts and tool to subsample illumina data and do SNP calling on them (SNP calling based off the software carpentry pipeline).

Steps:

1. create_shuffled_seqs.sh : Create a set of shuffled fastq files from the two ends of each sequenced sample.

2. batch_subsample.sh + subsampling_task_list.txt : This is a SLURM task list that will create the subsampled fastq files from 3x to 10x (note that we won't use them all at the end, but it's handy to have it all in one list. Only used 3x, 4x, 5x, and 9x).

3. split_paired.pl : This is a perl script that will split the shuffled files back into separate files. Technically, we could re-write the subsampler to take paired files separately.

4. batch_split.sh + split_task_list.txt : SLURM task list that will split all the files we created before into two pair-end files.

5. batch_trim_subsampled.sh + trim_task_list.txt : SLURM task list that will trim all the files created before. Note that we could kinda skip some of the jobs if the subsampling was done by using a fixed length for the subsampling script (modification that I might add later to save on tasks), but skewer is so quick that any decent cluster can chew through these tasks over a teabreak really. And this way it's all independent anyway.

6. batch_subtrimmed_alignment.sh + alignment_subtrimmed_tasks.txt : SLURM task list that will map the reads onto the reference genome. Again, all mappings done independently (as opposed to subsampling the alignments).

7. batch_bcf_call.sh + conversion_subtrimmed_tasks.txt : SLURM task list that does the finall SNP calling. It will convert the format of the sam files into bam, index them, then use mpileup to create the raw BCF files. Further filtering in these BCF files is recommended, but it'll be to user's taste. I use BCF instead of VCF just to save one conversion step when loading them into bcftools, but it's essentially the same.

The bcf files have been uploaded to the Dropbox (link via e-mail).

Additinal script(s):

1. subsample_fastq_dir.pl: A perl script that will take a directory full of fastq files and batch-subsample them to a certain percentage based on used input.
2. split_bamfiles_tasks.pl: This script will take a list of BAM files and a genome, and split those bamfiles into N user-defined non-overlapping alignments. The idea for this is to create subsets of BAM files that can be independently run through mpileup before doing SNP calling with bcftools. Because mpileup doesn't run in parallel, this will do the job. The output is a list of samtools tasks that can be sent to your favourite queue manager.
