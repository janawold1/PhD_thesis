#!/bin/bash -e

#SBATCH -J subsampling-kaki-trimmed_convert-array
#SBATCH -A moragar
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --partition=inv-blade-g8
#SBATCH --array=1-2796%128

task=`cat conversion_subtrimmed_tasks.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
source activate /dataset/Bird_genomes_Canterbury2017/ztmp/samtools
srun $task;

