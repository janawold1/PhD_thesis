#!/bin/bash -e

#SBATCH -J subsampling-kaki-trimmed-array
#SBATCH -A moragar
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G
#SBATCH --partition=inv-blade-g8
#SBATCH --array=1-696%32

task=`cat alignment_subtrimmed_tasks.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
source activate /dataset/Bird_genomes_Canterbury2017/ztmp/nextgenmap
srun $task;

