#!/bin/bash -e

#SBATCH -J trimming-kaky-array
#SBATCH -A moragar
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=inv-blade-g8
#SBATCH --array=1-29%16

task=`cat trim_task_list.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
source activate /dataset/Bird_genomes_Canterbury2017/ztmp/skewer
srun $task;

