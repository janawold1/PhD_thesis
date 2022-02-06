#!/bin/bash -e

#SBATCH -J subsampling-kaki-array
#SBATCH -A moragar
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=inv-blade-g8
#SBATCH --array=1-696%64

task=`cat subsampling_task_list.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
srun $task;

