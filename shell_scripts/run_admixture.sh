#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=0
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --array=1-10

module load nixpkgs/16.09
module load admixture/1.3.0

cd /home/biomatt/scratch/walleye_dna/admixture

admixture /home/biomatt/scratch/walleye_dna/populations_YP_filtered/pruned_mac3_gt0.7_ind0.7.bed -B1000 -j48 $SLURM_ARRAY_TASK_ID
