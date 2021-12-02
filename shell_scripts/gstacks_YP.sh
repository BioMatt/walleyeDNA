#!/bin/bash
#SBATCH --time=07:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=0
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

module load nixpkgs/16.09 gcc/7.3.0
module load stacks/2.3e

gstacks -I /home/biomatt/scratch/walleye_dna/bwa_filtered_reads -O /scratch/biomatt/walleye_dna/gstacks_filtered_out -M /home/biomatt/scratch/walleye_dna/combined_popmap2.txt -t $SLURM_CPUS_PER_TASK