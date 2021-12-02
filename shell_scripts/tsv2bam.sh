#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=0
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --constraint=skylake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

module load nixpkgs/16.09 gcc/7.3.0
module load stacks/2.3e

tsv2bam -P /home/biomatt/scratch/walleye_dna/ustacks_out_array -M /home/biomatt/scratch/walleye_dna/combined_popmap2.txt -R /home/biomatt/scratch/walleye_dna/processed_tags -t $SLURM_CPUS_PER_TASK