#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=40000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=4

module load nixpkgs/16.09 gcc/7.3.0
module load stacks/2.3e

populations --in_path /home/biomatt/scratch/walleye_dna/gstacks_filtered_out/ --popmap /home/biomatt/scratch/walleye_dna/combined_popmap2.txt --min-samples-per-pop 0.6 --write-random-snp --vcf --plink --threads $SLURM_CPUS_PER_TASK --ordered-export --out_path /home/biomatt/scratch/walleye_dna/populations_YP_filtered