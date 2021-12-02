#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=20000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL

module load nixpkgs/16.09 intel/2018.3
module load bwa/0.7.17

bwa index -p /home/biomatt/scratch/yellow_perch_genome/bwa_perch_index -a bwtsw /home/biomatt/scratch/yellow_perch_genome/GCF_004354835.1_PFLA_1.0_genomic.fna