#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=80000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL

module load nixpkgs/16.09  intel/2018.3
module load samtools/1.9
module load picard/2.18.9

java -jar $EBROOTPICARD/picard.jar FixVcfHeader \
I=/home/biomatt/scratch/walleye_dna/populations_YP_out/populations.snps.vcf \
O=/home/biomatt/scratch/walleye_dna/populations_YP_out/header_populations.snps.vcf.gz \
R=/home/biomatt/scratch/yellow_perch_genome/GCF_004354835.1_PFLA_1.0_genomic.fasta