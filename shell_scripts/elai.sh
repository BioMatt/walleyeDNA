#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=40000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-100

# Not used in the final manuscript

module load gcc/7.3.0 elai/1.0
cd /home/biomatt/scratch/walleye_dna/populations_YP_filtered/elai_out_mg1

cat /home/biomatt/scratch/walleye_dna/populations_YP_filtered/chromosome_list.txt | while read j; do \
chr=$(echo ${j}); \
elai -g /home/biomatt/scratch/walleye_dna/populations_YP_filtered/elai_inputs/LkWinnipeg_${chr}.recode.geno.txt -p 10 -g /home/biomatt/scratch/walleye_dna/populations_YP_filtered/elai_inputs/LkManitoba_${chr}.recode.geno.txt -p 11 -g /home/biomatt/scratch/walleye_dna/populations_YP_filtered/elai_inputs/Kansas_${chr}.recode.geno.txt -p 12 -g /home/biomatt/scratch/walleye_dna/populations_YP_filtered/elai_inputs/admixed_${chr}.recode.geno.txt -p 1 -pos /home/biomatt/scratch/walleye_dna/populations_YP_filtered/elai_inputs/admixed_${chr}.recode.pos.txt -R $RANDOM -s 30 -o ${chr}_run$SLURM_ARRAY_TASK_ID -C 3 -c 15 -mg 1; \
done