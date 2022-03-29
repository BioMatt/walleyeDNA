#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --account=def-kmj477
#SBATCH --mem=1000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-346

module load nixpkgs/16.09 intel/2018.3
module load samtools/1.10

list=/home/biomatt/scratch/walleye_dna/wall_bamlist.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${list}" 
str=$($string) 

var=$(echo $str | awk -F"\t" '{print $1}') 
set -- $var 

bam_file=$1 

echo ${bam_file}

samtools flagstat ${bam_file}
