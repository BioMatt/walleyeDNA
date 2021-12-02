#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=20000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-346

module load nixpkgs/16.09 intel/2018.3
module load bwa/0.7.17
module load samtools/1.10

list=/home/biomatt/scratch/walleye_dna/processed_reads.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${list}" 
str=$($string) 

var=$(echo $str | awk -F"\t" '{print $1, $2}') 
set -- $var 
forward=$1 
reverse=$2

echo ${forward}
echo ${reverse}

wall_id=$(basename ${forward} .1.fq.gz)

echo ${wall_id}

samtools sort -@ 7 /home/biomatt/scratch/walleye_dna/bwa_aligned_reads/${wall_id}.bam > /home/biomatt/scratch/walleye_dna/bwa_sorted_reads/${wall_id}.bam