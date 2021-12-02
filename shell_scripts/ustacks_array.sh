#!/bin/bash
#SBATCH --account=def-kmj477
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --time=02:00:00
#SBATCH --mem=15000M
#SBATCH --array=1-384

module load nixpkgs/16.09 gcc/7.3.0
module load stacks/2.3e

prefix="/home/biomatt/scratch/walleye_dna/processed_tags/"
output="/home/biomatt/scratch/walleye_dna/ustacks_out_array"

list=/home/biomatt/scratch/walleye_dna/combined_popmap.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${list}"
str=$($string)

var=$(echo $str | awk -F"\t" '{print $1, $2}')
set -- $var

sample=$1
pop=$2

echo ${sample}
echo ${pop}

ustacks -t gzfastq \
-f "$prefix$sample$f.1.fq.gz" \
-o "$output" \
-i $SLURM_ARRAY_TASK_ID \
--disable-gapped -m 3 -M 3 -H -p 32 \
--max_locus_stacks 4 \
--model_type bounded \
--bound_high 0.05