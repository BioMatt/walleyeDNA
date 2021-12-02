#!/bin/bash
#SBATCH --account=def-kmj477
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --time=16:00:00
#SBATCH --mem=50000M
#SBATCH --array=1-4

module load nixpkgs/16.09 gcc/7.3.0
module load stacks/2.3e

list=/home/biomatt/scratch/walleye_dna/raw_wall_DNA.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${list}"
str=$($string)

var=$(echo $str | awk -F"\t" '{print $1, $2, $3}')
set -- $var

forward=$1
reverse=$2
barcodes=$3

echo ${forward}
echo ${reverse}
echo ${barcodes}

process_radtags -1 ${forward} -2 ${reverse} -b ${barcodes} -e pstI -i gzfastq -c -q -r --filter_illumina --bestrad -t 140 -o /home/biomatt/scratch/walleye_dna/processed_tags