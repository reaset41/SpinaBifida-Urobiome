#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=breseq.txt

module load breseq 

while read -r col2 col6 col7 col8; do
    echo "Input: $col7 $col8"
    echo "Reference: $col6"
    echo "Output: $col2"

breseq -r $col6 -o $col2 $col7 $col8 -j 16

done < <(awk 'NR > 1 {print $2, $6, $7, $8}' BreSeq_Strains.txt)