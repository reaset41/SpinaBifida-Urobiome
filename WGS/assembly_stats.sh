#!/usr/bin/env bash

echo -e "File\tContigs\tGenomeLength\tN50" > assembly_stats.tsv

# Make sure files exist
set -- *.fasta
[ -e "$1" ] || { echo "No .fasta files found."; exit 1; }

for f in *.fasta; do

  # Contig count
  contigs=$(grep -c "^>" "$f")

  # Extract contig lengths and store them
  lengths=$(cat "$f" | \
    awk '
      /^>/ { if (seqlen) print seqlen; seqlen=0; next }
      { seqlen += length($0) }
      END { print seqlen }
    ')

  # Compute genome length
  genome_length=$(echo "$lengths" | awk '{sum+=$1} END{print sum}')

  # Compute N50
  n50=$(echo "$lengths" | \
    sort -nr | \
    awk '
      {
        sum += $1
        lengths[NR] = $1
      }
      END {
        half = sum / 2
        for (i=1; i<=NR; i++) {
          acc += lengths[i]
          if (acc >= half) {
            print lengths[i]
            break
          }
        }
      }
    ')

  echo -e "${f}\t${contigs}\t${genome_length}\t${n50}" >> assembly_stats.tsv

done