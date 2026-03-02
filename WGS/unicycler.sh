
#!/usr/bin/env bash

# Config
READS_DIR="reads"              
OUT_DIR="WGS_assemblies"       
THREADS=16
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

# Loop over R1 files and derive the sample name
for r1 in "${READS_DIR}"/*_R1.fastq.gz; do
  # Skip if no matches
  [[ -e "$r1" ]] || continue

  # Derive R2 path by replacing R1 pattern
  r2="${r1/_R1/_R2}"

  # Sample name stripping common suffixes
  base=$(basename "$r1")
  sample="${base%_R1.fastq.gz}"

  # Sanity check
  if [[ ! -f "$r2" ]]; then
    echo "[WARN] Missing R2 for $r1 — skipping."
    continue
  fi

  sample_out="${OUT_DIR}/${sample}"
  mkdir -p "$sample_out"

  echo "Assembling ${sample}..."
  unicycler \
    -1 "$r1" \
    -2 "$r2" \
    -o "$sample_out" \
    --threads "$THREADS" \
    --min_fasta_length 500 \
    --keep 0 \
    --mode normal \
    --verbosity 1 \
    > "${LOG_DIR}/${sample}.log"

done
