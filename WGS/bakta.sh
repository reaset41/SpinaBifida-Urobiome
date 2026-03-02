
#!/usr/bin/env bash

THREADS=16

for file in *.fasta; do
    sample=$(basename "$file" .fasta)
    outdir="bakta/$sample"

    # Skip if output directory already exists
    if [ -d "$outdir" ]; then
        echo "Skipping $sample — output directory already exists."
        continue
    fi

    echo "Running Bakta on $sample"
	bakta --db /data/path/bakta_db/db \
	--force \
	--threads "$THREADS" \
	--output "$outdir" "$file"

done
