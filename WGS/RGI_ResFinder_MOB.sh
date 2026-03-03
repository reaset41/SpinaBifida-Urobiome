
#!/usr/bin/env bash

# the following code chunks were used for running RGI with the CARD database, ResFinder, and MOB suite

# RGI 

rgi load --card_json card.json --local

for file in *.fasta; do
    sample=$(basename "$file" .fasta)
    outdir="./RGI/$sample"

    # Skip if output directory already exists
    if [ -d "$outdir" ]; then
        echo "Skipping $sample — output directory already exists."
        continue
    fi
    
    echo "Running RGI on $sample"
    rgi main --input_sequence "$file" --output_file ./RGI/"$outdir" --local --clean
done

# ResFinder

for file in *.fasta; do
    sample=$(basename "$file" .fasta)
    outdir="./resfinder/$sample"

    # Skip if output directory already exists
    if [ -d "$outdir" ]; then
        echo "Skipping $sample — output directory already exists."
        continue
    fi

    echo "Running ResFinder on $sample"
    run_resfinder.py \
        -ifa "$file" \
        -o "$outdir" \
        -db_res /home/path/to/resfinder_db \
        -acq \
        -l 0.6 \
        -t 0.9 \
        --species "Escherichia"
done


# MOB suite 

for file in *.fasta; do
sample=$(basename "$file" .fasta)
    outdir="./MOB/$sample"

    # Skip if output directory already exists
    if [ -d "$outdir" ]; then
        echo "Skipping $sample — output directory already exists."
        continue
    fi
    
    echo "Running MOB on $sample"
    mob_recon --infile "$file" --outdir "$outdir"
done