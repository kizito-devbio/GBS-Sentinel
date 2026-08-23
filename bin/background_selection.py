#!/usr/bin/env python3
import sys
import os
import subprocess
from Bio import Entrez, SeqIO
import pandas as pd
from datetime import datetime

def main(sample, genome_fasta, taxonomy_file, output_dir, downsample_max, email):
    Entrez.email = email
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] Step: Background Selection for sample {sample}. What is done: Downloading and filtering genus backgrounds from NCBI.")
    with open(taxonomy_file, 'r') as f:
        taxonomy = f.read().strip()
    genus = taxonomy.split(' ')[0]
    print(f"[{timestamp}] Genus identified: {genus}")

    subprocess.run(['datasets', 'download', 'genome', 'taxon', genus, '--assembly-level', 'complete', '--filename', 'backgrounds.zip', '--no-progressbar'], check=True)
    subprocess.run(['unzip', '-o', 'backgrounds.zip', '-d', 'backgrounds'], check=True)
    print(f"[{timestamp}] Download and unzip complete. Number of initial backgrounds: " + str(len(os.listdir('backgrounds/ncbi_dataset/data'))))

    filtered_paths = []
    for root, dirs, files in os.walk('backgrounds/ncbi_dataset/data'):
        for file in files:
            if file.endswith('.fna'):
                path = os.path.join(root, file)
                try:
                    length = sum(len(rec.seq) for rec in SeqIO.parse(path, 'fasta'))
                    if length >= 1470000:
                        filtered_paths.append(path)
                except Exception as e:
                    print(f"[{timestamp}] Error parsing {path}: {e}", file=sys.stderr)
                    sys.exit(1)
    print(f"[{timestamp}] Number of sequences after length filter (>=70% of GBS genome): {len(filtered_paths)}")

    catalog_path = 'backgrounds/ncbi_dataset/fetch_catalog.json'
    downsampled = []
    if os.path.exists(catalog_path):
        try:
            metadata = pd.read_json(catalog_path, lines=True)
            num_initial = len(metadata)
            print(f"[{timestamp}] Initial number of organisms in metadata: {num_initial}")
            metadata['collection_date'] = pd.to_datetime(metadata['collection_date'], errors='coerce')
            metadata['bin'] = metadata['isolation_source'].fillna('unknown') + '_' + metadata['collection_date'].dt.year.astype(str).fillna('unknown') + '_' + metadata['geographical_location'].fillna('unknown')
            downsampled_df = metadata.groupby('bin').head(downsample_max)
            downsampled = downsampled_df['accession'].tolist()
            num_passed = len(downsampled)
            print(f"[{timestamp}] Downsample threshold: {downsample_max} per bin. Number of organisms that passed downsample: {num_passed}")
        except Exception as e:
            print(f"[{timestamp}] Error processing metadata: {e}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    for acc in downsampled:
        src_dir = f"backgrounds/ncbi_dataset/data/{acc}"
        if os.path.exists(src_dir):
            fasta = os.path.join(src_dir, f"{acc}.fna")
            if os.path.exists(fasta):
                os.symlink(fasta, os.path.join(output_dir, f"{acc}.fasta"))
                print(f"[{timestamp}] Selected background organism: {acc}")
    print(f"[{timestamp}] Background selection complete for genus {genus}. Final number of downsampled backgrounds: " + str(len(downsampled)))

if __name__ == '__main__':
    if len(sys.argv) != 7:
        sys.exit("Usage: background_selection.py <sample> <genome_fasta> <taxonomy_file> <output_dir> <downsample_max> <email>")
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6])
