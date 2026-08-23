process DOWNLOAD_HUMAN_REF {
    tag "download_human_ref"
    label 'med_compute'

    // REMOVED the 'conda-forge::' prefixes to let Conda find tools in their natural channels
    conda "curl=8.5.0 bwa=0.7.17 grep=3.11"
    container 'kizitodevbio/strepto-pipeline:latest'

    publishDir params.human_ref_dir, mode: 'copy'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'human_setup_log.txt'

    input:
    val human_accession

    output:
    path "human_genome.fna", emit: human_ref
    path "human_genome.fna.*", emit: index_files
    path "human_setup_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    timestamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\$timestamp] Starting Human Reference Setup." > human_setup_log.txt

    mkdir -p human_ref

    if [ -f "${params.human_ref_dir}/human_genome.fna" ]; then
        echo "[\$timestamp] Local reference found. Copying..." >> human_setup_log.txt
        cp "${params.human_ref_dir}/human_genome.fna" human_genome.fna
        cp ${params.human_ref_dir}/human_genome.fna.* . 2>/dev/null || true
    else
        echo "[\$timestamp] Reference not found. Downloading FULL GRCh38..." >> human_setup_log.txt
        curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz" -o human.fna.gz
        gunzip -c human.fna.gz > human_genome.fna
        rm human.fna.gz
    fi

    if [ ! -f human_genome.fna.amb ]; then
        echo "[\$timestamp] Generating BWA Index..." >> human_setup_log.txt
        bwa index human_genome.fna
    else
        echo "[\$timestamp] BWA Index exists. Skipping." >> human_setup_log.txt
    fi
    """
}
