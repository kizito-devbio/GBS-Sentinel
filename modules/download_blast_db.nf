/*
 * BLAST Database Setup — download and index Streptococcus reference for taxonomy.
 *
 * Output: formatted BLAST nucleotide database, setup log
 */

process DOWNLOAD_BLAST_DB {
    label 'low_compute'
    tag "download_blast_db"

    conda "${projectDir}/envs/blast_db.yml"
    container 'kizitodevbio/strepto-pipeline:v1.2.0'

    publishDir params.blast_db_dir, mode: 'copy'
    publishDir "${params.outdir}/Logs", mode: 'copy', pattern: 'blast_setup_log.txt'

    output:
    path "blast_db", emit: blast_db
    path "blast_setup_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="blast_setup_log.txt"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       BLAST Database Setup"
        echo "Started:     \$TS"
        echo "Output dir:  ${params.blast_db_dir}"
        echo "Software:    BLAST+ \$(makeblastdb -version 2>&1 | head -1)"
        echo "----------------------------------------"
    } > "\$LOG"

    mkdir -p blast_db
    LOCAL_REF="${projectDir}/references/streptococcus_genus_db/blast_db/strep_ref.fna"

    if [ -f "\$LOCAL_REF" ]; then
        echo "Using cached reference." >> "\$LOG"
        cp "\$LOCAL_REF" blast_db/strep_ref.fna
    else
        echo "Downloading S. agalactiae NEM316 reference..." >> "\$LOG"
        curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/265/GCF_000007265.1_ASM726v1/GCF_000007265.1_ASM726v1_genomic.fna.gz" -o strep.fna.gz
        gunzip -c strep.fna.gz > blast_db/strep_ref.fna
        rm -f strep.fna.gz
    fi

    makeblastdb -in blast_db/strep_ref.fna -dbtype nucl -out blast_db/strep_db >> "\$LOG" 2>&1

    ELAPSED=\$(( \$(date +%s) - START ))
    echo "----------------------------------------" >> "\$LOG"
    echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    echo "Elapsed:     \${ELAPSED}s" >> "\$LOG"
    echo "Status:      SUCCESS" >> "\$LOG"
    echo "========================================" >> "\$LOG"
    """
}
