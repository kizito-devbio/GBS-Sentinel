/*
 * Background Genome Selection — collects sample genomes and reference for context.
 *
 * Input:  collected genome FASTA files
 * Output: background genome collection, step log
 */

process BACKGROUND_SELECTION {
    tag "background_selection"
    label 'med_compute'

    conda 'conda-forge::wget=1.21.4'
    container 'kizitodevbio/strepto-pipeline:v1.2.0'

    publishDir "${params.outdir}/Reports", mode: 'copy', pattern: "background_genomes/**"
    publishDir "${params.outdir}/Logs", mode: 'copy', pattern: "background_selection_log.txt"

    input:
    path genomes

    output:
    path "background_genomes/*.fna", emit: background_fasta
    path "background_selection_log.txt", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="background_selection_log.txt"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "========================================"
        echo "Stage:       Background Genome Selection"
        echo "Started:     \$TS"
        echo "Input:       ${genomes.size()} genome file(s)"
        echo "Output dir:  ${params.outdir}/Reports"
        echo "Reference:   ${params.background_ref_url}"
        echo "----------------------------------------"
    } > "\$LOG"

    mkdir -p background_genomes

    for genome_file in ${genomes}; do
        cp "\$genome_file" background_genomes/\$(basename "\$genome_file")
        echo "Added:       \$(basename "\$genome_file")" >> "\$LOG"
    done

    REF_GENOME="background_genomes/S_agalactiae_ref.fna"
    if [ ! -f "\$REF_GENOME" ]; then
        echo "Downloading reference genome..." >> "\$LOG"
        wget -q "${params.background_ref_url}" -O "\$REF_GENOME.gz"
        gunzip -f "\$REF_GENOME.gz"
    else
        echo "Reference genome already present." >> "\$LOG"
    fi

    ELAPSED=\$(( \$(date +%s) - START ))
    echo "----------------------------------------" >> "\$LOG"
    echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    echo "Elapsed:     \${ELAPSED}s" >> "\$LOG"
    echo "Status:      SUCCESS" >> "\$LOG"
    echo "========================================" >> "\$LOG"
    """
}
