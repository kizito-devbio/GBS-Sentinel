/*
 * Core Genome Extraction — pangenome analysis and core gene alignment via Panaroo.
 *
 * Input:  collected Prokka GFF annotation files (≥2 samples required)
 * Output: core gene alignment, Panaroo results directory, log
 */

process CORE_GENOME {
    tag { "core_genome" }
    label 'high_compute'

    conda 'conda-forge::python=3.9 bioconda::panaroo=1.5.0 bioconda::gclib'
    container 'kizitodevbio/strepto-pipeline:v1.2.0'

    publishDir "${params.outdir}/CoreGenome", mode: 'copy'
    publishDir "${params.outdir}/Logs", mode: 'copy', pattern: "core_genome.log"

    input:
    path gff_files

    output:
    path "core_alignment.aln", emit: alignment
    path "core_genome.log", emit: log
    path "panaroo_results", emit: all_results

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="core_genome.log"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")
    NUM_GFFS=\$(ls *.gff 2>/dev/null | wc -l || echo 0)

    {
        echo "========================================"
        echo "Stage:       Core Genome Extraction"
        echo "Started:     \$TS"
        echo "Input:       \$NUM_GFFS GFF file(s)"
        echo "Output dir:  ${params.outdir}/CoreGenome"
        echo "Software:    Panaroo \$(panaroo --version 2>&1 | head -1)"
        echo "Parameters:  clean-mode=strict, threads=${task.cpus}"
        echo "----------------------------------------"
    } > "\$LOG"

    if [ "\$NUM_GFFS" -lt 2 ]; then
        echo "ERROR: Core genome analysis requires ≥2 annotated genomes (found \$NUM_GFFS)." >> "\$LOG"
        echo "Skipping core genome step — phylogeny will not be produced." >> "\$LOG"
        echo "NO_ALIGNMENT" > core_alignment.aln
        mkdir -p panaroo_results
        exit 0
    fi

    panaroo \\
        -i *.gff \\
        -o panaroo_results \\
        --clean-mode strict \\
        -a core \\
        --threads ${task.cpus} >> "\$LOG" 2>&1

    if [ -f "panaroo_results/core_gene_alignment.aln" ]; then
        cp "panaroo_results/core_gene_alignment.aln" "core_alignment.aln"
        ALIGN_LEN=\$(grep -v '^>' core_alignment.aln | tr -d '\\n' | wc -c)
        echo "Alignment:   \$ALIGN_LEN bp" >> "\$LOG"
    else
        echo "ERROR: Panaroo did not produce core_gene_alignment.aln" >> "\$LOG"
        exit 1
    fi

    ELAPSED=\$(( \$(date +%s) - START ))
    echo "----------------------------------------" >> "\$LOG"
    echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")" >> "\$LOG"
    echo "Elapsed:     \${ELAPSED}s" >> "\$LOG"
    echo "Status:      SUCCESS" >> "\$LOG"
    echo "========================================" >> "\$LOG"
    """
}
