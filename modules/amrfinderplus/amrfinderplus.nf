/*
 * AMRFinderPlus — standardized antimicrobial-resistance / resistome detection.
 *
 * Primary resistome module for GBS-Sentinel.
 *
 * Organism:
 *   Streptococcus agalactiae
 *
 * Requirements:
 *   - AMRFinderPlus executable
 *   - prepared AMRFinderPlus database containing AMRProt.fa
 *     and its BLAST index files
 *
 * The module deliberately validates the database before execution.
 * A missing or unusable database causes the process to fail rather than
 * producing an empty TSV that could be mistaken for "no AMR detected".
 *
 * Input:
 *   tuple(sample, fasta_file)
 *
 * Output:
 *   tuple(sample, sample_amrfinder.tsv)
 *   sample_amrfinder.log
 */

process AMRFINDERPLUS {

    tag "$sample"
    label 'amrfinder_high_memory'


    container 'kizitodevbio/strepto-pipeline:v1.2.0'

    publishDir "${params.outdir}/AMR",
        mode: 'copy',
        pattern: "*_amrfinder.tsv"

    publishDir "${params.outdir}/Logs",
        mode: 'copy',
        pattern: "*_amrfinder.log"

    input:

    tuple val(sample), path(fasta_file)

    output:

    tuple val(sample), path("${sample}_amrfinder.tsv"), emit: results
    path "${sample}_amrfinder.log", emit: log

    script:

    """
    #!/usr/bin/env bash

    set -euo pipefail

    LOG="${sample}_amrfinder.log"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    AMRFINDER_DB_ENV="\${AMRFINDER_DB:-}"

    {
        echo "========================================"
        echo "Stage:       AMRFinderPlus (Resistome)"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Input:       ${fasta_file}"
        echo "Software:    \$(amrfinder --version 2>&1 | head -1 || echo 'amrfinder unavailable')"
        echo "========================================"
    } > "\$LOG"


    # ============================================================
    # 1. Validate AMRFinderPlus executable
    # ============================================================

    if ! command -v amrfinder >/dev/null 2>&1; then

        echo "FATAL: AMRFinderPlus executable not found in PATH." \
            >> "\$LOG"

        exit 1

    fi


    # ============================================================
    # 2. Locate a VALID prepared AMRFinderPlus database
    #
    # A directory is valid only if AMRProt.fa.phr exists.
    #
    # This prevents the parent directory
    # /opt/amrfinder/data
    # from being selected when the real database is inside
    # a dated directory such as:
    #
    # /opt/amrfinder/data/2026-08-07.1/
    # ============================================================

    DB_DIR=""


    if [ -n "\$AMRFINDER_DB_ENV" ] && \
       [ -f "\$AMRFINDER_DB_ENV/AMRProt.fa.phr" ]; then

        DB_DIR="\$AMRFINDER_DB_ENV"

    fi


    if [ -z "\$DB_DIR" ] && \
       [ -f "/opt/amrfinder/data/latest/AMRProt.fa.phr" ]; then

        DB_DIR="/opt/amrfinder/data/latest"

    fi


    if [ -z "\$DB_DIR" ] && \
       [ -f "/opt/amrfinder/data/AMRProt.fa.phr" ]; then

        DB_DIR="/opt/amrfinder/data"

    fi


    if [ -z "\$DB_DIR" ]; then

        while IFS= read -r candidate; do

            if [ -f "\$candidate/AMRProt.fa.phr" ]; then

                DB_DIR="\$candidate"
                break

            fi

        done < <(
            find /opt/amrfinder/data \
                -mindepth 1 \
                -maxdepth 2 \
                -type d \
                -print 2>/dev/null
        )

    fi


    # ============================================================
    # 3. Fail explicitly if the database is unavailable
    # ============================================================

    if [ -z "\$DB_DIR" ]; then

        echo "FATAL: No prepared AMRFinderPlus database was found." \
            >> "\$LOG"

        echo "Expected AMRProt.fa.phr under one of:" >> "\$LOG"
        echo "  /opt/amrfinder/data/latest" >> "\$LOG"
        echo "  /opt/amrfinder/data" >> "\$LOG"
        echo "  /opt/amrfinder/data/<versioned-directory>" >> "\$LOG"

        echo "Available AMRFinderPlus data:" >> "\$LOG"

        find /opt/amrfinder/data \
            -maxdepth 3 \
            -type f \
            -name "AMRProt.fa*" \
            -print 2>/dev/null \
            >> "\$LOG" || true

        exit 1

    fi


    echo "Using database: \$DB_DIR" >> "\$LOG"


    # ============================================================
    # 4. Validate the critical AMRFinderPlus database files
    # ============================================================

    if [ ! -s "\$DB_DIR/AMRProt.fa" ]; then

        echo "FATAL: AMRProt.fa is missing from \$DB_DIR" \
            >> "\$LOG"

        exit 1

    fi


    if [ ! -s "\$DB_DIR/AMRProt.fa.phr" ]; then

        echo "FATAL: AMRProt.fa.phr is missing from \$DB_DIR" \
            >> "\$LOG"

        exit 1

    fi


    if [ ! -s "\$DB_DIR/AMRProt.fa.pin" ]; then

        echo "FATAL: AMRProt.fa.pin is missing from \$DB_DIR" \
            >> "\$LOG"

        exit 1

    fi


    if [ ! -s "\$DB_DIR/AMRProt.fa.psq" ]; then

        echo "FATAL: AMRProt.fa.psq is missing from \$DB_DIR" \
            >> "\$LOG"

        exit 1

    fi


    echo "AMRFinderPlus database validation: OK" >> "\$LOG"


    # ============================================================
    # 5. Run AMRFinderPlus with organism-specific detection
    # ============================================================

    rm -f "${sample}_amrfinder.tsv"


    echo "Running organism-specific AMRFinderPlus analysis..." \
        >> "\$LOG"


    set +e

    amrfinder \
        -n "${fasta_file}" \
        -o "${sample}_amrfinder.tsv" \
        --database "\$DB_DIR" \
        --organism Streptococcus_agalactiae \
        --threads ${task.cpus} \
        >> "\$LOG" 2>&1

    AMR_EXIT=\$?

    set -e


    # ============================================================
    # 6. Handle AMRFinderPlus failure
    # ============================================================

    if [ "\$AMR_EXIT" -ne 0 ]; then

        echo "FATAL: AMRFinderPlus failed with exit code \$AMR_EXIT" \
            >> "\$LOG"

        if [ -f "${sample}_amrfinder.tsv" ]; then

            echo "AMRFinderPlus output file was created but the analysis failed." \
                >> "\$LOG"

        fi

        exit "\$AMR_EXIT"

    fi


    # ============================================================
    # 7. Require an output file
    # ============================================================

    if [ ! -s "${sample}_amrfinder.tsv" ]; then

        echo "FATAL: AMRFinderPlus completed without producing an output TSV." \
            >> "\$LOG"

        exit 1

    fi


    # ============================================================
    # 8. Count real AMRFinderPlus detections
    #
    # Header is one line.
    # Data lines are non-empty non-comment lines after the header.
    # ============================================================

    HIT_COUNT=\$(
        awk '
            NR > 1 &&
            \$0 !~ /^#/ &&
            NF > 0 {
                count++
            }
            END {
                print count + 0
            }
        ' "${sample}_amrfinder.tsv"
    )


    echo "AMR determinants reported: \$HIT_COUNT" >> "\$LOG"


    # ============================================================
    # 9. Completion
    # ============================================================

    ELAPSED=\$(( \$(date +%s) - START ))

    {
        echo "----------------------------------------"
        echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")"
        echo "Elapsed:     \${ELAPSED}s"
        echo "AMR hits:    \$HIT_COUNT"
        echo "Status:      SUCCESS"
        echo "========================================"

    } >> "\$LOG"
    """
}
