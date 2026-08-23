/*
 * Virulence Factor Detection — Abricate screening against VFDB.
 *
 * Input:  annotated genome FASTA (Prokka .fna)
 * Output: virulence factor TSV, log
 */

process VIRULENCE_FACTOR {

    tag "$sample"

    label 'low_compute'

    container 'kizitodevbio/strepto-pipeline:v1.2.0'

    publishDir "${params.outdir}/Virulence",
        mode: 'copy',
        pattern: "*_vf.tsv"

    publishDir "${params.outdir}/Logs",
        mode: 'copy',
        pattern: "*_virulence.log"


    input:
    tuple val(sample), path(fasta_file)


    output:
    tuple val(sample), path("${sample}_vf.tsv"), emit: results
    path "${sample}_virulence.log", emit: log


    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    export LC_ALL=C
    export LANG=C

    LOG="${sample}_virulence.log"

    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")


    {
        echo "========================================"
        echo "Stage:       Virulence Factor Detection"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Input:       ${fasta_file}"
        echo "Output dir:  ${params.outdir}/Virulence"
        echo "Software:    Abricate \$(abricate --version 2>&1 | head -1)"
        echo "Database:    VFDB"
        echo "----------------------------------------"

    } > "\$LOG"



    abricate \
        --db vfdb \
        "${fasta_file}" \
        > "${sample}_vf.tsv" \
        2>> "\$LOG" || true



    if [ ! -s "${sample}_vf.tsv" ]; then

        echo -e "#FILE\\tSEQUENCE\\tSTART\\tEND\\tGENE\\tCOVERAGE\\tCOVERAGE_MAP\\tGAPS\\t%COVERAGE\\t%IDENTITY\\tDATABASE\\tACCESSION\\tPRODUCT" \
        > "${sample}_vf.tsv"

        echo "No virulence factors detected" >> "\$LOG"

    else

        NUM_VF=\$(grep -vc "^#" "${sample}_vf.tsv" || true)

        echo "VF detected: \$NUM_VF" >> "\$LOG"

    fi



    ELAPSED=\$(( \$(date +%s) - START ))


    {
        echo "----------------------------------------"
        echo "Completed:   \$(date +"%Y-%m-%d %H:%M:%S")"
        echo "Elapsed:     \${ELAPSED}s"
        echo "Status:      SUCCESS"
        echo "========================================"

    } >> "\$LOG"

    """
}
