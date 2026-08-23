/*
 * MLST Extraction — Multi-Locus Sequence Typing using mlst (PubMLST schemes).
 *
 * Input:  genome FASTA
 * Output: MLST assignment TSV, log
 */


process MLST_EXTRACTION {

    tag "$sample"

    label 'low_compute'


    container 'kizitodevbio/strepto-pipeline:v1.2.0'


    publishDir "${params.outdir}/MLST",
        mode: 'copy',
        pattern: "*_mlst.tsv"


    publishDir "${params.outdir}/Logs",
        mode: 'copy',
        pattern: "*_mlst.log"



    input:

    tuple val(sample), path(fasta_file)



    output:

    tuple val(sample), path("${sample}_mlst.tsv"), emit: results

    path "${sample}_mlst.log", emit: log



    script:

    """
    #!/usr/bin/env bash

    set -euo pipefail

    export LC_ALL=C
    export LANG=C



    LOG="${sample}_mlst.log"

    START=\$(date +%s)

    TS=\$(date +"%Y-%m-%d %H:%M:%S")



    {

        echo "========================================"

        echo "Stage:       MLST Extraction"

        echo "Sample:      ${sample}"

        echo "Started:     \$TS"

        echo "Input:       ${fasta_file}"

        echo "Output dir:  ${params.outdir}/MLST"

        echo "Software:    mlst \$(mlst --version 2>&1 | head -1)"

        echo "Parameters:  scheme=${params.mlst_scheme}"

        echo "----------------------------------------"


    } > "\$LOG"



    mlst \
        --scheme ${params.mlst_scheme} \
        "${fasta_file}" \
        > "${sample}_mlst.tsv" \
        2>> "\$LOG" || true



    if [ ! -s "${sample}_mlst.tsv" ]; then

        echo -e "FILE\\tSCHEME\\tST\\tALLELES" \
        > "${sample}_mlst.tsv"

        echo "WARNING: No MLST assignment detected" >> "\$LOG"

    else

        echo "MLST assignment completed" >> "\$LOG"

    fi



    {

        echo "MLST result:"

        cat "${sample}_mlst.tsv"

    } >> "\$LOG"



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

