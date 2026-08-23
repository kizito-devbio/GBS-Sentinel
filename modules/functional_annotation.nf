/*
 * Functional Annotation — Prokka gene prediction and Abricate AMR screening (CARD).
 *
 * Input:
 *   tuple(sample, genome FASTA)
 *
 * Output:
 *   Prokka annotation directory
 *   AMR TSV
 *   annotation log
 */


process FUNCTIONAL_ANNOTATION {

    tag { sample }

    label 'med_compute'


    container 'kizitodevbio/strepto-pipeline:v1.2.0'


    publishDir "${params.outdir}/Annotation",
        mode: 'copy',
        pattern: "prokka_*"


    publishDir "${params.outdir}/AMR",
        mode: 'copy',
        pattern: "amr_results/*_amr.tsv"


    publishDir "${params.outdir}/Logs",
        mode: 'copy',
        pattern: "*_annotation.log"



    input:

    tuple val(sample), path(fasta_file)



    output:

    tuple val(sample),
          path("amr_results/${sample}_amr.tsv"),
          path("prokka_${sample}"),
          emit: annotation_results


    path "${sample}_annotation.log",
         emit: log



    script:

    """

    set -euo pipefail


    LOG="${sample}_annotation.log"

    START=\$(date +%s)

    TS=\$(date +"%Y-%m-%d %H:%M:%S")


    PROKKA_OUT="prokka_${sample}"

    AMR_OUT="amr_results"



    {

        echo "========================================"

        echo "Stage:       Functional Annotation"

        echo "Sample:      ${sample}"

        echo "Started:     \$TS"

        echo "Input:       ${fasta_file}"

        echo "Software:    Prokka \$(prokka --version 2>&1 | head -1)"

        echo "Software:    Abricate \$(abricate --version 2>&1 | head -1)"

        echo "Parameters:  kingdom=Bacteria genus=Streptococcus database=CARD"

        echo "----------------------------------------"


    } > "\$LOG"


mkdir -p "\$PROKKA_OUT" "\$AMR_OUT"

set +e

prokka \
    --outdir "\$PROKKA_OUT" \
    --prefix "${sample}" \
    --kingdom Bacteria \
    --genus Streptococcus \
    --force \
    --cpus ${task.cpus} \
    "${fasta_file}" >> "\$LOG" 2>&1

PROKKA_EXIT=\$?

set -e

if [ \$PROKKA_EXIT -ne 0 ]; then
    echo "WARNING: Prokka exited with code \$PROKKA_EXIT (checking for usable output anyway)" >> "\$LOG"
fi


PROKKA_FNA="\$PROKKA_OUT/${sample}.fna"


if [ -f "\$PROKKA_FNA" ]; then

    abricate \
        --db card \
        "\$PROKKA_FNA" \
        > "\$AMR_OUT/${sample}_amr.tsv" \
        2>> "\$LOG"


    AMR_COUNT=\$(grep -vc "^#" "\$AMR_OUT/${sample}_amr.tsv" || true)

    echo "AMR genes: \$AMR_COUNT" >> "\$LOG"


else

    echo "ERROR: Prokka did not produce ${sample}.fna" >> "\$LOG"

    exit 1

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


