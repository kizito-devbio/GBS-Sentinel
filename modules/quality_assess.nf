/*
 * Assembly Quality Assessment — QUAST metrics for assembled genomes.
 *
 * Input:
 *   tuple(sample, assembled_fasta)
 *
 * Output:
 *   tuple(sample, fasta, n50, total_length, contigs)
 *   QUAST report
 *   log file
 */

process QUALITY_ASSESS {

    tag { sample }

    label 'med_compute'

    publishDir "${params.outdir}/Assembly",
        mode: 'copy',
        pattern: "*_quast_report.tsv"

    publishDir "${params.outdir}/Logs",
        mode: 'copy',
        pattern: "*_qa_log.txt"


    input:

    tuple val(sample), path(assembly)


    output:

    tuple val(sample),
          path(assembly),
          path("n50.txt"),
          path("total_len.txt"),
          path("contigs.txt"),
          emit: metrics

    path "*_quast_report.tsv",
         emit: report

    path "*_qa_log.txt",
         emit: log


    script:

    """

    LOG="${sample}_qa_log.txt"

    START=\$(date +%s)

    {

    echo "========================================"
    echo "Stage:       Assembly Quality Assessment"
    echo "Sample:      ${sample}"
    echo "Input:       ${assembly}"
    echo "Software:    QUAST"
    echo "========================================"

    } > \$LOG


    quast.py ${assembly} \
        -o quast_out \
        --min-contig 500 \
        --no-plots \
        --no-html \
        >> \$LOG 2>&1


    cp quast_out/report.tsv ${sample}_quast_report.tsv


    grep "^N50" ${sample}_quast_report.tsv \
        | cut -f2 > n50.txt || echo 0 > n50.txt


    grep "^Total length" ${sample}_quast_report.tsv \
        | cut -f2 > total_len.txt || echo 0 > total_len.txt


    grep "^# contigs" ${sample}_quast_report.tsv \
        | cut -f2 > contigs.txt || echo 0 > contigs.txt



    echo "N50: \$(cat n50.txt)" >> \$LOG
    echo "Total length: \$(cat total_len.txt)" >> \$LOG
    echo "Contigs: \$(cat contigs.txt)" >> \$LOG


    ELAPSED=\$(( \$(date +%s) - START ))


    echo "Elapsed: \${ELAPSED}s" >> \$LOG
    echo "Status: SUCCESS" >> \$LOG

    """
}
