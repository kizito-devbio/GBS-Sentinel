/*
 * Taxonomic Confirmation — BLASTn against Streptococcus reference database.
 *
 * Input:
 *   tuple(sample, genome)
 *   BLAST database directory
 *
 * Output:
 *   taxonomy assignment text file
 *   taxonomy log
 *
 * Purpose:
 *   Confirms Streptococcus taxonomy using a local BLAST database.
 */

process TAXONOMY {

    tag { sample }

    label 'low_compute'

    container 'kizitodevbio/strepto-pipeline:v1.2.0'

    stageInMode 'copy'

    publishDir "${params.outdir}/Annotation",
        mode: 'copy',
        pattern: "*_taxonomy.txt"

    publishDir "${params.outdir}/Logs",
        mode: 'copy',
        pattern: "*_taxonomy.log"

    input:
    tuple val(sample), path(genome)
    path blast_db_folder

    output:
    tuple val(sample),
          path(genome),
          path("${sample}_taxonomy.txt"),
          emit: results

    path "${sample}_taxonomy.log",
         emit: log

    shell:
    '''

    set -euo pipefail

    LOG="!{sample}_taxonomy.log"

    START=$(date +%s)

    DATE=$(date +"%Y-%m-%d %H:%M:%S")

    {

        echo "========================================"

        echo "Stage:       Taxonomic Confirmation"

        echo "Sample:      !{sample}"

        echo "Started:     $DATE"

        echo "Input:       !{genome}"

        echo "Database:    !{blast_db_folder}"

        echo "Software:    $(blastn -version 2>&1 | head -1)"

        echo "Parameters:  identity >=95%, max_hsps=1"

        echo "----------------------------------------"

    } > "$LOG"

    DB_INDEX="!{blast_db_folder}/strep_db"

    if [ ! -f "${DB_INDEX}.nhr" ] && \
       [ ! -f "${DB_INDEX}.00.nhr" ]; then

        echo "ERROR: BLAST database index not found:" >> "$LOG"

        echo "${DB_INDEX}" >> "$LOG"

        exit 1

    fi

    blastn \
        -query "!{genome}" \
        -db "$DB_INDEX" \
        -perc_identity 95 \
        -max_hsps 1 \
        -num_alignments 5 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart ssend evalue bitscore stitle" \
        -out "!{sample}_blast_raw.tsv" \
        >> "$LOG" 2>&1

    if [ -s "!{sample}_blast_raw.tsv" ]; then

        SPECIES_NAME=$(
            sort -k12,12nr "!{sample}_blast_raw.tsv" |
            awk -F'\t' 'NR==1 {
                if ($13 != "")
                    print $13;
                else
                    print $2;
            }'
        )

        echo "$SPECIES_NAME" > "!{sample}_taxonomy.txt"

        echo "Taxonomy: $SPECIES_NAME" >> "$LOG"

    else

        echo "No significant Streptococcus BLAST hit" > "!{sample}_taxonomy.txt"

        echo "WARNING: No BLAST hits above threshold" >> "$LOG"

    fi

    ELAPSED=$(( $(date +%s) - START ))

    {

        echo "----------------------------------------"

        echo "Completed: $(date +"%Y-%m-%d %H:%M:%S")"

        echo "Elapsed:   ${ELAPSED}s"

        echo "Status:    SUCCESS"

        echo "========================================"

    } >> "$LOG"

    '''
}
