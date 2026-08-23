process DECONTAM {
    tag "${sample}"
    label 'med_compute'

    // Publish clean reads to the sample folder and logs to the log folder
    publishDir "${params.outdir}/decontamination/${sample}", mode: 'copy', pattern: "${sample}_clean_{1,2}.fastq.gz"
    publishDir "${params.outdir}/logs/decontamination", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample), path(reads)
    path human_ref
    path human_indices

    output:
    tuple val(sample), path("${sample}_clean_{1,2}.fastq.gz"), emit: cleaned
    path "${sample}_decontamination.log", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOGFILE="${sample}_decontamination.log"
    DATE_CMD='date "+%Y-%m-%d %H:%M:%S"'

    log() {
        echo "[ \$( \$DATE_CMD ) ] \$1" | tee -a "\${LOGFILE}"
    }

    log "START decontamination for sample: ${sample}"
    log "Using reference: ${human_ref}"
    log "Input reads: ${reads[0]} and ${reads[1]}"

    # Run BWA MEM: 
    # -f 12 keeps only reads where both pairs are unmapped to the human ref
    # -F 256 excludes secondary alignments
    bwa mem -t ${task.cpus} ${human_ref} ${reads[0]} ${reads[1]} 2>> "\${LOGFILE}" | \
    samtools view --threads ${task.cpus} -b -f 12 -F 256 - > unmapped.bam

    if [ ! -s unmapped.bam ]; then
        log "ERROR: unmapped.bam is empty or missing. Check human reference/indices."
        exit 1
    fi

    log "Extracting clean paired-end reads from BAM"
    samtools fastq --threads ${task.cpus} \
        -1 "${sample}_clean_1.fastq.gz" \
        -2 "${sample}_clean_2.fastq.gz" \
        -0 /dev/null -s /dev/null -n \
        unmapped.bam 2>> "\${LOGFILE}"

    # Verify output files exist and are not empty
    if [ ! -s "${sample}_clean_1.fastq.gz" ] || [ ! -s "${sample}_clean_2.fastq.gz" ]; then
        log "ERROR: Output clean FASTQs are empty. Decontamination might have removed all reads."
        exit 1
    fi

    # Count reads (Total lines / 4)
    n_r1=\$(zcat "${sample}_clean_1.fastq.gz" | wc -l | awk '{print \$1/4}')
    n_r2=\$(zcat "${sample}_clean_2.fastq.gz" | wc -l | awk '{print \$1/4}')

    log "Read pairs remaining after decontamination: \${n_r1}"
    
    rm -f unmapped.bam
    log "FINISHED decontamination for sample: ${sample}"
    """
}
