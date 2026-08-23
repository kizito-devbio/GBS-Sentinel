process ALIGNMENT {
    label 'med_compute'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'step_log.txt'

    input:
    path combined_core

    output:
    // This 'emit' fixes the path value error for the next step
    path "aligned.fasta", emit: alignment
    path "step_log.txt",  emit: log

    script:
    """
    timestamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\$timestamp] Step: Multiple Sequence Alignment using MAFFT." > step_log.txt
    
    mafft --auto $combined_core > aligned.fasta || { echo "MAFFT failed"; exit 1; }
    
    num_aligned=\$(grep -c '^>' aligned.fasta)
    echo "[\$timestamp] Number of aligned sequences: \$num_aligned" >> step_log.txt
    """
}
