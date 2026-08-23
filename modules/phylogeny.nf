/*
 * Phylogenetic Tree Construction — maximum-likelihood tree via IQ-TREE.
 *
 * Input:
 *   Core genome alignment (FASTA/ALN)
 *
 * Output:
 *   Newick tree file
 *   Phylogeny log
 *
 * Logic:
 *   - <3 genomes: Skip phylogeny gracefully
 *   - 3 genomes: IQ-TREE maximum likelihood without bootstrap
 *   - >=4 genomes: IQ-TREE ModelFinder + 1000 ultrafast bootstrap
 */


process PHYLOGENY {

    tag "phylogeny"
    label 'high_compute'

    conda 'bioconda::iqtree=2.2.0'
    container 'kizitodevbio/strepto-pipeline:v1.2.0'

    publishDir "${params.outdir}/Phylogeny",
        mode: 'copy'

    publishDir "${params.outdir}/Logs",
        mode: 'copy',
        pattern: "phylogeny.log"

    input:
    path aligned

    output:
    path "gbs_phylogeny_tree.nwk",
        emit: tree

    path "phylogeny.log",
        emit: log


    script:
    """

    #!/usr/bin/env bash
    set -euo pipefail


    LOG="phylogeny.log"

    START=\$(date +%s)

    TS=\$(date +"%Y-%m-%d %H:%M:%S")


    {
        echo "========================================"
        echo "Stage:       Phylogenetic Tree Construction"
        echo "Started:     \$TS"
        echo "Input:       ${aligned}"
        echo "Output dir:  ${params.outdir}/Phylogeny"
        echo "Software:    IQ-TREE \$(iqtree2 --version 2>&1 | head -1)"
        echo "----------------------------------------"
    } > "\$LOG"



    ############################################################
    # Validate alignment
    ############################################################


    if [ ! -s "${aligned}" ]; then

        echo "ERROR: Alignment file is missing or empty." >> "\$LOG"

        echo "()" > gbs_phylogeny_tree.nwk

        exit 0

    fi



    if grep -q "NO_ALIGNMENT" "${aligned}" 2>/dev/null; then

        echo "WARNING: No valid alignment detected." >> "\$LOG"

        echo "()" > gbs_phylogeny_tree.nwk

        exit 0

    fi



    ############################################################
    # Count genomes
    ############################################################


    SEQ_COUNT=\$(grep -c "^>" "${aligned}" || true)

    echo "Number of genomes: \$SEQ_COUNT" >> "\$LOG"



    ############################################################
    # Select phylogeny strategy
    ############################################################


    if [ "\$SEQ_COUNT" -lt 3 ]; then


        echo "WARNING: Fewer than three genomes available." >> "\$LOG"

        echo "Phylogenetic reconstruction skipped." >> "\$LOG"

        echo "()" > gbs_phylogeny_tree.nwk



    elif [ "\$SEQ_COUNT" -lt 4 ]; then


        echo "Strategy: Maximum likelihood without bootstrap" >> "\$LOG"


        iqtree2 \\
            -s "${aligned}" \\
            -m MFP \\
            -nt AUTO \\
            --redo \\
            --prefix phylogeny_out \\
            >> "\$LOG" 2>&1



    else


        echo "Strategy: Maximum likelihood with 1000 ultrafast bootstrap replicates" >> "\$LOG"


        iqtree2 \\
            -s "${aligned}" \\
            -m MFP \\
            -bb 1000 \\
            -nt AUTO \\
            --redo \\
            --prefix phylogeny_out \\
            >> "\$LOG" 2>&1


    fi



    ############################################################
    # Check tree output
    ############################################################


    if [ -f "phylogeny_out.treefile" ]; then


        cp phylogeny_out.treefile gbs_phylogeny_tree.nwk


        NUM_TIPS=\$(grep -oE '[A-Za-z0-9_.-]+:[0-9.eE+-]+' gbs_phylogeny_tree.nwk | wc -l || true)


        echo "Tree tips:   \$NUM_TIPS" >> "\$LOG"


    else


        echo "No IQ-TREE tree generated (insufficient genomes)." >> "\$LOG"


    fi



    ############################################################
    # Finish log
    ############################################################


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

