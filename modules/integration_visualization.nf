/*
 * GBS-Sentinel Integration & Visualization — final workflow stage.
 *
 * Consumes real outputs from:
 *
 *   - AMRFinderPlus
 *   - Virulence-factor analysis
 *   - MLST
 *   - GBS-SBG capsular serotyping
 *   - Core-genome phylogeny
 *
 * Generates:
 *
 *   - Integrated surveillance summary table
 *   - AMR visualizations
 *   - Virulence visualizations
 *   - MLST visualizations
 *   - Capsular serotype visualizations
 *   - ST × serotype relationships
 *   - ST × AMR relationships
 *   - ST × virulence relationships
 *   - Serotype × AMR relationships
 *   - Serotype × virulence relationships
 *   - Integrated AMR/virulence burden plots
 *   - Annotated phylogenies
 *
 * All figures are derived from actual pipeline outputs.
 * No synthetic biological values are generated.
 */

process INTEGRATION_VISUALIZATION {

    tag "visualization"

    label 'low_compute'

    container 'kizitodevbio/strepto-pipeline:v1.2.0'


    // ------------------------------------------------------------------
    // Publish visualization outputs
    // ------------------------------------------------------------------

    publishDir "${params.outdir}/Figures",
        mode: 'copy',
        pattern: "*.png"

    publishDir "${params.outdir}/Figures",
        mode: 'copy',
        pattern: "*.svg"

    publishDir "${params.outdir}/Figures",
        mode: 'copy',
        pattern: "*.pdf"

    publishDir "${params.outdir}/Figures",
        mode: 'copy',
        pattern: "*.tsv"

    publishDir "${params.outdir}/Reports",
        mode: 'copy',
        pattern: "Reports/**"

    publishDir "${params.outdir}/Reports",
        mode: 'copy',
        pattern: "surveillance_summary.tsv"


    // ------------------------------------------------------------------
    // Inputs
    // ------------------------------------------------------------------

    input:

    path tree

    path amr_files

    path vf_files

    path mlst_files

    path serotype_files


    // ------------------------------------------------------------------
    // Outputs
    // ------------------------------------------------------------------

    output:

    path "*.png",
        emit: figures_png

    path "*.svg",
        emit: figures_svg

    path "*.pdf",
        emit: figures_pdf

    path "*.tsv",
        emit: tables

    path "surveillance_summary.tsv",
        emit: surveillance_summary

    path "Reports/**",
        emit: reports

    path "integration_viz.log",
        emit: log


    // ------------------------------------------------------------------
    // Execution
    // ------------------------------------------------------------------

    script:
    """

    #!/usr/bin/env bash

    set -euo pipefail


    LOG="integration_viz.log"
    START=\$(date +%s)


    # ------------------------------------------------------------------
    # Headless plotting configuration
    # ------------------------------------------------------------------

    export QT_QPA_PLATFORM=offscreen
    export MPLBACKEND=Agg
    export DISPLAY=:99


    {
        echo "============================================================"
        echo "GBS-Sentinel Integration & Visualization"
        echo "Started:       \$(date)"
        echo "Container:     kizitodevbio/strepto-pipeline:v1.2.0"
        echo "============================================================"
    } > "\$LOG"


    # ------------------------------------------------------------------
    # Prepare input directories
    # ------------------------------------------------------------------

    mkdir -p amr_inputs
    mkdir -p vf_inputs
    mkdir -p mlst_inputs
    mkdir -p serotype_inputs
    mkdir -p Reports


    # ------------------------------------------------------------------
    # Collect AMRFinderPlus results
    # ------------------------------------------------------------------

    echo "Collecting AMRFinderPlus inputs..." >> "\$LOG"

    for f in ${amr_files}; do

        cp "\$f" amr_inputs/

        echo "AMRFinderPlus: \$(basename "\$f")" >> "\$LOG"

    done


    # ------------------------------------------------------------------
    # Collect virulence results
    # ------------------------------------------------------------------

    echo "Collecting virulence inputs..." >> "\$LOG"

    for f in ${vf_files}; do

        cp "\$f" vf_inputs/

        echo "Virulence: \$(basename "\$f")" >> "\$LOG"

    done


    # ------------------------------------------------------------------
    # Collect MLST results
    # ------------------------------------------------------------------

    echo "Collecting MLST inputs..." >> "\$LOG"

    for f in ${mlst_files}; do

        cp "\$f" mlst_inputs/

        echo "MLST: \$(basename "\$f")" >> "\$LOG"

    done


    # ------------------------------------------------------------------
    # Collect GBS-SBG serotyping results
    # ------------------------------------------------------------------

    echo "Collecting GBS-SBG serotyping inputs..." >> "\$LOG"

    for f in ${serotype_files}; do

        cp "\$f" serotype_inputs/

        echo "Serotype: \$(basename "\$f")" >> "\$LOG"

    done


    # ------------------------------------------------------------------
    # Collect phylogenetic tree
    # ------------------------------------------------------------------

    cp ${tree} tree_input.nwk

    echo "Phylogenetic tree staged as: tree_input.nwk" >> "\$LOG"

    # ------------------------------------------------------------------
    # Verify the visualization environment
    # ------------------------------------------------------------------

    echo "Checking Python visualization environment..." >> "\$LOG"

    python3 - <<'PY' >> "\$LOG" 2>&1

import matplotlib

matplotlib.use("Agg")

import pandas
import seaborn

from ete3 import Tree

tree = Tree("tree_input.nwk")

print("Python visualization stack: OK")
print("pandas:", pandas.__version__)
print("seaborn:", seaborn.__version__)
print("ETE3: OK")
print("Tree leaves:", len(tree))

PY


    # ------------------------------------------------------------------
    # Verify the supplied biological result files
    # ------------------------------------------------------------------

    echo "Checking supplied result files..." >> "\$LOG"

    echo "AMRFinderPlus files:" >> "\$LOG"
    find amr_inputs -maxdepth 1 -type f -print >> "\$LOG"

    echo "Virulence files:" >> "\$LOG"
    find vf_inputs -maxdepth 1 -type f -print >> "\$LOG"

    echo "MLST files:" >> "\$LOG"
    find mlst_inputs -maxdepth 1 -type f -print >> "\$LOG"

    echo "Serotyping files:" >> "\$LOG"
    find serotype_inputs -maxdepth 1 -type f -print >> "\$LOG"


    # ------------------------------------------------------------------
    # Run integrated figure generation
    # ------------------------------------------------------------------

    echo "Running GBS-Sentinel integrated figure generation..." \
        >> "\$LOG"

    python3 ${projectDir}/bin/generate_figures.py \
        --amr-dir amr_inputs \
        --vf-dir vf_inputs \
        --mlst-dir mlst_inputs \
        --serotype-dir serotype_inputs \
        --tree tree_input.nwk \
        --outdir . \
        >> "\$LOG" 2>&1


    # ------------------------------------------------------------------
    # Validate canonical surveillance output
    # ------------------------------------------------------------------

    if [ ! -s surveillance_summary.tsv ]; then

        echo "ERROR: surveillance_summary.tsv was not generated." \
            >> "\$LOG"

        exit 1

    fi


    # ------------------------------------------------------------------
    # Copy canonical surveillance summary into Reports
    # ------------------------------------------------------------------

    cp surveillance_summary.tsv \
        Reports/surveillance_summary.tsv


    # ------------------------------------------------------------------
    # Record final integrated component summary
    # ------------------------------------------------------------------

    cat > Reports/visualization_summary.txt <<EOF

GBS-Sentinel
Integrated Surveillance Visualization Summary
==============================================

Status: SUCCESS

Biological data integrated:
- AMRFinderPlus resistance determinants
- AMR class/subclass information
- Virulence-factor detections
- MLST sequence types
- GBS-SBG capsular serotypes
- Core-genome phylogeny

Primary integrated table:
- surveillance_summary.tsv

Visualization classes:
- AMR determinant frequency
- AMR determinant heatmap
- AMR class distribution
- AMR detection-method distribution
- AMR co-occurrence
- Virulence determinant frequency
- Virulence determinant heatmap
- Virulence co-occurrence
- MLST sequence-type distribution
- GBS capsular serotype distribution
- MLST × serotype
- MLST × AMR burden
- MLST × virulence burden
- Serotype × AMR burden
- Serotype × virulence burden
- Integrated AMR/virulence burden
- AMR burden versus virulence burden
- Surveillance dashboard
- Phylogeny with ST labels
- Phylogeny with serotype labels
- Phylogeny with AMR burden
- Phylogeny with virulence burden
- Integrated surveillance phylogeny

All figures are generated from actual pipeline outputs.
No synthetic biological observations are generated.

EOF


    # ------------------------------------------------------------------
    # Completion logging
    # ------------------------------------------------------------------

    ELAPSED=\$(( \$(date +%s) - START ))

    {
        echo "------------------------------------------------------------"
        echo "Completed: \$(date)"
        echo "Elapsed:   \${ELAPSED}s"
        echo "Status:    SUCCESS"
        echo "------------------------------------------------------------"
    } >> "\$LOG"

    """
}
