/*
 * Capsular Serotyping (CPS) — GBS serotype prediction via GBS-SBG.
 *
 * GBS-SBG (GBS Serotyping by Genome Sequencing) is a GBS-specific
 * assembly-based serotyping method using the curated GBS-SBG reference
 * database and BLASTN.
 *
 * Recognised GBS serotypes are represented by the bundled GBS-SBG
 * reference database.
 *
 * IMPORTANT:
 * The Nextflow process interface and downstream TSV schema are unchanged.
 *
 * Input:
 *   tuple(sample, fasta_file)
 *
 * Output:
 *   tuple(sample, sample_serotype.tsv)
 *
 * TSV schema:
 *   sample_id, serotype, best_match, confidence, status
 *
 * The GBS-SBG executable and reference are bundled in the Docker image:
 *
 *   /opt/GBS-SBG/GBS-SBG.pl
 *   /opt/GBS-SBG/GBS-SBG.fasta
 *
 * GBS-SBG uses BLASTN internally and creates/reuses its BLAST database
 * from the supplied local reference.
 */

process SEROTYPE_CPS {

    tag "$sample"
    label 'low_compute'

    container 'kizitodevbio/strepto-pipeline:v1.2.0'

    publishDir "${params.outdir}/Serotyping",
        mode: 'copy',
        pattern: "*_serotype.tsv"

    publishDir "${params.outdir}/Logs",
        mode: 'copy',
        pattern: "*_serotype.log"

    input:
    tuple val(sample), path(fasta_file)

    output:
    tuple val(sample), path("${sample}_serotype.tsv"), emit: results
    path "${sample}_serotype.log", emit: log

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    LOG="${sample}_serotype.log"
    START=\$(date +%s)
    TS=\$(date +"%Y-%m-%d %H:%M:%S")

    GBS_SBG="/opt/GBS-SBG/GBS-SBG.pl"
    GBS_SBG_REF="/opt/GBS-SBG/GBS-SBG.fasta"

    {
        echo "========================================"
        echo "Stage:       Capsular Serotyping (CPS)"
        echo "Sample:      ${sample}"
        echo "Started:     \$TS"
        echo "Input:       ${fasta_file}"
        echo "Method:      GBS-SBG"
        echo "Method type: Assembly-based BLASTN"
        echo "Software:    GBS-SBG"
        echo "BLAST+:      \$(blastn -version 2>&1 | head -1)"
        echo "Reference:   \$GBS_SBG_REF"
        echo "Parameters:  -best; GBS-SBG default call/uncertainty thresholds"
        echo "----------------------------------------"
    } > "\$LOG"

    # ============================================================
    # Preserve the existing downstream output schema
    # ============================================================

    printf "sample_id\\tserotype\\tbest_match\\tconfidence\\tstatus\\n" \
        > "${sample}_serotype.tsv"


    # ============================================================
    # 1. Validate bundled GBS-SBG installation
    # ============================================================

    if [ ! -f "\$GBS_SBG" ]; then

        echo "ERROR: GBS-SBG executable not found: \$GBS_SBG" >> "\$LOG"

        printf "${sample}\\tNo confident serotype\\tNA\\tNA\\tTOOL_ERROR\\n" \
            >> "${sample}_serotype.tsv"

        echo "Status: TOOL_ERROR" >> "\$LOG"

    elif [ ! -x "\$GBS_SBG" ]; then

        echo "ERROR: GBS-SBG executable is not executable: \$GBS_SBG" \
            >> "\$LOG"

        printf "${sample}\\tNo confident serotype\\tNA\\tNA\\tTOOL_ERROR\\n" \
            >> "${sample}_serotype.tsv"

        echo "Status: TOOL_ERROR" >> "\$LOG"

    elif [ ! -s "\$GBS_SBG_REF" ]; then

        echo "ERROR: GBS-SBG reference not found: \$GBS_SBG_REF" \
            >> "\$LOG"

        printf "${sample}\\tNo confident serotype\\tNA\\tNA\\tTOOL_ERROR\\n" \
            >> "${sample}_serotype.tsv"

        echo "Status: TOOL_ERROR" >> "\$LOG"

    elif ! command -v blastn >/dev/null 2>&1; then

        echo "ERROR: BLASTN was not found in the Docker image." \
            >> "\$LOG"

        printf "${sample}\\tNo confident serotype\\tNA\\tNA\\tTOOL_ERROR\\n" \
            >> "${sample}_serotype.tsv"

        echo "Status: TOOL_ERROR" >> "\$LOG"

    else

        echo "GBS-SBG executable: \$GBS_SBG" >> "\$LOG"
        echo "GBS-SBG reference:  \$GBS_SBG_REF" >> "\$LOG"


        # ========================================================
        # 2. Run GBS-SBG
        #
        # GBS-SBG:
        #   - accepts an assembled FASTA
        #   - uses BLASTN
        #   - uses the supplied local GBS-SBG reference
        #   - makes/reuses BLAST databases internally
        #   - -best reports only the primary serotype call and
        #     associated uncertainty information
        #
        # All GBS-SBG output is written to STDOUT.
        # ========================================================

        SBG_OUTPUT=""

        set +e

        SBG_OUTPUT=\$(
            "\$GBS_SBG" \
                "${fasta_file}" \
                -name "${sample}" \
                -best \
                -blastn "\$(command -v blastn)" \
                -ref "\$GBS_SBG_REF" \
                2>> "\$LOG"
        )

        SBG_EXIT=\$?

        set -e

        echo "GBS-SBG exit status: \$SBG_EXIT" >> "\$LOG"


        # ========================================================
        # 3. Handle execution failure
        # ========================================================

        if [ "\$SBG_EXIT" -ne 0 ]; then

            echo "ERROR: GBS-SBG failed." >> "\$LOG"

            if [ -n "\$SBG_OUTPUT" ]; then
                echo "GBS-SBG STDOUT:" >> "\$LOG"
                echo "\$SBG_OUTPUT" >> "\$LOG"
            fi

            printf "${sample}\\tNo confident serotype\\tNA\\tNA\\tTOOL_ERROR\\n" \
                >> "${sample}_serotype.tsv"

            echo "Status: TOOL_ERROR" >> "\$LOG"


        elif [ -z "\$SBG_OUTPUT" ]; then

            echo "ERROR: GBS-SBG returned no output." >> "\$LOG"

            printf "${sample}\\tNo confident serotype\\tNA\\tNA\\tNO_MATCH\\n" \
                >> "${sample}_serotype.tsv"

            echo "Status: NO_MATCH" >> "\$LOG"


        else

            echo "GBS-SBG raw output:" >> "\$LOG"
            echo "\$SBG_OUTPUT" >> "\$LOG"


            # ====================================================
            # 4. Parse the actual GBS-SBG output format
            #
            # Non-verbose output:
            #
            # # Name    Serotype    Uncertainty
            # sample   III         Coverage:0.98;...
            #
            # For a clean confident call, the third field may be empty.
            #
            # For a non-typeable result:
            #
            # sample   NT          MaxCov:...;MaxID:...
            # ====================================================

            RESULT_LINE=\$(
                printf '%s\\n' "\$SBG_OUTPUT" |
                awk -F '\\t' '
                    \$1 !~ /^#/ &&
                    NF >= 2 {
                        print
                        exit
                    }
                '
            )


            if [ -z "\$RESULT_LINE" ]; then

                echo "ERROR: No interpretable GBS-SBG result line." \
                    >> "\$LOG"

                printf "${sample}\\tNo confident serotype\\tNA\\tNA\\tNO_MATCH\\n" \
                    >> "${sample}_serotype.tsv"

                echo "Status: NO_MATCH" >> "\$LOG"

            else

                # ------------------------------------------------
                # Actual fields from GBS-SBG:
                #
                # Field 1 = sample name
                # Field 2 = serotype
                # Field 3 = uncertainty information
                # ------------------------------------------------

                RESULT_SAMPLE=\$(
                    printf '%s\\n' "\$RESULT_LINE" |
                    cut -f1
                )

                RESULT_SEROTYPE=\$(
                    printf '%s\\n' "\$RESULT_LINE" |
                    cut -f2
                )

                RESULT_UNCERTAINTY=\$(
                    printf '%s\\n' "\$RESULT_LINE" |
                    cut -f3-
                )


                # ------------------------------------------------
                # Normalize confidence/uncertainty
                # ------------------------------------------------

                if [ -z "\$RESULT_UNCERTAINTY" ]; then
                    CONFIDENCE="HIGH"
                else
                    CONFIDENCE="\$RESULT_UNCERTAINTY"
                fi


                # ------------------------------------------------
                # Handle non-typeable result
                # ------------------------------------------------

                if [ "\$RESULT_SEROTYPE" = "NT" ] || \
                   [ -z "\$RESULT_SEROTYPE" ]; then

                    SEROTYPE="No confident serotype"
                    BEST_MATCH="NT"
                    STATUS="NO_MATCH"

                    if [ -z "\$RESULT_UNCERTAINTY" ]; then
                        CONFIDENCE="NA"
                    fi

                else

                    BEST_MATCH="\$RESULT_SEROTYPE"

                    SEROTYPE="\$RESULT_SEROTYPE"
                    SEROTYPE="\${SEROTYPE#GBS-SBG:}"

                    STATUS="OK"

                fi

                
                # ------------------------------------------------
                # Write the established GBS-Sentinel result schema
                # ------------------------------------------------

                printf "%s\\t%s\\t%s\\t%s\\t%s\\n" \
                    "${sample}" \
                    "\${SEROTYPE}" \
                    "\${BEST_MATCH}" \
                    "\${CONFIDENCE}" \
                    "\${STATUS}" \
                    >> "${sample}_serotype.tsv"


                echo "Serotype:           \$SEROTYPE" >> "\$LOG"
                echo "Best match:         \$BEST_MATCH" >> "\$LOG"
                echo "Confidence/status:  \$CONFIDENCE" >> "\$LOG"
                echo "Status:             \$STATUS" >> "\$LOG"

            fi

        fi

    fi


    # ============================================================
    # Completion
    # ============================================================

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
