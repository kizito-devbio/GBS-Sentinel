#!/usr/bin/env nextflow

/*
 * GBS-Sentinel
 *
 * Modular Nextflow workflow for genomic surveillance of
 * Streptococcus agalactiae (Group B Streptococcus, GBS).
 *
 * Core surveillance outputs:
 *   - Genome quality assessment
 *   - MLST sequence typing
 *   - GBS-specific capsular serotyping (GBS-SBG)
 *   - Antimicrobial-resistance detection (AMRFinderPlus)
 *   - Virulence-factor detection
 *   - Functional annotation
 *   - Core-genome analysis
 *   - Phylogenetic inference
 *   - Integrated surveillance visualization
 *
 * Foundation:
 *   GBS-Genomics-Pipeline
 *
 * GBS-Sentinel extensions:
 *   GBS-SBG capsular serotyping
 *   AMRFinderPlus resistome detection
 *
 * Input pathways:
 *   --curated_dir  Pre-assembled FASTA files
 *   --raw_dir      Paired-end FASTQ files
 *
 * Supported genome formats:
 *   *.fa
 *   *.fna
 *   *.fasta
 *
 * Supported paired-end read formats:
 *   *_1 / *_2
 *   *_R1 / *_R2
 *   *.fastq
 *   *.fastq.gz
 *   *.fq
 *   *.fq.gz
 *
 * Human-genome decontamination is not part of the active workflow.
 * Previously disabled/frozen modules remain under modules/frozen/.
 */

nextflow.enable.dsl = 2


// ============================================================================
// MODULE IMPORTS
// ============================================================================

// Existing pipeline modules
include { FASTP_QC }                         from './modules/qc.nf'
include { ASSEMBLY }                   from './modules/assembly.nf'
include { QUALITY_ASSESS }             from './modules/quality_assess.nf'
include { BACKGROUND_SELECTION }        from './modules/background_selection.nf'
include { FUNCTIONAL_ANNOTATION }      from './modules/functional_annotation.nf'
include { VIRULENCE_FACTOR }           from './modules/virulence_factor.nf'
include { MLST_EXTRACTION }             from './modules/mlst_extraction.nf'
include { CORE_GENOME }                from './modules/core_genome.nf'
include { PHYLOGENY }                  from './modules/phylogeny.nf'
include { INTEGRATION_VISUALIZATION }  from './modules/integration_visualization.nf'

// GBS-Sentinel surveillance extensions
include { SEROTYPE_CPS }               from './modules/serotyping/serotype_cps.nf'
include { AMRFINDERPLUS }              from './modules/amrfinderplus/amrfinderplus.nf'


// ============================================================================
// HELP
// ============================================================================

def printHelp() {

    log.info """

    GBS-Sentinel
    ======================================================================

    Modular genomic surveillance workflow for
    Streptococcus agalactiae (Group B Streptococcus, GBS)

    Surveillance components:
      * Genome quality assessment
      * MLST sequence typing
      * GBS-specific capsular serotyping (GBS-SBG)
      * Antimicrobial-resistance detection (AMRFinderPlus)
      * Virulence-factor profiling
      * Functional annotation
      * Core-genome analysis
      * Phylogenetic analysis
      * Integrated surveillance visualization

    INPUT
    ----------------------------------------------------------------------

    Required (provide exactly one):

      --raw_dir <directory>
          Paired-end FASTQ/FASTQ.GZ files

      --curated_dir <directory>
          Pre-assembled FASTA genomes

    OUTPUT
    ----------------------------------------------------------------------

      --outdir <directory>
          Output directory
          Default: ./results

    OPTIONS
    ----------------------------------------------------------------------

      --mlst_scheme <string>
          MLST scheme
          Default: sagalactiae

      --min_n50 <integer>
          Minimum assembly N50 for raw-read assemblies
          Default: 10000

      --max_cpus <integer>
          Maximum CPU cores

      --max_memory <string>
          Maximum memory
          Default: 6 GB

    PROFILES
    ----------------------------------------------------------------------

      -profile docker
      -profile singularity
      -profile conda
      -profile cluster
      -profile test

    EXAMPLES
    ----------------------------------------------------------------------

    Curated genomes:

      nextflow run pipeline.nf \
          -profile docker \
          --curated_dir curated_data \
          --outdir results

    Raw reads:

      nextflow run pipeline.nf \
          -profile docker \
          --raw_dir raw_reads \
          --outdir results

    Resume:

      nextflow run pipeline.nf \
          -profile docker \
          --curated_dir curated_data \
          --outdir results \
          -resume

    ======================================================================
    """

}


// ============================================================================
// PARAMETER VALIDATION
// ============================================================================

def validateParams() {

    if (params.help) {
        printHelp()
        exit 0
    }

    if (!params.curated_dir && !params.raw_dir) {

        error """
        No input supplied.

        Provide exactly one of:

          --curated_dir <directory>

        or

          --raw_dir <directory>
        """
    }

    if (params.curated_dir && params.raw_dir) {

        error """
        Both --curated_dir and --raw_dir were supplied.

        Provide only one input pathway.
        """
    }

    if (params.curated_dir &&
        !file(params.curated_dir).exists()) {

        error """
        Curated genome directory does not exist:

          ${params.curated_dir}
        """
    }

    if (params.raw_dir &&
        !file(params.raw_dir).exists()) {

        error """
        Raw-read directory does not exist:

          ${params.raw_dir}
        """
    }
}


// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {

    // ------------------------------------------------------------------------
    // Parameter validation
    // ------------------------------------------------------------------------

    validateParams()


    // ------------------------------------------------------------------------
    // Pipeline banner
    // ------------------------------------------------------------------------

    log.info """

    ╔══════════════════════════════════════════════════════════════════════╗
    ║                            GBS-Sentinel                              ║
    ║                                                                      ║
    ║  Genomic Surveillance · Capsular Serotyping · Resistome Monitoring  ║
    ║  Streptococcus agalactiae (Group B Streptococcus)                    ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Input:
      ${params.curated_dir
          ? "Curated genomes (${params.curated_dir})"
          : "Raw reads (${params.raw_dir})"}

    Output:
      ${params.outdir}

    Profile:
      ${workflow.profile ?: 'default'}

    Container:
      kizitodevbio/strepto-pipeline:v1.2.0

    Serotyping:
      GBS-SBG

    AMR:
      AMRFinderPlus 4.2.7

    """


    // ========================================================================
    // 1. INPUT ACQUISITION
    // ========================================================================

    Channel ch_genomes


    // ------------------------------------------------------------------------
    // Curated genome pathway
    // ------------------------------------------------------------------------

    if (params.curated_dir) {

        ch_genomes = Channel
            .fromPath(
                "${params.curated_dir}/*.{fa,fna,fasta}",
                checkIfExists: true
            )
            .ifEmpty {
                error """
                No curated FASTA genomes were found in:

                  ${params.curated_dir}

                Expected extensions:

                  *.fa
                  *.fna
                  *.fasta
                """
            }
            .map { fasta ->

                tuple(
                    fasta.baseName,
                    fasta
                )
            }

    } else {

        // --------------------------------------------------------------------
        // Raw paired-end FASTQ pathway
        // --------------------------------------------------------------------

        def raw_channels = [

            Channel.fromFilePairs(
                "${params.raw_dir}/*_{1,2}.fastq",
                flat: true
            ),

            Channel.fromFilePairs(
                "${params.raw_dir}/*_{1,2}.fastq.gz",
                flat: true
            ),

            Channel.fromFilePairs(
                "${params.raw_dir}/*_R{1,2}.fastq",
                flat: true
            ),

            Channel.fromFilePairs(
                "${params.raw_dir}/*_R{1,2}.fastq.gz",
                flat: true
            ),

            Channel.fromFilePairs(
                "${params.raw_dir}/*_{1,2}.fq",
                flat: true
            ),

            Channel.fromFilePairs(
                "${params.raw_dir}/*_{1,2}.fq.gz",
                flat: true
            ),

            Channel.fromFilePairs(
                "${params.raw_dir}/*_R{1,2}.fq",
                flat: true
            ),

            Channel.fromFilePairs(
                "${params.raw_dir}/*_R{1,2}.fq.gz",
                flat: true
            )
        ]


        ch_raw_reads = raw_channels
            .inject { acc, next_channel ->
                acc.mix(next_channel)
            }
            .ifEmpty {

                error """
                No paired-end FASTQ files were found in:

                  ${params.raw_dir}

                Supported naming patterns include:

                  sample_1.fastq
                  sample_2.fastq

                  sample_1.fastq.gz
                  sample_2.fastq.gz

                  sample_R1.fastq
                  sample_R2.fastq

                  sample_R1.fastq.gz
                  sample_R2.fastq.gz

                  sample_1.fq
                  sample_2.fq

                  sample_1.fq.gz
                  sample_2.fq.gz

                  sample_R1.fq
                  sample_R2.fq

                  sample_R1.fq.gz
                  sample_R2.fq.gz
                """
            }


        // --------------------------------------------------------------------
        // Raw reads → QC
        // --------------------------------------------------------------------

        qc_out = FASTP_QC(
            ch_raw_reads
        )


        // --------------------------------------------------------------------
        // QC → assembly
        // --------------------------------------------------------------------

        assembly_out = ASSEMBLY(
            qc_out.trimmed
        )


        // --------------------------------------------------------------------
        // Assembly → quality assessment
        // --------------------------------------------------------------------

        qa_out = QUALITY_ASSESS(
            assembly_out.assembled
        )


        // --------------------------------------------------------------------
        // Assembly quality filtering
        // --------------------------------------------------------------------

        ch_genomes = qa_out.metrics

            .filter {
                meta,
                fasta,
                n50_file,
                total_len,
                contigs ->

                n50_file
                    .text
                    .trim()
                    .toInteger() >= params.min_n50
            }

            .map {
                meta,
                fasta,
                n50,
                total_len,
                contigs ->

                tuple(
                    meta,
                    fasta
                )
            }
    }


    // ========================================================================
    // 2. BACKGROUND GENOME SELECTION
    // ========================================================================

    /*
     * The previous BLAST taxonomy database stage has been removed from the
     * active GBS-Sentinel execution path.
     *
     * This does not affect:
     *
     *   - MLST
     *   - GBS-SBG serotyping
     *   - AMRFinderPlus
     *   - Virulence analysis
     *   - Functional annotation
     *   - Core-genome analysis
     *   - Phylogeny
     *
     * Background selection remains available for phylogenetic context.
     */

    BACKGROUND_SELECTION(
        ch_genomes
            .map { sample, fasta -> fasta }
            .collect()
    )


    // ========================================================================
    // 3. FUNCTIONAL ANNOTATION
    // ========================================================================

    /*
     * Genome
     *   ↓
     * Prokka / existing functional annotation
     *   ↓
     * Annotated genome directory
     */

    ch_annot_results = FUNCTIONAL_ANNOTATION(
        ch_genomes
    ).annotation_results


    // ========================================================================
    // 4. PREPARE ANNOTATED GENOMES FOR LEGACY DOWNSTREAM TOOLS
    // ========================================================================

    /*
     * Existing modules MLST and Virulence Factor analysis consume the
     * Prokka-derived FASTA. Their existing interfaces are preserved.
     */

    ch_fna_for_tools = ch_annot_results

        .map {
            sample,
            amr,
            prokka_dir ->

            tuple(
                sample,
                file("${prokka_dir}/${sample}.fna")
            )
        }


    // ========================================================================
    // 5. MLST
    // ========================================================================

    /*
     * Genome
     *   ↓
     * MLST
     *   ↓
     * Sample-specific sequence type
     */

    ch_mlst_out = MLST_EXTRACTION(
        ch_fna_for_tools
    ).results


    // ========================================================================
    // 6. GBS-SBG CAPSULAR SEROTYPING
    // ========================================================================

    /*
     * GBS-specific capsular serotyping.
     *
     * Genome
     *   ↓
     * GBS-SBG
     *   ↓
     * Serotype / uncertainty information
     *
     * The GBS-SBG executable and reference are bundled in Docker v1.2.0.
     *
     * The existing SEROTYPE_CPS process interface is preserved so the
     * downstream output remains:
     *
     *   sample_id
     *   serotype
     *   best_match
     *   confidence
     *   status
     */

    ch_serotype_out = SEROTYPE_CPS(
        ch_genomes
    ).results


    // ========================================================================
    // 7. AMRFINDERPLUS
    // ========================================================================

    /*
     * Primary AMR surveillance result.
     *
     * Genome
     *   ↓
     * AMRFinderPlus 4.2.7
     *   ↓
     * Resistance determinants
     *   ↓
     * Gene / class / subclass / method information
     */

    ch_amrfinder_out = AMRFINDERPLUS(
        ch_genomes
    ).results


    // ========================================================================
    // 8. VIRULENCE FACTORS
    // ========================================================================

    /*
     * Existing virulence-factor workflow remains unchanged.
     */

    ch_vf_out = VIRULENCE_FACTOR(
        ch_fna_for_tools
    ).results


    // ========================================================================
    // 9. CORE GENOME ANALYSIS
    // ========================================================================

    /*
     * Prokka GFFs
     *   ↓
     * Panaroo/core-genome processing
     *   ↓
     * Core-genome alignment
     */

    ch_gffs = ch_annot_results

        .map {
            sample,
            amr,
            prokka_dir ->

            file(
                "${prokka_dir}/${sample}.gff"
            )
        }

        .collect()


    ch_core = CORE_GENOME(
        ch_gffs
    )


    // ========================================================================
    // 10. PHYLOGENETIC INFERENCE
    // ========================================================================

    /*
     * Core-genome alignment
     *   ↓
     * Phylogenetic inference
     *   ↓
     * Newick tree
     */

    ch_tree = PHYLOGENY(
        ch_core.alignment
    )


    // ========================================================================
    // 11. INTEGRATED SURVEILLANCE VISUALIZATION
    // ========================================================================

    /*
     * The visualization stage now receives ALL surveillance dimensions:
     *
     *   1. Core-genome phylogeny
     *   2. AMRFinderPlus
     *   3. Virulence
     *   4. MLST
     *   5. GBS-SBG serotyping
     *
     * Importantly, AMRFinderPlus is passed directly rather than using
     * the older AMR output embedded in FUNCTIONAL_ANNOTATION.
     *
     * This ensures the final surveillance figures represent the actual
     * AMRFinderPlus results generated by the GBS-Sentinel workflow.
     */

    INTEGRATION_VISUALIZATION(

        // ---------------------------------------------------------------
        // Core-genome phylogeny
        // ---------------------------------------------------------------

        ch_tree
            .tree
            .collect(),

        // ---------------------------------------------------------------
        // AMRFinderPlus
        // ---------------------------------------------------------------

        ch_amrfinder_out
            .map {
                sample,
                result ->

                result
            }
            .collect(),

        // ---------------------------------------------------------------
        // Virulence
        // ---------------------------------------------------------------

        ch_vf_out
            .map {
                sample,
                result ->

                result
            }
            .collect(),

        // ---------------------------------------------------------------
        // MLST
        // ---------------------------------------------------------------

        ch_mlst_out
            .map {
                sample,
                result ->

                result
            }
            .collect(),

        // ---------------------------------------------------------------
        // GBS-SBG serotyping
        // ---------------------------------------------------------------

        ch_serotype_out
            .map {
                sample,
                result ->

                result
            }
            .collect()
    )
}

