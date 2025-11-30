#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir params.outdir, mode: 'symlink'

    input:
        path input_bam

    output:
        path "${input_bam}.bai"

    script:
    """
    samtools index ${input_bam}
    """

}

workflow {

    // Create input channel

    // Create index file for input BAM file

}
