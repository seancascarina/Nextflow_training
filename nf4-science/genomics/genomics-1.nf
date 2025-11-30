#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input
// ${projectDir} is a built-in Nextflow variable that points to the directory where the current Nextflow workflow script (genomics-1.nf) is located.

// Run for a single file========
// params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
// Run for an array of files=========
/* params.reads_bam = ["${projectDir}/data/bam/reads_mother.bam",
                    "${projectDir}/data/bam/reads_father.bam",
                    "${projectDir}/data/bam/reads_son.bam"] */
// Run for a text file of file paths===========
params.reads_bam = "${projectDir}/data/sample_bams.txt"

params.outdir = "results_genomics"

// Accessory files
params.reference        = "${projectDir}/data/ref/ref.fasta"
params.reference_index  = "${projectDir}/data/ref/ref.fasta.fai"
params.reference_dict   = "${projectDir}/data/ref/ref.dict"
params.intervals        = "${projectDir}/data/ref/intervals.bed"

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    // Don't use mode: 'symLink' in your final workflow, this is just for demo purposes.
    // Using 'symLink' avoids copying large files, but can lead to issues if the original files are moved or deleted.
    publishDir params.outdir, mode: 'symlink'

    input:
        path input_bam

    output:
        // Output tuples so that the order is preserved for inputting into GATK_HAPLOTYPECALLER
        tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index ${input_bam}
    """

}


process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        // I'm assuming this unpacks the tuple that was passed to this process
        tuple path(input_bam), path(input_bam_index)
        // path input_bam
        // path input_bam_index // does not appear in script call as GATK knows where to look, but still need to be provided here for GATK to work with Nextflow.
        path ref_fasta
        path ref_index // does not appear in script call as GATK knows where to look, but still need to be provided here for GATK to work with Nextflow.
        path ref_dict // does not appear in script call as GATK knows where to look, but still need to be provided here for GATK to work with Nextflow.
        path interval_list

    output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx

    script:
    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.vcf \
        -L ${interval_list}
    """
}

workflow {

    // Create input channel
    // reads_ch = Channel.fromPath(params.reads_bam)
    reads_ch = Channel.fromPath(params.reads_bam).splitText()   // Splits text from a .txt file containing filepaths of bam files.

    // Load the file paths for the accessory files (reference and intervals)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

    // .view() statements to print values to terminal for manual inspection
    // reads_ch.view()
    // SAMTOOLS_INDEX.out.view()

    // Call variants from the indexed BAM files
    GATK_HAPLOTYPECALLER(
        // reads_ch,
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )
}
