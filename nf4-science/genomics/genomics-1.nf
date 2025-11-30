#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input
// ${projectDir} is a built-in Nextflow variable that points to the directory where the current Nextflow workflow script (genomics-1.nf) is located.
params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
params.outdir = "results_genomics"

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
        path "${input_bam}.bai"

    script:
    """
    samtools index ${input_bam}
    """

}


process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        path input_bam
        path input_bam_index // does not appear in script call as GATK knows where to look, but still need to be provided here for GATK to work with Nextflow.
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
    reads_ch = Channel.fromPath(params.reads_bam)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

}
