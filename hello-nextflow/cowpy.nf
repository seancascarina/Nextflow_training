#! /usr/bin/env nextflow

process cowpy {

    publishDir 'results', mode: 'copy'

    input:
        path input_file
        val character

    output:
        cowpy-${input_file}

    script:
    """
    cat $input_file | cowpy -c $character > cowpy-${input_file}
    """
}