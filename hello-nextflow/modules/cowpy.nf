#! /usr/bin/env nextflow

process cowpy {

    publishDir 'results', mode: 'copy'
    container "community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273"
    // conda "conda-forge::cowpy=1.1.5"  //Not used here, but demonstrates how to specify a conda package instead of a container.

    input:
        path input_file
        val character

    output:
        path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c $character > cowpy-${input_file}
    """
}