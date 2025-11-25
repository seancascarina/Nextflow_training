#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    output:
        path 'output.txt'

    script:
    """
    echo 'Hello World!' > different_output_name.txt
    """
}

workflow {

    // emit a greeting
    sayHello()
}
