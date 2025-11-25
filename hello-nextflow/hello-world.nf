#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        // "val" here indicates that this is a single value input (variable)
        val greeting

    output:
        path 'output.txt'

    script:
    // Use our val variable in the echo statement
    """
    echo '$greeting' > output.txt
    """
}

// Set default value for greeting parameter (in case user does not specify one on command line)


workflow {

    // emit a greeting
    sayHello(params.greeting)   // uses the greeting parameter
}
