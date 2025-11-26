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

// Set default value for greeting parameter (in case 
// user does not specify one on command line)
params.greeting = 'Hello, world!'

workflow {

    // explicitly create a channel with a greeting message
    greeting_ch = channel.of('Hello Channels!')

    // emit a greeting
    sayHello(greeting_ch)   // uses the greeting channel
}
