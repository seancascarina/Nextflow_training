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
        path "$greeting-output.txt"     // double quotes are necessary here for variable interpolation. Single quotes are treated as a string literal.

    script:
    // Use our val variable in the echo statement
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

// Set default value for greeting parameter (in case 
// user does not specify one on command line)
params.greeting = 'Hello, world!'

workflow {

    // explicitly create a channel with a greeting message
    greeting_ch = channel.of('Hello', 'Bonjour', 'Hola', 'Ciao', 'Hallo')

    // emit a greeting
    sayHello(greeting_ch)   // uses the greeting channel
}
