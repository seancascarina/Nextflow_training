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

    // explicitly create a channel with an array of different greetings
    greetings = ['Hello', 'Bonjour', 'Hola', 'Ciao', 'Hallo']
    greeting_ch = channel.of(greetings)
                        // This view operator with curly brackets is like list comprehension...it runs with a temporary local variable name and does it for every item in the items. Before the flatten, there is only one item, so the greetings array is printed.
                         .view {greeting -> "Before flatten: $greeting"}  // view operator to print the channel contents. view is like Python print()
                         .flatten()  // flatten the channel to emit each greeting separately
                         // In this view operator, each greeting is a separate item (from .flatten()), so this runs for each greeting.
                         .view {greeting -> 'After flatten: ' + greeting}   // view operator to print the channel contents
                        // Can use either "$greeting" variable interpolation or "+ greeting" string concatenation.

    // emit a greeting
    sayHello(greeting_ch)   // uses the greeting channel
}
