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


process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
        // "val" here indicates that this is a single value input (variable)
        path input_file

    output:
        path "UPPER-${input_file}.txt"     // double quotes are necessary here for variable interpolation. Single quotes are treated as a string literal.

    script:
    // Use our val variable in the echo statement
    """
    cat $input_file | tr '[a-z]' '[A-Z]' > UPPER-${input_file}.txt
    """
}


// Set default value for greeting parameter (in case 
// user does not specify one on command line)
params.greeting = 'greetings.csv'

workflow {

    // explicitly create a channel with an array of different greetings
    // greetings = ['Hello', 'Bonjour', 'Hola', 'Ciao', 'Hallo']
    greeting_ch = channel.fromPath(params.greeting)
                        .view{csv -> "Before splitCsv: $csv"}
                        .splitCsv()
                        .view{csv -> "After splitCSV: $csv"}
                        .map{item -> item[0]}  // get first column from CSV.
                        .view{csv -> "After map: $csv"}

    // emit a greeting
    sayHello(greeting_ch)   // uses the greeting channel

    convertToUpper(sayHello.out)  // use the output files from sayHello as input to convertToUpper
}
