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


process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        // "val" here indicates that this is a single value input (variable)
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt"     // double quotes are necessary here for variable interpolation. Single quotes are treated as a string literal.

    script:
    // Use our val variable in the echo statement
    """
    cat $input_files > "COLLECTED-${batch_name}-output.txt"
    """
}


// Set default value for greeting parameter (in case 
// user does not specify one on command line)
params.greeting = 'greetings.csv'
params.batch_name = 'test-batch'

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

    // .collect() operator gathers all output items from the channel into a single list
    collectGreetings(convertToUpper.out.collect(), params.batch_name)  // use the output files from convertToUpper as input to collectGreetings
    // view operators to print the channel contents before and after collect()
    convertToUpper.out.view{ output_files -> "Before collect: $output_files" }
    convertToUpper.out.collect().view{ output_files -> "After collect: $output_files" }
}
