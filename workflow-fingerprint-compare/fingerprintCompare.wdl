version 1.0

workflow fingerprintCompare {

    input {
        File finList
    }

    call generateMatrix { input: finList = finList }

    parameter_meta {
        finList: "Input .txt file containing list of .fin files"
    }

    meta {
        author: "Michelle Feng"
        email: "mfeng@oicr.on.ca"
        description: "Workflow to generate jaccard matrices for all projects."
        dependencies: [
            {
            name: "perl/5.30",
            url: "https://www.perl.org/"
            }
        ]
        output_meta: {
            output: "generated jaccard matrix"
        }
    }

    output {
        File output = generateMatrix.fin
    }
}

task generateMatrix {

    input {
        File finList
        String modules = "sample-fingerprinting/0"
        Int timeout = 12
        Int memory = 21
    }

    parameter_meta {
        finList: "The .txt file containing a list of .fin files on which a jaccard matrix will be generated."
    }

    command <<<
        jaccard_coeff_matrix_mc --list ~{finList} > matrix.txt
    >>>

    runtime {
        modules: "~{modules}"
        timeout: "~{timeout}"
        memory: "~{memory}G"
    }

    output {
        File fin = "matrix.txt"
    }

    meta {
        output_meta: {
            fin: "A good description of what an output file is:"
        }
    }
}
