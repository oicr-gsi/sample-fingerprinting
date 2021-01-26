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
    }

    parameter_meta {
        finList: "input .fin file"
    }

    command <<<
        module load sample-fingerprinting/0
        jaccard_coeff_matrix_mc --list ~{finList} > matrix.txt
    >>>

    output {
        File fin = "matrix.txt"
    }
}
