version 1.0

workflow fingerprintCompare {

    input {
        File finFile
    }

    call generateMatrix { input: finFile = finFile }

    parameter_meta {
        finFile: "Input .fin file containing list of .fin files"
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
            outputFin: "generated jaccard matrix"
        }
    }

    output {
        File outputFin = generateMatrix.fin
    }
}

task generateMatrix {

    input {
        File finFile
    }

    parameter_meta {
        finFile: "input .fin file"
    }

    command <<<
        module load sample-fingerprinting/0
        jaccard_coeff_matrix_mc --list ~{finFile} > matrix.txt
    >>>

    output {
        File fin = "matrix.txt"
    }
}
