workflows = "$projectDir/workflows"
src_dir = "$projectDir/bin"
ass_A = "short"
ass_B = "skesa"
pangenome = "pangenome"
ass_A_type = "uni"
ass_B_type = "ske"
ass_A_gfa = "short.gfa"
ass_B_gfa = "skesa.gfa"
ass_A_gfa_gz = "short.gfa.gz"
ass_B_gfa_gz = "skesa.gfa.gz"
pan_type = "pan"



if (params.input == null & params.dataset == null) {
    println "Please provide input folder or dataset folder"
    System.exit(1)
}


include { FASTA_PREPROCESS      } from "./workflows/fasta_preprocess"
include { BLAST                 } from "./workflows/blast"
include { PANGENOME             } from "./workflows/pangenome"
include { PLASBIN as PFLOW_PAN  } from "./workflows/plasbinflow.nf"
include { PLASBIN as PFLOW_ASM    } from "./workflows/plasbinflow.nf"
include { GPLAS as GPLAS_PAN    } from "./workflows/gplas.nf"
include { GPLAS as GPLAS_ASM      } from "./workflows/gplas.nf"
include { PLOTS                 } from "./workflows/plots.nf"


workflow PANGEBIN {
    take:
    input

    main:

    FASTA_PREPROCESS(input)
    PANGENOME(
        FASTA_PREPROCESS.out.id
        FASTA_PREPROCESS.out.fasta_A
        FASTA_PREPROCESS.out.fasta_B
    )

    PFLOW_PAN(PANGENOME.out.augmented)

    RESULTS(PFLOW_PAN.out)
    
    emit:
    bins = PFLOW_PAN.out.bins
    res = RESULTS.out

}

workflow ASSEMBLER_A {
    take:
    input
    
    main:

    FASTA_PREPROCESS(input)
    PFLOW_ASM(FASTA_PREPROCESS.out.gfa_A_gz)

    RESULTS(PFLOW_ASM.out)

    emit:
    RESULTS.out

}

workflow ASSEMBLER_B {
    take:
    input
    
    main:

    FASTA_PREPROCESS(input)
    PFLOW_ASM(FASTA_PREPROCESS.out.gfa_B_gz)

    RESULTS(PFLOW_ASM.out)

    emit:
    RESULTS.out

}

workflow GPLAS {

    take:
    input
    
    main:
    FASTA_PREPROCESS(input)
    PANGENOME(
        FASTA_PREPROCESS.out.id
        FASTA_PREPROCESS.out.fasta_A
        FASTA_PREPROCESS.out.fasta_B
    )


    GPLAS_PAN(PANGENOME.out.augmented)
    GPLAS_ASM(FASTA_PREPROCESS.out.gfa_A_gz)
}


work


workflow {
    // GET INPUT FOLDER WITH SAMPLES, assemlby one and two
    if (params.dataset != null ) {
        samples = Channel.fromPath("$params.dataset/*", type: 'dir')
    } else if (params.input != null) {
        samples = Channel.fromPath("$params.input", type: 'dir')
    }

    input = samples.map{dir -> [dir.getName(), file("$dir/${ass_A_gfa_gz}"), file("$dir/${ass_B_gfa_gz}")]}
    
    PANGEBIN(input)

    PANGEBIN.out.bins | view

    // ASSEMBLER_A(input)

    // ASSEMBLER_B(input)

    // GPLAS(input)

}
