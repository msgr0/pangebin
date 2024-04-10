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


include { PREPROCESS            } from "./workflows/preprocess"
include { BLAST                 } from "./workflows/blast"
include { PANGENOME             } from "./workflows/pangenome"
include { PLASBIN               } from "./workflows/plasbinflow.nf"
include { PLASEVAL              } from "./workflows/plaseval.nf"
// include { PLASBIN as PFLOW_PAN    } from "./workflows/plasbinflow.nf"
// include { PLASBIN as PFLOW_ASM    } from "./workflows/plasbinflow.nf"
include { GPLAS as GPLAS_PAN    } from "./workflows/gplas.nf"
include { GPLAS as GPLAS_ASM      } from "./workflows/gplas.nf"
// include { PLOTS                 } from "./workflows/plots.nf"


workflow PANGEBIN {
    take:
    pangenome
    fasta

    main:
    PLASBIN(pangenome, "pangenome", 1000, 0.5, 1500)
    PLASEVAL(PLASBIN.out.bins, pangenome, fasta, "pangenome")

    emit:
    res = PLASBIN.out.bins
    eval_0 = PLASEVAL.out.mode0
    eval_1 = PLASEVAL.out.mode1


}

workflow ASSEMBLER_A {
    take:
    pangenome
    asm
    fasta

    main:
    PLASBIN(asm, "unicycler", 2000, 0.7, 1500)
    PLASEVAL(PLASBIN.out.bins, pangenome, fasta, "unicycler")

    emit:
    res = PLASBIN.out.bins
    eval_0 = PLASEVAL.out.mode0

}

workflow ASSEMBLER_B {
    take:
    pangenome
    asm
    fasta

    main:
    PLASBIN(asm, "skesa", 2000, 0.7, 1500)
    PLASEVAL(PLASBIN.out.bins, pangenome, fasta, "skesa")

    emit:
    res = PLASBIN.out.bins
    eval_0 = PLASEVAL.out.mode0


}




workflow {
    // GET INPUT FOLDER WITH SAMPLES, assemlby one and two
    if (params.dataset != null ) {
        samples = Channel.fromPath("$params.dataset/*", type: 'dir')
    } else if (params.input != null) {
        samples = Channel.fromPath("$params.input", type: 'dir')
    }

    input_ch = samples.map{dir -> tuple([id: dir.getName()], file("$dir/${ass_A_gfa_gz}"), file("$dir/${ass_B_gfa_gz}"))}

    input_ch | view

    PREPROCESS(input_ch)
    PANGENOME(
        PREPROCESS.out.assemblies,
        PREPROCESS.out.pangenome
    )
    
    // PANGEBIN(
    //     PANGENOME.out.augmented,
    //     PREPROCESS.out.pangenome.map{meta, fasta_AB, fasta_AB_gz -> tuple(meta, fasta_AB)}
    // )
    
    // ASSEMBLER_A(
    //     PANGENOME.out.augmented,
    //     PREPROCESS.out.assemblies.map{it, A, B -> tuple(it, A)},
    //     PREPROCESS.out.fastas.map{it, A, B -> tuple(it, A)}
    // )

    // ASSEMBLER_B(
    //     PANGENOME.out.augmented,
    //     PREPROCESS.out.assemblies.map{it, A, B -> tuple(it, B)},
    //     PREPROCESS.out.fastas.map{it, A, B -> tuple(it, B)}
    // )
    
    // PANGEBIN.out.res.map{meta, bins -> bins} | view
    // PANGEBIN.out.eval_0 | view
    // PANGEBIN.out.eval_1 | view

    // ASSEMBLER_A.out.res.map{meta, bins -> bins} | view
    // ASSEMBLER_A.out.eval_0 | view

    // ASSEMBLER_B.out.res.map{meta, bins -> bins} | view
    // ASSEMBLER_B.out.eval_0 | view


    GPLAS_PAN(
        PANGENOME.out.augmented
    )

    GPLAS_ASM(
        PREPROCESS.out.assemblies.map{it, A, B -> tuple(it, A)}
    )

    GPLAS_PAN.out.bins | view
    GPLAS_ASM.out.bins | view


}
