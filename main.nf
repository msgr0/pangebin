workflows = "$projectDir/workflows"
src_dir = "$projectDir/bin"

nameA = "short"
nameB = "skesa"
nameP = "pangenome"
typeA = "uni"
typeB = "ske"
typeP = "pan"

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

params.db = null

// if (params.input == null & params.dataset == null) {
//     println "Please provide input folder or dataset folder"
//     System.exit(1)
// }


include { PREPROCESS            } from "./workflows/preprocess"
// include { BLAST                 } from "./workflows/blast"
include { PANGENOME             } from "./workflows/pangenome"
include { PLASBIN               } from "./workflows/plasbinflow.nf"
// include { PLASEVAL              } from "./workflows/plaseval.nf"
// include { PLASBIN as PFLOW_PAN    } from "./workflows/plasbinflow.nf"
// include { PLASBIN as PFLOW_ASM    } from "./workflows/plasbinflow.nf"
// include { GPLAS as GPLAS_PAN    } from "./workflows/gplas.nf"
// include { GPLAS as GPLAS_ASM      } from "./workflows/gplas.nf"
// include { PLOTS                 } from "./workflows/plots.nf"


// workflow ASSEMBLER_A { //runs plasbin-flow on Unicyler
//     take:
//     pangenome
//     asm
//     fasta

//     main:
//     PLASBIN(asm, "unicycler", 2000, 0.7, 1500)
//     PLASEVAL(PLASBIN.out.bins, pangenome, fasta, "unicycler")

//     emit:
//     res = PLASBIN.out.bins
//     eval_0 = PLASEVAL.out.mode0

// }

// workflow ASSEMBLER_B { //runs plasbin-flow on Skesa
//     take:
//     pangenome
//     asm
//     fasta

//     main:
//     PLASBIN(asm, "skesa", 2000, 0.7, 1500)
//     PLASEVAL(PLASBIN.out.bins, pangenome, fasta, "skesa")

//     emit:
//     res = PLASBIN.out.bins
//     eval_0 = PLASEVAL.out.mode0


// }
process PUBLISH {
    publishDir "${params.out}", mode: 'copy', overwrite: true, pattern: "${name}", saveAs: {filename -> 
            if (meta.tool == "gplas") { "gplas.${filename}" }
            else if (meta.tool == "gplas2") { "gplas2.${filename}" }
            else if (meta.tool = "pangebin") { "pangebin.${filename}"}
            else if (meta.tool = "pbfu") { "pbf-u.${filename}"}
            else if (meta.tool = "pbfs") { "pbf-s.${filename}"}
            else { "other/${filename}" }
        }
    
    input:
    tuple val(meta), path(item)

    output:
    tuple val(meta), path(item), val(name)
    
    script:
    name = item.getName()
    """
    """
}



workflow PANGEBIN {
    take:
    pangenome
    // fasta

    main:
    PLASBIN(pangenome, "pangenome", 1000, 0.5, 1500)
    // PLASEVAL(PLASBIN.out.bins, pangenome, fasta, "pangenome")

    emit:
    res = PLASBIN.out.bins
    // eval_0 = PLASEVAL.out.mode0
    // eval_1 = PLASEVAL.out.mode1

}


workflow {
    // GET INPUT FOLDER WITH SAMPLES, assemlby one and two
    // if (params.dataset != null ) {
    //     samples = Channel.fromPath("$params.dataset/*", type: 'dir')
    // } else if (params.input != null) {
    //     samples = Channel.fromPath("$params.input", type: 'dir')
    // }

    samples_u = Channel.fromPath("$params.db/robertson-benchmark_*-*-u.gfa.gz", type: 'file')
    samples_s = Channel.fromPath("$params.db/robertson-benchmark_*-*-s.gfa.gz",
    type: 'file')

    
    input_s = samples_s.map{file -> tuple([id: file.getName()[20..36], species: file.getName()[20..23]], file)}
    input_u = samples_u.map{file -> tuple([id: file.getName()[20..36], species:
    file.getName()[20..23]], file)}
    input_ch = input_u.join(input_s)
    // input_ch | view

    PREPROCESS(input_ch)
    PANGENOME(
        PREPROCESS.out.assemblies,
        PREPROCESS.out.pangenome
    )
    
    PANGEBIN(
        PANGENOME.out.augmented//,
        // PREPROCESS.out.pangenome.map{meta, fasta_AB, fasta_AB_gz -> tuple(meta, fasta_AB)}
    )

    PUBLISH(PANGEBIN.out.res)
}
    
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


    // GPLAS_PAN(
    //     PANGENOME.out.augmented
    // )

    // GPLAS_ASM(
    //     PREPROCESS.out.assemblies.map{it, A, B -> tuple(it, A)}
    // )

    // GPLAS_PAN.out.bins | view
    // GPLAS_ASM.out.bins | view