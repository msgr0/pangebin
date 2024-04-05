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
include { PLASBIN as PFLOW_A    } from "./workflows/plasbinflow.nf"
include { PLASBIN as PFLOW_B    } from "./workflows/plasbinflow.nf"
include { GPLAS as GPLAS_PAN    } from "./workflows/gplas.nf"
include { GPLAS as GPLAS_A      } from "./workflows/gplas.nf"
include { PLOTS                 } from "./workflows/plots.nf"


workflow PANGEBIN {
    take:
    input

    main:
    FASTA_PREPROCESS(m)
    PANGENOME(
        FASTA_PREPROCESS.out.id
        FASTA_PREPROCESS.out.fasta_A
        FASTA_PREPROCESS.out.fasta_B
    )

    PFLOW_PAN(PANGENOME.out.augmented)
    GPLAS_PAN(PANGENOME.out.augmented)
    GPLAS_A(input.gfa_A_gz)
    PFLOW_A(input.gfa_A_gz)
    PFLOW_B(input.gfa_B_gz)

    PLOT(PFLOW_PAN.out.results, GPLAS_PAN.out.results)
    
    emit:
    plot = PLOT.out
    pflow = PFLOW_PAN.out.bin
    gplas = GPLAS_PAN.out.bin
}

workflow {
    // GET INPUT FOLDER WITH SAMPLES, assemlby one and two
    if (params.dataset != null ) {
        samples = Channel.fromPath("$params.dataset/*", type: 'dir')
    } else if (params.input != null) {
        samples = Channel.fromPath("$params.input", type: 'dir')
    }

    input = samples.map{dir -> [dir.getName(), file("$dir/${ass_A_gfa_gz}"), file("$dir/${ass_B_gfa_gz}")]}
    
    PANGEBIN(input)

    PANGEBIN.pangebin | view

}





// workflow {
//     // sample_ids = Channel.of(test_sample)
//     sample_folders = Channel.fromPath("${test_dataset}/*", type: 'dir')
//     sample_ids = sample_folders.map{it-> it.getName()}
//     genomic_files = ncbi_retrive(sample_ids)

//     fastas = sample_folders.map{dir ->  [dir.getName(), file("$dir/${ass_A}.gfa.gz"), file("$dir/${ass_B}.gfa.gz")]}
//     fastas_prep = fastas_preprocess(fastas)
//     blast_info = blast(fastas_prep.map{it -> [it[0], it[-1]]}.join(genomic_files))
//     pangenome = make_pangenome(fastas_prep.map{it->[it[0], it[-2]]})

    
     
// }