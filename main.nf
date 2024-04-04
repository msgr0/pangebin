workflows = "$projectDir/workflows"
src_dir = "$projectDir/bin"
test_sample = "ecol-SAMN04014855"
test_dataset = "$projectDir/test_data"
ass_A = "short"
ass_B = "skesa"
ass_A_type = "uni"
ass_B_type = "ske"

params.input = ""
src_dir = "$projectDir/src"
test_data = "$projectDir/test_data"

if (params.input == null & params.dataset == null) {
    println "Please provide input folder or dataset folder"
    System.exit(1)
}


include {FASTA_PREPROCESS} from "./workflows/fasta_preprocess"
include {BLAST} from "./workflows/blast"
include {PANGENOME} from "./workflows/pangenome"
include {CHECK_INPUT} from "./workflows/check_input"
include {PLASBIN as PLASBIN_FLOW_PANGENOME} from "./workflows/plasbinflow.nf"
include {PLASBIN as PLASBIN_FLOW_ASSEMBLY_ONE} from "./workflows/plasbinflow.nf"
include {PLASBIN as PLASBIN_FLOW_ASSEMBLY_TWO} from "./workflows/plasbinflow.nf"



workflow PIPE_PANGENOME {
    FASTA_PREPROCESS()
    BLAST()
    PANGENOME()

}

workflow PIPE_COMPARE {
    FASTA_PREPROCESS()
    BLAST()
    PANGENOME()
    PLASBIN_FLOW_PANGENOME()
    PLASBIN_FLOW_ASSEMBLY()
    PLASBIN_FLOW_ASSEMBLY()
}

workflow {
    
    CHECK_INPUT()


    // INPUT FOLDER WITH SAMPLES
    samples = Channel.fromPath("")
    if (params.pangenome)  // pangenome only
    PIPE_PANGENOME()
    else if (params.compare) // comparision
    PIPE_COMPARE()
}





workflow {
    // sample_ids = Channel.of(test_sample)
    sample_folders = Channel.fromPath("${test_dataset}/*", type: 'dir')
    sample_ids = sample_folders.map{it-> it.getName()}
    genomic_files = ncbi_retrive(sample_ids)

    fastas = sample_folders.map{dir ->  [dir.getName(), file("$dir/${ass_A}.gfa.gz"), file("$dir/${ass_B}.gfa.gz")]}
    fastas_prep = fastas_preprocess(fastas)
    blast_info = blast(fastas_prep.map{it -> [it[0], it[-1]]}.join(genomic_files))
    pangenome = make_pangenome(fastas_prep.map{it->[it[0], it[-2]]})

    
     
}