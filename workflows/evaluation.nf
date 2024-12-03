#!/usr/bin/env nextflow

/*
EVALUATION pre-preprocessing
. unicycler and skesa bins
. sample name -> to retrive ground truth
. compute unicycler and skesa ground truth LABELING and BINS
. evaluate unicycler and skesa labeling and binning

. choose MLPLASMID or PLSDENSITY
. panassembly and panassembly-100 bins
. computa unicycler and skesa gt LABELING and BINS in respect of PANASSEMBLY fragments
. compute unicycler and skesa -100 gt LABELING and BINS in respect of PANASSEMBLY-100 fragments
. evaluate pan-assembly and pan-assembly 100 labeling and binning against (unicycler,skesa) and (uniycler,skesa -100)


LABELING:
. evaluate precision recall and f1 of how many fragments/contigs * lengths have been correctly evaluated as plasmid

BINNING:
. evaluate precision recall and f1 using PLASEVAL and as describe in plaseval, using the cumulative-length value

*/

/*

current evaluation
MLPLASMID + pan-assembly/pan-assembly-100 vs Unicycler vs Skesa with super-bins

next:
plasmidness, extended+super bins
*/


process labeling {
    cache false

    input:
    tuple val(meta), path(prediction), path(gt)

    output:
    tuple val(meta), path(stats)

    script:

    output = prediction.getBaseName()
    id = output.split("\\.")[0]
    asm = output.split("\\.")[1]
    thr = output.split("\\.")[2]
    tool = output.split("\\.")[3]

    reference = gt.getBaseName().split("\\.")[2]

    description = "sample ${id}, ${asm} graph (cut ${thr}, ref ${reference}) with ${tool}"
    stats = "${output}.${reference}.lab.txt"

    """
    python $projectDir/bin/evaluate_bins.py --bin ${prediction} --csv ${gt} --sample ${output} --output ${stats} --description '${description}'
    """
}


process binning { 
    cache false

    input:
    tuple val(meta), path(prediction), path(gt)

    output:
    tuple val(meta), path(stats)
    // tuple val(meta), path(plots), emit: plots

    script:
    output = prediction.getBaseName()
    reference = gt.getBaseName().split("\\.")[2]
    stats = "${output}.${reference}.bin.txt"
    // awk -i inplace '{\$0=gensub(/\s*\S+/,\\"\\",3)}1' ${gt} 
    // awk -i inplace '{\$0=gensub(/\s*\S+/,\\"\\",4)}1' ${gt} 
    """
    #!/bin/bash
    echo 'plasmid\tcontig\tnan\tnan\tcontig_len' | cat - ${gt} > temp && mv temp ${gt}
    python $projectDir/PlasEval/src/plaseval.py eval --pred ${prediction} --gt ${gt} --out ${stats} 
    """

}

workflow EVALUATION {
    take:
    input_ch

    main:
    input_ch
    labeling_ch = labeling(input_ch)
    binning_ch = binning(input_ch)

    emit:
    
    labeling = labeling_ch
    binning = binning_ch
}
  
