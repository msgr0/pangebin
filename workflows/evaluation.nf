#!/usr/bin/env nextflow


workflow EVALUATE {
    take:
    input

    main:
    EVAL(input)
    inpath = Channel.fromPath('${params.input}', type: 'dir')

    emit:
    stats = EVAL.out.stats

}

process EVAL {

    input:
    tuple val(meta), path(prediction), path(gt), val(tool)

    output:
    tuple val(meta), path(stats), emit: stats
    // tuple val(meta), path(plots), emit: plots

    script:

    output = "${meta.id}.${tool}"
    description = "${meta.id}.${meta.thr} with ${tool}"
    // plots = "${output}.scores.pdf"
    stats = "${output}.stats.txt"

    """
    python $projectDir/bin/evaluate_bins.py --bin ${prediction} --csv ${gt} --sample ${meta.id} --output ${output} --description '${description}'
    """
}