#!/usr/bin/env nextflow

process PUBLISH {
    publishDir "${params.output}/${outpath}", mode: 'copy', overwrite: true, pattern: "${item}", saveAs: {-> it[0..-3] + ".${outpath}." + it[-3..-1]}
    
    input:
    tuple val(meta), path(item)
    path(outpath)
    
    output:
    tuple val(meta), path(item)

    """
    echo
    """
}
