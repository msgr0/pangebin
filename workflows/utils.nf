#!/usr/bin/env nextflow

process PUBLISH {
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: "${item}"
    // , saveAs: {-> it[0..-3] + "." + ${meta.thr} + it[-3..-1]}
    //, saveAs: {-> it[0..-1]} 
    // + it[-3..-1]}
    
    input:
    tuple val(meta), path(item)
    
    
    output:
    tuple val(meta), path(item)

    """
    echo
    """
}
