#!/usr/bin/env nextflow

workflow CHECK_INPUT {
    if (params.dataset != null) {
        sample_directories = Channel.fromPath("$params.dataset/*", type: 'dir')
    }
    else if (params.input != null){
        sample_directories = Channel.fromPath("$params.input", type: 'dir')
    }
    else if (params.plasgraph != null){
        sample_files = Channel.fromFilePairs("$params.plasgraph/*-{u,s}.gfa.gz", flat: true)
    }
    else {
        println "Please provide input folder or dataset folder"
        System.exit(1)
    }
}