#!/usr/bin/env nextflow

def convert_species(spec) {
  switch(spec) {
    case 'efea':
      species = 'Enterococcus faecium'
      break
    case 'kpne':
      species = 'Klebsiella pneumoniae'
      break
    case 'abau':
      species = 'Acinetobacter baumannii'
      break
    case 'ecol':
      species = 'Escherichia coli'
      break
    default:
      species = 'Other species'
  }
  return species
}

process MLPLASMIDS {
    // errorStrategy 'ignore'


    input:
    tuple val(meta), path(mixed), path(uni), path(ske)

    output:
    tuple val(meta), path(mixedpred), emit: mixed
    tuple val(meta), path(unipred), emit: uni
    tuple val(meta), path(skepred), emit: ske

    
    script:
    mixedpred = "${meta.id}.mlplasmid.mix.pred"
    unipred = "${meta.id}.mlplasmid.uni.pred"
    skepred = "${meta.id}.mlplasmid.ske.pred"

    mlplas_threshold = '0.5'
    species = convert_species("${meta.species}")
    """
    #!/bin/bash

      Rscript $projectDir/bin/run_mlplasmids.R ${mixed} ${mixedpred} ${mlplas_threshold} '${species}' TRUE
      Rscript $projectDir/bin/run_mlplasmids.R ${uni} ${unipred} ${mlplas_threshold} '${species}' TRUE
      Rscript $projectDir/bin/run_mlplasmids.R ${ske} ${skepred} ${mlplas_threshold} '${species}' TRUE
    """
}