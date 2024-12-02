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
    cache 'lenient'

    input:
    tuple val(meta), path(gfa), path(fasta)

    output:
    tuple val(meta), path(pred_pbf), emit: pbf
    
    script:
    mlplas_threshold = '0.5'
    species = convert_species("${meta.species}")
    pred = fasta.baseName + "mlplas.tab"
    """
    #!/bin/bash
    Rscript $projectDir/bin/run_mlplasmids.R ${fasta} ${pred} ${mlplas_threshold} '${species}' TRUE
    python $projectDir/bin/mlpl.asm.py --pred ${pred} --graph ${gfa}  --output ${pred_gplas}
    python $projectDir/bin/mlpl.asm.py --pred ${pred} --graph ${gfa}  --pbf ${pred_pbf}
    """
}