// input: pan-assembly graph  + mlplasmid prediction/plasmid database

include { PUBLISH } from "./utils"



process computeScores {
  cache 'lenient'
  input:
  tuple val(meta), path(gfa), path(fasta)
  val(mode)

  output:
  tuple val(meta), path(plscore), path(gc)
  

  script:
  if (mode == "pan") {
    asm = 'pan'
  }
  else if( mode == 'asm') {
    asm = meta.asm
  }
  base = "${meta.id}.${asm}.${meta.thr}"
  input_csv = base + ".input.csv"
  plscore = "${meta.id}.${asm}.gd.tsv"
  gc = "${meta.id}.${asm}.gc.tsv"
  putils = "$projectDir/PlasBin-flow-pangenome/code/plasbin_utils.py"
  pls_db_file = "$projectDir/PlasBin-flow-pangenome/database/genes.fasta"

  """
  #!/bin/bash
  echo "sample,gfa" > "${input_csv}"
  echo "${meta.id}.${asm},${gfa}.gz" >> "${input_csv}"
  
  bgzip -k ${gfa}

  python ${putils} preprocessing --input_file ${input_csv} --out_dir . --tmp_dir ./tmp --out_file out_file.csv --db_file ${pls_db_file}
  """
}


process model {
    cache 'lenient'
    maxForks 1

    input:
    tuple val(meta), path(gfa), path(plscore), path(gc)
    val(mode)
    
    output:
    tuple val(meta), path(res), emit: res
    tuple val(meta), path(bins), emit: bins
    tuple val(meta), path(bins_mod1), emit: mod1

    script:

    if (mode == "pan") {
      asm = "pan"
      args = ""
      seed_len = 1000
      seed_score = 0.5
      min_pls_len = 1000
    }
    else if (mode == "asm") {
      asm = meta.asm
      args = "--nopenalty"
      seed_len = 2650
      seed_score = 0.58
      min_pls_len = 3000
    }


    base = "${meta.id}.${asm}.${meta.thr}"
    bins = base + ".pbf.bins.tsv"
    res = base + ".pbf.pred.tab"
    pflow = "$projectDir/PlasBin-flow-pangenome/code/plasbin_flow.py"
    mod1 = base + 'pbf.mod1.tab'

    """
    bgzip -k ${gfa}
    python ${pflow} ${args} -alpha4 1 -ag ${gfa}.gz -gc ${gc} -out_dir . -out_file ${bins}  -score ${plscore} -assembler pangenome -seed_len ${seed_len}  -seed_score ${seed_score} -min_pls_len ${min_pls_len}
    python $projectDir/bin/evaluation/transform_pbf_pred.py --input ${bins} --gfa ${gfa} --output ${res} 

    """
}

// process extendBins {
//   input

//   output

  
//   script:
//    """
//     python $projectDir/bin/extend_bins.py --pred ${res} --out ${mod1} --naive -n 1 ${gfa}
// """
// }

workflow MODEL {
  
  take:
  mode
  gfa_ch
  fasta_ch

  main:

  scores_ch = computeScores(gfa_ch.join(fasta_ch), mode)
  bins = model(gfa_ch.join(scores_ch), mode)
  // extendBins(bins)

  emit:
  res = model.out.res
  bins = model.out.bins

}
