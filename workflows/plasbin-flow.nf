// input: unicycler & skesa graph  + mlplasmid prediction/plasmid database

include { MLPLASMIDS } from "./mlplasmids"
include { PUBLISH } from "./utils"

binfolder = "~/bin"
projectDir = "~/pangebin"
putils = "$projectDir/PlasBin-flow-pangenome/code/plasbin_utils.py"
pflow = "$projectDir/PlasBin-flow-pangenome/code/plasbin_flow.py"
pdatabase = "$projectDir/Plasbin-flow-pangenome/database/genes.fasta"


seed_len = 2000
seed_score = 0.7
min_pls_len = 1500

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
      species = 'Escherichia coli'
  }
  return species
}

process MLPLAS_PRED_ASM {
    cache 'lenient'

    input:
    tuple val(meta), path(gfa), path(fasta)
    
    output:
    // tuple val(meta), path(pred_gplas), emit: gplas
    tuple val(meta), path(pred_pbf), emit: pbf

    script:
    mlplas_threshold = '0.5'
    species = convert_species("${meta.species}")

    pred = fasta.baseName + "mlplas.tab"
    pred_pbf = fasta.baseName + ".mlplas.pbf.tab"


    """
    #!/bin/bash
    Rscript $projectDir/bin/run_mlplasmids.R ${fasta} ${pred} ${mlplas_threshold} '${species}' TRUE
    python $projectDir/bin/mlpl.asm.py --pred ${pred} --graph ${gfa}  --pbf ${pred_pbf}

    """
    // python $projectDir/bin/mlpl.asm.py --pred ${pred} --graph ${gfa}  --output ${pred_gplas}
}

// process MLPLAS_PRED_PAN {
//     cache 'lenient'
//     errorStrategy 'ignore'
    
//     input:
//     tuple val(meta), path(gfa), path(fasta)
    
//     output:
//     // tuple val(meta), path(pred_gplas), emit: gplas
//     tuple val(meta), path(pred_pbf), emit: pbf

//     script:
//     mlplas_threshold = '0.5'
//     species = convert_species("${meta.species}")

//     pred = fasta.baseName + ".mlplas.tab"
//     // pred_gplas = fasta.baseName + ".mlplas.gplas.tab"
//     pred_pbf = fasta.baseName + ".mlplas.pbf.tab"
    

//     """
//     #!/bin/bash
//     Rscript $projectDir/bin/run_mlplasmids.R ${fasta} ${pred} ${mlplas_threshold} '${species}' TRUE
//     python $projectDir/bin/mlpl.asmtopan.py --pred ${pred} --graph ${gfa} --pbf ${pred_pbf}

//     """
//     // python $projectDir/bin/mlpl.asmtopan.py --pred ${pred} --graph ${gfa} --output ${pred_gplas}
// }

// process PLASMIDNESS {
//     cache 'lenient'

//     input:
//     tuple val(meta), path(gfa), path(fasta)
    
//     output:
//     tuple val(meta), path(plscore), path(gc_content) emit: pred

//     script:
//     base = "${meta.species}.${meta.id}${asm}.${meta.thr}"
//     input_csv = base + ".input.csv"
//     plscore = "${meta.id}.${asm}.gd.tsv"
//     gc_content = "${meta.id}.${asm}.gc.tsv"
//     pls_db_file = "$projectDir/Plasbin-flow-pangenome/database/genes.fasta"


//     """
//     #!/bin/bash
//     echo "sample,gfa" > "${input_csv}"
//     echo "${meta.id}.${asm},${gfa}.gz" >> "${input_csv}"
//     python plasbin_utils.py preprocessing --input_file ${input_csv} --out_dir . --tmp_dir ./tmp --out_file out_file.csv --db_file ${pls_db_file}
//     """
// }

// process GENE_DB {
// 	cache 'lenient'

// 	input:
// 	tuple val(meta), path()

// 	output:

// 	script:
// 	"""
// 	#!/bin/bash
// 	"""
// }

// process PBF {
//     cache 'lenient'
//     maxForks 4

//     input:
//     tuple val(meta), path(gfa), path(plscore), path(gc)
    
//     output:
//     tuple val(meta), path(res), emit: res
//     tuple val(meta), path(bins), emit: bins

//     script:
//     asm = ""

//     base = "${meta.species}.${meta.id}${asm}.${meta.thr}"
//     bins = base + ".pbf.bins.tsv"
//     res = base + ".pbf.pred.tab"

//     """
//     bgzip -k ${gfa}
//     # --nopenalty
//     python ${pflow} -alpha4 1 -ag ${gfa}.gz -gc ${gc} -out_dir . -out_file ${bins}  -score ${pred} -assembler pangenome -seed_len ${seed_len}  -seed_score ${seed_score} -min_pls_len ${min_pls_len}
//     python $projectDir/bin/evaluation/transform_pbf_pred.py --input ${bins} --gfa ${gfa} --output ${res} 
//     """

// }

process PBF_MLPLAS {
    cache 'lenient'
    errorStrategy 'ignore'
    maxForks 4

    input:
    tuple val(meta), path(gfa), path(pred)

    output:
    tuple val(meta), path(res), emit: res
    tuple val(meta), path(bins), emit: bins

    script:

    asm = ""

    base = "${meta.species}.${meta.id}${asm}.${meta.thr}"
    input_csv = base + ".input.csv"
    gc_content = "${meta.id}.${asm}.gc.tsv"
    bins = base + ".pbf.bins.tsv"
    res = base + ".pbf.pred.tab"

    """
    bgzip -k ${gfa}
    echo "sample,gfa,pls_score" > "${input_csv}"
    echo "${meta.id}.${asm},${gfa}.gz,${pred}" >> "${input_csv}"
    python ${putils} gc_probabilities --input_file ${input_csv} --out_dir . --tmp_dir ./tmp

    python ${pflow} --nopenalty -alpha4 1 -ag ${gfa}.gz -gc ${gc_content} -out_dir . -out_file ${bins}  -score ${pred} -assembler unicycler -seed_len ${seed_len}  -seed_score ${seed_score} -min_pls_len ${min_pls_len}  

    python $projectDir/bin/evaluation/transform_pbf_pred.py --input ${bins} --gfa ${gfa} --output ${res} 
    """
}

// workflow PREDICTION {
// 	MLPLAS_PRED_PAS(panassembly_gfa.join(mixed_fasta))
//     mlplas_gplas_ch = MLPLAS_PRED_PAN.out.gplas
//     mlplas_pbf_ch = MLPLAS_PRED_PAN.out.pbf
// }

workflow {

	// panassembly_gfa = Channel.fromPath("${params.input}/*S*panasm.gfa", type: 'file', checkIfExists: true)
	// panassembly_gfa = panassembly_gfa.map{file ->
  //       def fmeta = [:];
  //       fmeta.id = file.toString().split('\\.')[1];
  //       fmeta.species = file.toString().split('\\.')[0].split('/')[-1];
  //       fmeta.thr = file.toString().split('\\.')[2];
  //       [fmeta, file]
  //   }
  //   // panassembly_gfa | view

	// mixed_fasta = Channel.fromPath("${params.input}/*S*\\.fa", type: 'file', checkIfExists: true)
	// mixed_fasta = mixed_fasta.map{file ->
  //       def fmeta = [:];
  //       fmeta.id = file.toString().split('\\.')[1];
  //       fmeta.species = file.toString().split('\\.')[0].split('/')[-1];
  //       fmeta.thr = file.toString().split('\\.')[2];
  //       [fmeta, file]
  //   }
  //   panassembly_gfa_fasta = panassembly_gfa.join(mixed_fasta)

  ske_ch = Channel.fromFilePairs("${params.input}/*-s.0.{gfa,fasta}", type: 'file', checkIfExists: true)
  uni_ch = Channel.fromFilePairs("${params.input}/*-u.0.{gfa,fasta}", type: 'file', checkIfExists: true)


  ske_ch = ske_ch.map{meta, filepair ->
        def fmeta = [:];
        file = filepair[0].toString().split("/")[-1];
        fmeta.id = file.split('-')[1] + '-s';
        fmeta.species = file.split('-')[0];
        fmeta.thr = file.split('-')[2].split('\\.')[1];
        [fmeta, filepair[1], filepair[0]]
  }

  uni_ch = uni_ch.map{meta, filepair ->
        def fmeta = [:];
        file = filepair[0].toString().split("/")[-1];
        fmeta.id = file.split('-')[1] + '-u';
        fmeta.species = file.split('-')[0];
        fmeta.thr = file.split('-')[2].split('\\.')[1];
        [fmeta, filepair[1], filepair[0]]
  }

  // ske_ch | view
  // uni_ch | view


  MLPLAS_PRED_ASM(ske_ch.mix(uni_ch))
  mlplas_pbf_ch = MLPLAS_PRED_ASM.out.pbf
  ske_ch_pred = ske_ch.map{meta, gfa, fasta -> [meta, gfa]}.join(mlplas_pbf_ch)
  uni_ch_pred = uni_ch.map{meta, gfa, fasta -> [meta, gfa]}.join(mlplas_pbf_ch)

  // ske_ch_pred | view
  // uni_ch_pred | view
  // assemblies_ch = ske_ch.map{meta, gfa, fasta -> [meta, gfa]}.join(mlplas_pbf_ch).mix(uni_ch.map{meta, gfa, fasta -> [meta, gfa]}.join(mlplas_pbf_ch))
  PBF_MLPLAS(ske_ch_pred.mix(uni_ch_pred))
  // PBF_MLPLAS(assemblies_ch)
  PUBLISH(PBF_MLPLAS.out.res.mix(PBF_MLPLAS.out.bins))
  
  
  // // unicycler mlplasmid
  // MLPLAS_PRED_ASM()

  // // skesa mlplasmid
  //   // 1
  //   MLPLAS_PRED_PAN(panassembly_gfa_fasta)
  //   mlplas_pbf_ch = MLPLAS_PRED_PAN.out.pbf
  //   PBF_MLPLAS(panassembly_gfa.join(mlplas_pbf_ch))
  //   PUBLISH(PBF_MLPLAS.out.res.mix(PBF_MLPLAS.out.bins))
  //   // 1u
  //   // 1s



  //   // 2
  //   PLASMIDNESS(panassembly_gfa_fasta)
  //   plscores_pbf_ch = PLASMIDNESS.out.pbf
  //   PBF(panassembly_gfa.join(plscores_pbf_ch))
  //   PUBLISH(PBF.out.res.mix(PBF.out.bins))
  //   // 2u
  //   // 2s
    
    
  //   // 1,2 for cut 100 and cut 0
    


}
