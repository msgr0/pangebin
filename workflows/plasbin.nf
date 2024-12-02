#!/usr/bin/env nextflow

process model {
    
    cache 'lenient'
    maxForks 1

    input:
    tuple val(meta), path(gfa), path(gc), path(plscore)
    val(mode)
    
    output:
    tuple val(meta), path(res), emit: res
    tuple val(meta), path(bins), emit: bins
    // tuple val(meta), path(bins_mod1), emit: mod1

    script:

    if (mode == "pan") {
        assembler = "pangenome"
        asm = "pan"
        args = ""
        seed_len = params.pangenome_seedlen
        seed_score = params.pangenome_seedscore
        min_pls_len = params.pangenome_minplaslen
    }

    else if (mode == "uni" || mode == "ske") {
        asm = mode
        assembler = "unicycler"
        args = "--nopenalty"
        seed_len = params.assembly_seedlen
        seed_score = params.assembly_seedscore
        min_pls_len = params.assembly_minplaslen
    }

    base = "${meta.id}.${mode}.${meta.thr}"
    bins = base + ".pbf.bins.tsv"
    res = base + ".pbf.pred.tab"
    pflow = "$projectDir/PlasBin-flow-pangenome/code/plasbin_flow.py"
    // mod1 = base + 'pbf.mod1.tab'

    """
    bgzip -k ${gfa}
    python ${pflow} ${args} -alpha4 1 -ag ${gfa}.gz -gc ${gc} -out_dir . -out_file ${bins}  -score ${plscore} -assembler ${assembler} -seed_len ${seed_len}  -seed_score ${seed_score} -min_pls_len ${min_pls_len}
    python $projectDir/bin/evaluation/transform_pbf_pred.py --input ${bins} --gfa ${gfa} --output ${res} 

    """
}

// process MOD_BINS {
//     cache false

//     input:
//     tuple val(meta), path(pred), path(graph)
//     val(flag)
//     val(name)

//     output:
//     tuple val(meta), path(pred_mod)

//     script:
//     pred_mod = pred.baseName + ".${name}.tsv"

//     """
//     python $projectDir/bin/extend_bins.py --pred ${pred} --out ${pred_mod} ${flag} ${graph}
//     """

// }

// process extendBins {
//     input:

//     output:

  
//     script:
//     """
//     python $projectDir/bin/extend_bins.py --pred ${res} --out ${mod1} --naive -n 1 ${gfa}
//   	"""
// }

// process MOD_BINS {
//     cache false

//     input:
//     tuple val(meta), path(pred), path(graph)
//     val(flag)
//     val(name)

//     output:
//     tuple val(meta), path(pred_mod)

//     script:
//     pred_mod = pred.baseName + ".${name}.tsv"

//     """
//     python $projectDir/bin/extend_bins.py --pred ${pred} --out ${pred_mod} ${flag} ${graph}
//     """

// }

workflow MODEL {
  
	take:

	mode
	gfa_ch

	gc_ch
	gd_ch

	main:

	bins_ch = model(gfa_ch.join(gc_ch).join(gd_ch), mode)
    
    // naiveBins_ch        = Channel.empty()
    // graphOverlapBins_ch = Channel.empty()

    // if (mode == "pan") {
    //     naiveBins_ch        = extend_naive(bins_ch)
    //     graphOverlapBins_ch = extend_overlap(bins_ch)
    // }
    

	emit:
	res = model.out.res
	bins = model.out.bins
    // naiveB = naiveBins_ch
    // overlapB = graphOverlapBins_ch
}
