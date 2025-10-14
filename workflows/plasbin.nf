#!/usr/bin/env nextflow
def ptylist = [0, 1, 0.5, 2]
def metaname(meta) {
    def name = meta.id
    if (meta.containsKey('species')) {
        name += "-${meta.species}"
    }
    if (meta.containsKey('asm')) {
        name += "-${meta.asm}"
    }
    if (meta.containsKey('cutlen')) {
        name += "-c${meta.cutlen}"
    }
    if (meta.containsKey('pctid')) {
        name += "-${meta.pctid}p"
    }
    if (meta.containsKey('thr')) {
        name += "-t${meta.thr}"
    }
    if (meta.containsKey('pty')) {
        name += "-p${meta.pty}"
    }
    if (meta.containsKey('t')) {
        name += "-${meta.t}"
    }    
    return name
}

process model {
    executor 'slurm'
    cpus 16
    memory '32 GB'    
    cache 'lenient'
    time { 4.hour + 4.hour * task.attempt }
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(meta), path(gfa), path(gc), path(plscore)
    val(type)

    
    output:
    tuple val(meta), path(bins)

    script:
    arg = ""
    if (type == "pan") {
        assembler = "pangenome"
        // asm = "pan"
        if (meta.pty == 0) {
            args = "--nopenalty"
        }
        else {
            args = ""
        }
        seed_len = params.pangenome_seedlen
        seed_score = params.pangenome_seedscore
        min_pls_len = params.pangenome_minplaslen
    }

    else if (type == "asm") {
        assembler = meta.asm == "s" ? "skesa" : "unicycler"
        // asm = meta.asm == "s" ? "ske" : "uni"
        args = "--nopenalty"
        seed_len = params.assembly_seedlen
        seed_score = params.assembly_seedscore
        min_pls_len = params.assembly_minplaslen
    }

    base = metaname(meta)
    bins = "${base}.pbf.bins.tsv"
    pflow = "$projectDir/PlasBin-flow-pangenome/code/plasbin_flow.py"

    """
    bgzip -k ${gfa}
    python ${pflow} ${args} -alpha4 ${meta.pty} -ag ${gfa}.gz -gc ${gc} -out_dir . -out_file ${bins}  -score ${plscore} -assembler ${assembler} -seed_len ${seed_len}  -seed_score ${seed_score} -min_pls_len ${min_pls_len}
    """
}

process transform {
    cache false
    input:
    tuple val(meta), path(bins), path(gfa)

    output:
    tuple val(meta), path(res)
    
    script:

    res = metaname(meta) + ".res.pbf.tab"
    // res = bins.getBaseName()[0..-6] + ".pred.tab"
    
    """
    python $projectDir/bin/evaluation/transform_pbf_pred.py --input ${bins} --gfa ${gfa} --output ${res} 
    """
} 


process modifyBins {
    cache false
    input:

    tuple val(meta), path(pred), path(graph)
    val(mode)

    output:
    tuple val(meta), path(modded)

    script:
    modded = metaname(name)
    if (mode == "naive") {
        modded += ".res.nve.tab"
        """
        #!/usr/bin/env bash
        python $projectDir/bin/extend_bins.py --pred ${pred} --out ${modded} --naive --n 1
        """
    }

    else if (mode == "overlap") {
        modded += ".res.ovl.tab"
        """
        #!/usr/bin/env bash
        python $projectDir/bin/extend_bins.py --pred ${pred} --out ${modded} --super --graph ${graph}
        """
    }        
}

workflow OVERLAP {
    take:

    gfa_ch
    pred_ch

    main:

    res_ch = modifyBins(pred_ch.join(gfa_ch), "overlap")

    emit:
    pred = res_ch
}

workflow NAIVE
{
    take: 

    gfa_ch
    pred_ch

    main:

    res_ch = modifyBins(pred_ch.join(gfa_ch), "naive")

    emit:
    pred = res_ch
}

workflow MODEL {
  
    take:

	// mode
    // scoring
	gfa_ch

	gc_ch
	gd_ch

    type

    main:
    if (type == "pan") {
        model_ch = gfa_ch.join(gc_ch).join(gd_ch)
        model_ch = model_ch.combine( Channel.fromList(ptylist) ).map{ meta, _gfa, _gc, _pl, pty -> meta += ['pty':pty]; [meta, _gfa, _gc, _pl]}

	    bins_ch = model(model_ch, type)
    }
    else {
	bins_ch = model(gfa_ch.join(gc_ch).join(gd_ch), type)
    }
    pred_ch = transform(bins_ch.join(gfa_ch))
   
    naiveBins_ch        = Channel.empty()
    graphOverlapBins_ch = Channel.empty()

    if (type == "pan") {
        NAIVE(gfa_ch, pred_ch)
        OVERLAP(gfa_ch, pred_ch)
        naiveBins_ch        = NAIVE.out.pred
        graphOverlapBins_ch = OVERLAP.out.pred
    }
    

    emit:

    res = pred_ch
    bins = bins_ch
    naive = naiveBins_ch
    overlap = graphOverlapBins_ch
}
