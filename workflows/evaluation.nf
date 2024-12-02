#!/usr/bin/env nextflow

/*
EVALUATION pre-preprocessing
. unicycler and skesa bins
. sample name -> to retrive ground truth
. compute unicycler and skesa ground truth LABELING and BINS
. evaluate unicycler and skesa labeling and binning

. choose MLPLASMID or PLSDENSITY
. panassembly and panassembly-100 bins
. computa unicycler and skesa gt LABELING and BINS in respect of PANASSEMBLY fragments
. compute unicycler and skesa -100 gt LABELING and BINS in respect of PANASSEMBLY-100 fragments
. evaluate pan-assembly and pan-assembly 100 labeling and binning against (unicycler,skesa) and (uniycler,skesa -100)


LABELING:
. evaluate precision recall and f1 of how many fragments/contigs * lengths have been correctly evaluated as plasmid

BINNING:
. evaluate precision recall and f1 using PLASEVAL and as describe in plaseval, using the cumulative-length value

*/

/*

current evaluation
MLPLASMID + pan-assembly/pan-assembly-100 vs Unicycler vs Skesa with super-bins

next:
plasmidness, extended+super bins
*/


process LABELING {

    input:
    tuple val(meta), path(prediction), path(gt), val(tool)

    output:
    tuple val(meta), path(stats), emit: stats
    // tuple val(meta), path(plots), emit: plots

    script:

    output = "${meta.id}.${tool}"
    description = "${meta.id}.${meta.thr} with ${tool}"
    stats = "${output}.stats.txt"

    """
    python $projectDir/bin/evaluate_bins.py --bin ${prediction} --csv ${gt} --sample ${meta.id} --output ${output} --description '${description}'
    """
}

process labeling {

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



process PLASEVAL {
    cache false

    input:
    tuple val(meta), path(prediction), path(gt), val(tool)

    output:
    tuple val(meta), path(stats), emit: stats
    // tuple val(meta), path(plots), emit: plots

    script:

    output = "${meta.id}.${tool}"
    // description = "${meta.id}.${meta.thr} with ${tool}"
    // plots = "${output}.scores.pdf"
    stats = "${output}.plaseval.txt"

    // awk -i inplace '{\$0=gensub(/\s*\S+/,\\"\\",3)}1' ${gt} 
    // awk -i inplace '{\$0=gensub(/\s*\S+/,\\"\\",4)}1' ${gt} 
    """
    #!/bin/bash
    echo 'plasmid\tcontig\tnan\tnan\tcontig_len' | cat - ${gt} > temp && mv temp ${gt}
    python $projectDir/PlasEval/src/plaseval.py eval --pred ${prediction} --gt ${gt} --out ${stats} 
    """
    // python $projectDir/bin/evaluate_bins.py --bin ${prediction} --csv ${gt} --sample ${meta.id} --output ${output} --description '${description}'

}


process binning { //plaseval
    cache false

    input:
    tuple val(meta), path(prediction), path(gt), val(tool)

    output:
    tuple val(meta), path(stats), emit: stats
    // tuple val(meta), path(plots), emit: plots

    script:

    output = "${meta.id}.${tool}"
    // description = "${meta.id}.${meta.thr} with ${tool}"
    // plots = "${output}.scores.pdf"
    stats = "${output}.plaseval.txt"

    // awk -i inplace '{\$0=gensub(/\s*\S+/,\\"\\",3)}1' ${gt} 
    // awk -i inplace '{\$0=gensub(/\s*\S+/,\\"\\",4)}1' ${gt} 
    """
    #!/bin/bash
    echo 'plasmid\tcontig\tnan\tnan\tcontig_len' | cat - ${gt} > temp && mv temp ${gt}
    python $projectDir/PlasEval/src/plaseval.py eval --pred ${prediction} --gt ${gt} --out ${stats} 
    """
    // python $projectDir/bin/evaluate_bins.py --bin ${prediction} --csv ${gt} --sample ${meta.id} --output ${output} --description '${description}'

}


workflow {

    fasta_p_ch = Channel.fromPath("${params.input}/*.fa", type: 'file', checkIfExists: true)
    fasta_s_ch = Channel.fromPath("${params.input}/*S*-s*.fasta", type: 'file', checkIfExists: true)
    fasta_u_ch = Channel.fromPath("${params.input}/*S*-u*.fasta", type: 'file', checkIfExists: true)


    panassembly_gfa = Channel.fromPath("${params.input}/*S*panasm.gfa", type: 'file', checkIfExists: true)
	panassembly_gfa = panassembly_gfa.map{file ->
        def fmeta = [:];
        fmeta.id = file.toString().split('\\.')[1];
        // fmeta.species = file.toString().split('\\.')[0].split('/')[-1];
        fmeta.thr = file.toString().split('\\.')[2];
        [fmeta, file]
    }
    
    fasta_p_ch = fasta_p_ch.map{file ->
        def fmeta = [:];
        fmeta.id = file.toString().split('/')[-1].split('\\.')[1];
        // fmeta.species = file.toString().split('/')[-1].split('\\.')[0];
        fmeta.thr = file.toString().split('/')[-1].split('\\.')[2];
        [fmeta, file]
    }
    fasta_s_ch = fasta_s_ch.map{file ->
        def fmeta = [:];
        fmeta.id = file.toString().split('/')[-1].split('\\-')[1];
        // fmeta.species = file.toString().split('/')[-1].split('\\-')[0];
        fmeta.thr = file.toString().split('/')[-1].split('\\.')[-2];
        [fmeta, file]
    }

    fasta_u_ch = fasta_u_ch.map{file ->
        def fmeta = [:];
        fmeta.id = file.toString().split('/')[-1].split('\\-')[1];
        // fmeta.species = file.toString().split('/')[-1].split('\\-')[0];
        fmeta.thr = file.toString().split('/')[-1].split('\\.')[-2];
        [fmeta, file]
    }

    meta_ch = fasta_p_ch.map{meta, files -> meta}

    NCBI(meta_ch)

    BLAST(panassembly_gfa.join(fasta_p_ch).join(fasta_u_ch).join(fasta_s_ch).join(NCBI.out.ref))

    // PUBLISH(BLAST.out.unigt.mix(BLAST.out.skegt))

    pan_uni = BLAST.out.unigt
    pan_ske = BLAST.out.skegt

    // pan_uni = Channel.fromPath("${params.gt}/*.uni.gt.tsv", type: 'file', checkIfExists: true)
    // pan_ske = Channel.fromPath("${params.gt}/*.ske.gt.tsv", type: 'file', checkIfExists: true)

    // pan_uni = pan_uni.map{ file ->  
    //     def fmeta = [:];
    //     fmeta.id = file.toString().split('/')[-1].split('\\.')[0];
    //     [fmeta, file]
    // }
    // pan_ske = pan_ske.map{ file ->  
    //     def fmeta = [:];
    //     fmeta.id = file.toString().split('/')[-1].split('\\.')[0];
    //     [fmeta, file]
    // }

    // asm_res = Channel.fromFilePairs("${params.asm}/*-{s,u}*.pred.tab", type: 'file', checkIfExists: true)

    // ske_res = asm_res.map{meta, pairs ->
    //     def fmeta = [:];
    //     file = pairs[0].toString().split('/')[-1]
    //     fmeta.id = file.split('\\.')[1].split('-')[0];
    //     fmeta.species = file.split('\\.')[0];
    //     fmeta.thr = file.split('\\.')[2];
    //     [fmeta, pairs[0]]
    // }
    // uni_res = asm_res.map{meta, pairs ->
    //     def fmeta = [:];
    //     file = pairs[0].toString().split('/')[-1]
    //     fmeta.id = file.split('\\.')[1].split('-')[0];
    //     fmeta.species = file.split('\\.')[0];
    //     fmeta.thr = file.split('\\.')[2];
    //     [fmeta, pairs[1]]
    // }

	pbf_res = Channel.fromPath("${params.pbf}/*.pred.tab", type: 'file', checkIfExists: true)
	pbf_res = pbf_res.map{file ->
        def fmeta = [:];
        fmeta.id = file.toString().split('/')[-1].split('\\.')[1];
        // fmeta.species = file.toString().split('/')[-1].split('\\.')[0];
        fmeta.thr = file.toString().split('/')[-1].split('\\.')[2];
        [fmeta, file]
    }
    // pbf_res = MOD_BINS(pbf_res, "--naive --n 1", "A")
    // pbf_res = MOD_BINS(pbf_res, "--naive --n 2", "B")
    pbf_res = MOD_BINS(pbf_res.join(panassembly_gfa), "--super --graph", "C")

    pbf_panu_ch = pbf_res.join(BLAST.out.panunigt.map{id, pan -> [id, pan, "pbf.psm.u"]})
    pbf_pans_ch = pbf_res.join(BLAST.out.panskegt.map{id, pan -> [id, pan, "pbf.psm.s"]})


    // pbf_ske_ch = ske_res.join(BLAST.out.skegt.map{id, ske -> [id, ske, "pbf.ske"]})
    // pbf_uni_ch = uni_res.join(BLAST.out.unigt.map{id, ske -> [id, ske, "pbf.uni"]})

    // evaluate_ch = (pbf_ske_ch).mix(pbf_uni_ch)
    evaluate_ch = (pbf_panu_ch).mix(pbf_pans_ch)

    // EVAL(evaluate_ch)
    PLASEVAL(evaluate_ch)
    PLASEVAL.out.stats | view

    // PUBLISH(EVAL.out.stats.mix(PLASEVAL.out.stats))
    PUBLISH(PLASEVAL.out.stats)
}

