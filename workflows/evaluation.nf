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


process NCBI {
    // publishDir "${params.input}/", mode: 'copy'
    errorStrategy 'ignore'
    maxForks 1
    input:
    val (meta)

    output:
    tuple val(meta), path(reference_ren), emit: ref

    script:
    name = "${meta.id}"
    referencegz = "${name}.fna.gz"
    reference = "${name}.fna"
    reference_ren = "${name}.ren.fna"

    """
    sleep \$((RANDOM % 5))
    python $projectDir/bin/evaluation/ncbi_link.py --input ${meta.id} --output ${referencegz}
    bgzip -d -c ${referencegz} > ${reference}

    python $projectDir/bin/evaluation/strip_plasmid_fasta.py --input ${reference} --output ${reference_ren}
    """

}



process EVAL {

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

process MOD_BINS {
    cache false

    input:
    tuple val(meta), path(pred), path(graph)
    val(flag)
    val(name)

    output:
    tuple val(meta), path(pred_mod)

    script:
    pred_mod = pred.baseName + ".${name}.tsv"

    """
    python $projectDir/bin/extend_bins.py --pred ${pred} --out ${pred_mod} ${flag} ${graph}
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

process BLAST {

    input:
    tuple val(meta), path(graph), path(mix), path(uni), path(ske), path(reference)

    output:
    tuple val(meta), path(pan_mix_gt), emit: panmixgt
    tuple val(meta), path(mix_gt), emit: mixgt

    tuple val(meta), path(pan_uni_gt), emit: panunigt
    tuple val(meta), path(uni_gt), emit: unigt

    tuple val(meta), path(pan_ske_gt), emit: panskegt
    tuple val(meta), path(ske_gt), emit: skegt


    script:
    outmix = meta.id + "." + "mix"
    // mapping_file = output + ".mapping.tsv"
    mix_gt = outmix + ".gt.tsv"
    pan_mix_gt = outmix + ".pan.gt.tsv"

    outuni = meta.id + "." + "uni"
    // mapping_file = output + ".mapping.tsv"
    uni_gt = outuni + ".gt.tsv"
    pan_uni_gt = outuni + ".pan.gt.tsv"

    outske = meta.id + "." + "ske"
    // mapping_file = output + ".mapping.tsv"
    ske_gt = outske + ".gt.tsv"
    pan_ske_gt = outske + ".pan.gt.tsv"
    
    """
    python $projectDir/bin/evaluation/build_truth.py --pangenome ${graph} --assembly ${mix} --reference ${reference} --output ${outmix}
    
    python $projectDir/bin/evaluation/build_truth.py --pangenome ${graph} --assembly ${uni} --reference ${reference} --output ${outuni}

    python $projectDir/bin/evaluation/build_truth.py --pangenome ${graph} --assembly ${ske} --reference ${reference} --output ${outske}
    """

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

process MOD_BINS {
    cache false

    input:
    tuple val(meta), path(pred), path(graph)
    val(flag)
    val(name)

    output:
    tuple val(meta), path(pred_mod)

    script:
    pred_mod = pred.baseName + ".${name}.tsv"

    """
    python $projectDir/bin/extend_bins.py --pred ${pred} --out ${pred_mod} ${flag} ${graph}
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
// workflow EVAL {
//     take:
//     mode
//     result_ch

//     main:
    
// }