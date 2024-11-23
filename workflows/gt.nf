#!/usr/bin/env nextflow

// include {PUBLISH; PUBLISH as REF} from "./utils.nf"

process NCBI {
    publishDir "${params.input}/", mode: 'copy'
    // errorStrategy 'ignore'

    input:
    val (meta)

    output:
    tuple val(meta), path(reference_ren), emit: ref

    script:
    name = "${meta.species}-${meta.id}"
    referencegz = "${name}.fna.gz"
    reference = "${name}.fna"
    reference_ren = "${name}.ren.fna"

    """
    python $projectDir/bin/evaluation/ncbi_link.py --input ${meta.id} --output ${referencegz}
    bgzip -d -c ${referencegz} > ${reference}

    python $projectDir/bin/evaluation/strip_plasmid_fasta.py --input ${reference} --output ${reference_ren}
    """
    

}

process BLAST {

    input:
    tuple val(meta), path(graph), path(mix), path(uni), path(ske)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path(pan_mix_gt), path(mix_gt), emit: pangt
    tuple val(meta), path(pan_uni_gt), path(uni_gt), emit: unigt
    tuple val(meta), path(pan_ske_gt), path(ske_gt), emit: skegt


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


workflow DOWNLOAD_GT {
    take:
    meta_ch

    main:
    NCBI(meta_ch)

    emit:
    reference = NCBI.out.ref
}

workflow GT {
    take:
    
    meta
    input_ch


    main:

    NCBI( meta )
    
    BLAST( input_ch, NCBI.out.ref )

    emit:
    all_contigs = BLAST.out.pangt
    uni_contigs = BLAST.out.unigt
    ske_contigs = BLAST.out.skegt
}
