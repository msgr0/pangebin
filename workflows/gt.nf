#!/usr/bin/env nextflow

process ncbi {
    storeDir "${params.input}/"

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
    #!/usr/bin/env bash

    python $projectDir/bin/evaluation/ncbi_link.py --input ${meta.id} --output ${referencegz}
    bgzip -d -c ${referencegz} > ${reference}

    python $projectDir/bin/evaluation/strip_plasmid_fasta.py --input ${reference} --output ${reference_ren}
    """
    

}

process blast {
    storeDir "${params.input}/"

    input:
    tuple val(meta), path(graph), path(mix), path(uni), path(ske)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path(pan_mix_gt), path(mix_gt), emit: pangt
    tuple val(meta), path(pan_uni_gt), path(uni_gt), emit: unigt
    tuple val(meta), path(pan_ske_gt), path(ske_gt), emit: skegt


    script:
    outmix = meta.id + "." + meta.thr + ".mix" 
    // mapping_file = output + ".mapping.tsv"
    mix_gt = outmix + ".gt.tsv"
    pan_mix_gt = outmix + ".pan.gt.tsv"

    outuni = meta.id + "." + meta.thr + ".uni"
    // mapping_file = output + ".mapping.tsv"
    uni_gt = outuni + ".gt.tsv"
    pan_uni_gt = outuni + ".pan.gt.tsv"

    ootske = meta.id + "." + meta.thr + ".ske"
    // mapping_file = output + ".mapping.tsv"
    ske_gt = outske + ".gt.tsv"
    pan_ske_gt = outske + ".pan.gt.tsv"
    
    """
    #!/usr/bin/env bash

    python $projectDir/bin/evaluation/build_truth.py --pangenome ${graph} --assembly ${mix} --reference ${reference} --output ${outmix}
    
    python $projectDir/bin/evaluation/build_truth.py --pangenome ${graph} --assembly ${uni} --reference ${reference} --output ${outuni}

    python $projectDir/bin/evaluation/build_truth.py --pangenome ${graph} --assembly ${ske} --reference ${reference} --output ${outske}
    """

}

workflow GT {
    take:
    
    meta
    input_ch


    main:

    ncbi( meta )
    
    blast( input_ch, ncbi.out.ref )

    emit:
    mixReference = blast.out.pangt
    uniReference = blast.out.unigt
    skeReference = blast.out.skegt
}
