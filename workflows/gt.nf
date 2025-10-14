#!/usr/bin/env nextflow

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
process ncbi {

    errorStrategy { sleep(Math.pow(2, task.attempt) * 5 as long); return 'retry' }
    maxRetries 5

    input:
    val (meta)

    output:
    tuple val(meta), path(reference_ren), emit: ref

    script:
    name = metaname(meta)
    referencegz = "${name}.fna.gz"
    reference = "${name}.fna"
    reference_ren = "${name}.ren.fna"

    """
    #!/usr/bin/env bash

    python $projectDir/bin/evaluation/ncbi_link.py --input ${meta.id} --output ${referencegz}
    bgzip -d -c ${referencegz} > ${reference}
    if [ -f ${reference} ] && [ ! -s ${reference} ]; then
        echo "File exists and it's empty"
        exit 1
    fi

    python $projectDir/bin/evaluation/strip_plasmid_fasta.py --input ${reference} --output ${reference_ren}
    """
    

}

process blast {
    errorStrategy 'terminate'
    cache false

    input:
    tuple val(meta), path(graph), path(mix), path(ske), path(uni)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path(pan_mix_gt), path(mix_gt), emit: pangt
    tuple val(meta), path(pan_ske_gt), path(ske_gt), emit: skegt
    tuple val(meta), path(pan_uni_gt), path(uni_gt), emit: unigt


    script:
    name = metaname(meta)
    outmix = name + ".mix" 
    // mapping_file = output + ".mapping.tsv"
    mix_gt = outmix + ".gt.tsv"
    pan_mix_gt = outmix + ".pan.gt.tsv"

    outuni = metaname(meta) + ".uni"
    // mapping_file = output + ".mapping.tsv"
    uni_gt = outuni + ".gt.tsv"
    pan_uni_gt = outuni + ".pan.gt.tsv"

    outske = metaname(meta) + ".ske"
    // mapping_file = output + ".mapping.tsv"
    ske_gt = outske + ".gt.tsv"
    pan_ske_gt = outske + ".pan.gt.tsv"
    
    """
    #!/usr/bin/env bash

    python $projectDir/bin/evaluation/build_truth.py --pangenome ${graph} --assembly ${mix} --reference ${reference} --output ${outmix} --mix
    
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
    skeReference = blast.out.skegt
    uniReference = blast.out.unigt
}
