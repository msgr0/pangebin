// include { NFCORE_PANGENOME } from "$projectDir/pangenome"
paramfile = "$projectDir/workflows/params/pangenome-params.json"
release = '1.1.2'
profile = 'conda'
src_dir = "$projectDir/bin"



process MAKE {
    maxForks 2
    cache 'lenient'

    input:
    tuple val(meta), path(fasta_AB), path(fasta_AB_gz)
    
    output:
    tuple val(meta), path(pangenome)

    script:
    pangenome = "${meta.id}.pangenome.gfa"
    out_core = "${meta.id}_nfcore"
    haplos = 2
    """
    nextflow run nf-core/pangenome -r ${release} -profile $profile -resume --input ${fasta_AB_gz} --n_haplotypes ${haplos} --outdir ${out_core} -params-file ${paramfile}
    cp ${out_core}/FINAL_GFA/${fasta_AB_gz}.gfaffix.unchop.Ygs.view.gfa ${pangenome}
    """
}


process CLEAN {
    conda 'gfapy'
    cache 'lenient'

    input:
    tuple val(meta), path(pangenome)
    output:
    tuple val(meta), path(clean_pangenome)

    script:
    clean_pangenome = "${meta.id}.cleaned.pangenome.gfa"
    """
    python $src_dir/gfa_cleaner.py --input ${pangenome} --output ${clean_pangenome}
    """

}

process ENHANCE {
    cache 'lenient'
    conda 'gfapy pandas numpy'
    input:
    tuple val(meta), path(clean_pangenome), path(ass_A), path(ass_B)

    output:
    tuple val(meta), path(augmented_pangenome)

    script:
    augmented_pangenome = "${meta.id}.augmented.pangenome.gfa"
    """
    python $src_dir/pangenome_enancher.py --input ${clean_pangenome} --output ${augmented_pangenome} --assembly ${ass_A} ${ass_B}
    """
}

workflow PANGENOME {
    take:
    assemblies
    pangenome

    main:
    pange_out = MAKE(pangenome)
    // INPUT_CHECK(pangenome_fa)


    // PGGB (
    //     INPUT_CHECK.out.fasta,
    //     INPUT_CHECK.out.fai,
    //     INPUT_CHECK.out.gzi
    // )
    cleaned = CLEAN(pange_out)
    cleaned | view

    augmented = ENHANCE(cleaned.join(assemblies))
    augmented | view
    emit:
    augmented
}



