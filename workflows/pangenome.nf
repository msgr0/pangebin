include { NFCORE_PANGENOME } from "$projectDir/pangenome"
paramfile = "$projectDir/workflows/params/pangenome-params.yaml"
profile = docker

process make_pangenome {
    maxForks 2

    input:
    tuple val(id), path(fasta_AB_gz)
    
    output:
    tuple val(id), path(pangenome)

    script:
    pangenome = "${id}.pangenome.gfa"
    out_core = "${id}_nfcore"
    haplos = 2
    """
    nextflow run nf-core/pangenome -r 1.1.2 -profile $profile --input ${fasta_AB_gz} --n_haplotypes ${haplos} --outdir ${out_core} -params-file ${paramfile}
    cp ${out_core}/FINAL_GFA/${fasta_AB_gz}.gfaffix.unchop.Ygs.view.gfa ${pangenome}
    """
}


process clean_pangenome {

    input:
    tuple val(id), path(pangenome)
    output:
    tuple val(id), path(clean_pangenome)

    script:
    clean_pangenome = "${id}.cleaned.pangenome.gfa"
    """
    python $src_dir/gfa_cleaner.py --input ${pangenome} --output ${clean_pangenome}
    """

}

process enhance_pangenome {
    input:
    tuple val(id), path(clean_pangenome)
    path(ass_A)
    path(ass_B)

    output:
    tuple val(id), path(augmented_pangenome)

    script:
    augmented_pangenome = "${id}.augmented.pangenome.gfa"
    """
    python $src_dir/pangenome_enancher.py --input ${clean_pangenome} --output ${augmented_pangenome} --assembly ${ass_A} ${ass_B}
    """
}

workflow PANGENOME {
    take:
    id
    fasta_AB_gz
    ass_A
    ass_B


    main:
    cleaned = make_pangenome(id, fasta_AB_gz)
    | clean_pangenome()

    augmented = enhance_pangenome(cleaned, ass_A, ass_B)

    emit:
    augmented

}




process make_pangenome {
    maxForks 1

    input:
    tuple val(id), path(fasta_AB_gz)

    output:
    tuple val(id), path(pangenome)

    script:


    pangenome = "${id}.pangenome.gfa"
    haplos = 2
    paramfile = "${src_dir}/param-pangenome.yaml"
    """
    nextflow run nf-core/pangenome -r 1.1.1 -profile conda --input ${fasta_AB_gz} --n_haplotypes ${haplos} --outdir ${id}-nfcore -params-file ${paramfile}
    mv ${id}-nfcore/FINAL_GFA/${fasta_AB_gz}.gfaffix.unchop.Ygs.view.gfa ${id}.pangenome.gfa
    """
    
}