include {NFCORE_PANGENOME} from "$projectDir/pangenome"

workflow PANGENOME {
    take:


    main:

    params {
        max_memory = "32GB"
        wfmash_segment_length = "300"
        seqwish_min_match_length = ''
        smoothxg_poa_params = 'asm5'
        wfmash_chunks = 1
        wfmash_map_pct_id = 95
        max_cpus = 16
        input = 
        n_haplotypes = 2
        outdir = "pangenome-out"
    }

    NFCORE_PANGENOME()

    emit:


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

