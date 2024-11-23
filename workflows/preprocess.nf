#!/usr/bin/env nextflow
/*
input:   - a pair of assembly graphs
output:  - cleaned assembly graphs, standardized on Unicycler
          (eventually cutted in 0/100)
          - pan-assembly graph
*/
include { PUBLISH } from "./utils"


paramfile = "$projectDir/pangenome-params.json"


process EXTRACT_GFA {
    cache 'lenient'

    input:
    tuple val(meta), path(gfa)
    
    output:
    tuple val(meta), path(gfa_std)

    script:
    if (meta.asm == 's') {
        param = '-c'
        asm = "ske"
    }
    else if (meta.asm == 'u') {
        param = ''
        asm = "uni" 
    }
    
    gfa_xtracted = gfa.getBaseName()
    gfa_std = gfa.getBaseName()[0..-5] + ".r.gfa"

    """
    bgzip -d ${gfa}
    python $projectDir/bin/rename_gfa.py -i ${gfa_xtracted} -o ${gfa_std} -p ${asm} ${param}
    """
}

process REMOVE_NODES {
    cache 'lenient'
    
    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path(trimmed), emit: gfa

    script:
    trimmed = gfa.getBaseName() + ".${meta.thr}.gfa"
    if (meta.thr == 0 )
    """
    mv ${gfa} ${trimmed}
    """
    else
    """
    python $projectDir/bin/remove_nodes.py -v 2 -i ${gfa} -o ${trimmed} -t ${meta.thr}
    """
}


process EXTRACT_FASTA {
    cache 'lenient'

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val (meta), path(fasta), emit: fasta

    script:
    fasta = gfa.getBaseName() + ".fasta"

    """
    awk '/^S/{print ">"\$2; print ""\$3}' ${gfa} > ${fasta}    
    """
}

process MIX_FASTA {
    cache 'lenient'

    input:
    tuple val(meta), path(fasta_a), path(fasta_b)
    output:
    tuple val(meta), path(fasta)

    script:
    fasta = "${meta.species}.${meta.id}.${meta.thr}.fa"
    """
    cat ${fasta_a} ${fasta_b} > ${fasta}
    """
}


process MAKE_PANGENOME {
    maxForks 4
    time '60m'
    cache 'lenient'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(mixed_fasta)
    
    output:
    tuple val(meta), path(pangenome)

    script:
    pangenome = "${meta.id}.pan.gfa"
    out_core = "${meta.id}_nfcore"
    haplos = 2
    paramfile = "$projectDir/pangenome-params.json"
    release = '1.1.2'
    profile = 'docker'

    """
    bgzip ${mixed_fasta}
    nextflow run nf-core/pangenome -r ${release} -profile $profile -resume --input ${mixed_fasta}.gz --n_haplotypes ${haplos} --outdir ${out_core} -params-file ${paramfile}
    cp ${out_core}/FINAL_GFA/${mixed_fasta}.gz.gfaffix.unchop.Ygs.view.gfa ${pangenome}
    """

}

process MAKE_PANASSEMBLY {
    cache 'lenient'
        
    input:
    tuple val(meta), path(pangenome), path(ske), path(uni)

    output:
    tuple val(meta), path(panassembly), emit: panassembly
    // tuple val(meta), path(panassemblycut), emit: panassembly

    script:
    panassembly = pangenome.baseName + "asm.gfa"
    pangenomecl = pangenome.baseName + "cl.gfa"
    """
    python $projectDir/bin/gfa_cleaner.py --input ${pangenome} --output ${pangenomecl}
    python $projectDir/bin/pangenome_enancher.py --input ${pangenomecl} --output ${panassembly} --assembly ${uni} ${ske}
    """
}


workflow PREPROCESS {
    take:
    gfagz

    main:
    extract_ch = EXTRACT_GFA(gfagz)
    extract_ch.map{meta, files -> meta += [thr:params.cutlen]; [meta, files]} | REMOVE_NODES | EXTRACT_FASTA

    skesa_gfa = REMOVE_NODES.out.gfa.filter { meta, file -> meta['asm'] == 's' }.map{
        meta, files -> meta = [id: meta.id, species: meta.species, thr:meta.thr]; [meta, files]
    }

    unicycler_gfa = REMOVE_NODES.out.gfa.filter { meta, file -> meta['asm'] == 'u' }.map{
        meta, files -> meta = [id: meta.id, species: meta.species, thr:meta.thr]; [meta, files]
    }

    skesa_fasta = EXTRACT_FASTA.out.fasta.filter { meta, file -> meta['asm'] == 's' }.map{
        meta, files -> meta = [id: meta.id, species: meta.species, thr:meta.thr]; [meta, files]
    } | view

    unicycler_fasta = EXTRACT_FASTA.out.fasta.filter { meta, file -> meta['asm'] == 'u' }.map{
        meta, files -> meta = [id: meta.id, species: meta.species, thr:meta.thr]; [meta, files]
    } | view

    mixed_fasta = skesa_fasta.combine(unicycler_fasta, by: 0) | MIX_FASTA | view

    pangenome_gfa = mixed_fasta | MAKE_PANGENOME

    panassembly_gfa = MAKE_PANASSEMBLY(pangenome_gfa.join(skesa_gfa).join(unicycler_gfa)) 

    emit:
    fasta_mix = mixed_fasta
    fasta_ske = skesa_fasta
    fasta_uni = unicycler_fasta

    gfa_pan = pangenome_gfa
    gfa_psm = panassembly_gfa
    gfa_ske = skesa_gfa
    gfa_uni = unicycler_gfa

}
