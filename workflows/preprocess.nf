// input:   - a pair of assembly graphs
// output:  - cleaned assembly graphs, standardized on Unicycler
//          (eventually cutted in 0/100)
//          - pan-assembly graph

include { PUBLISH } from "./utils"

params.thr = 0

nfcore_profile = "apptainer"
paramfile = "$projectDir/pangenome-params.json"
release = '1.1.2'

seed_len = 1000
seed_score = 0.5
min_pls_len = 1000

binfolder = "~/bin"

putils = "$projectDir/PlasBin-flow-pangenome/code/plasbin_utils.py"
pflow = "$projectDir/PlasBin-flow-pangenome/code/plasbin_flow.py"
pdatabase = "$projectDir/Plasbin-flow-pangenome/database/genes.fasta"

process EXTRACT_GFA {
    cache 'lenient'

    input:
    tuple val(meta), path(gfa)
    output:
    tuple val(meta), path(out_gfa)
    script:
    dec_gfa = gfa.baseName
    out_gfa = gfa.baseName

    if (meta.asm == 'ske') {
        param = '-c' // convert
    }
    else {
        param = '' 
    }
    """
    bgzip -d ${gfa}
    python $projectDir/bin/rename_gfa.py -i ${dec_gfa} -o ${out_gfa} -p ${meta.asm} ${param}
    """
}

process REMOVE_NODES {
    cache 'lenient'
    
    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path(trimmed), emit: gfa

    script:
    trimmed = gfa.baseName + ".${meta.thr}.gfa"
    """
    if [ ${meta.thr} == 0]; then
        mv ${gfa} ${trimmed}
    else
        python $projectDir/bin/remove_nodes.py -v 2 -i ${gfa} -o ${trimmed} -t ${meta.thr}
    fi
    """
}


process EXTRACT_FASTA {
    cache 'lenient'

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val (meta), path(fasta), emit: fasta

    script:
    fasta = gfa.baseName + ".fasta"
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
    fasta = "${meta.species}.${meta.id}.mix.fa"
    """
    cat ${fasta_a} ${fasta_b} > ${fasta}
    """
}

process MAKE_PANGENOME {
    cache 'lenient'

    maxForks 4
    time '45m'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(pangenome), emit: pangenome

    script:
    haplos = 2
    pangenome = "${meta.species}.${meta.id}.pan.gfa"

    """
    bgzip ${fasta}
    nextflow run nf-core/pangenome -r '${release}' -profile ${nfcore_profile} --input ${fasta}.gz --n_haplotypes ${haplos} --outdir . -params-file ${paramfile} -w ${task.workDir} 
    cp ./FINAL_GFA/${fasta}.gz.gfaffix.unchop.Ygs.view.gfa ${pangenome}
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
    panassembly = "${meta.species}.${meta.id}.pasm.gfa"
    pangenomecl = "${meta.species}.${meta.id}.pan.cl.gfa"
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
    extract_ch.map{meta, files -> meta += [thr:params.thr]; [meta, files]} | REMOVE_NODES | EXTRACT_FASTA

    skesa_gfa = REMOVE_NODES.out.gfa.filter { meta, file -> meta['asm'] == 'ske' }.map{
        meta, files -> meta = [id: meta.id, species: meta.species, thr:meta.thr]; [meta, files]
    }

    unicycler_gfa = REMOVE_NODES.out.gfa.filter { meta, file -> meta['asm'] == 'uni' }.map{
        meta, files -> meta = [id: meta.id, species: meta.species, thr:meta.thr]; [meta, files]
    }

    skesa_fasta = EXTRACT_FASTA.out.fasta.filter { meta, file -> meta['asm'] == 'ske' }.map{
        meta, files -> meta = [id: meta.id, species: meta.species, thr:meta.thr]; [meta, files]
    }

    unicycler_fasta = EXTRACT_FASTA.out.fasta.filter { meta, file -> meta['asm'] == 'uni' }.map{
        meta, files -> meta = [id: meta.id, species: meta.species, thr:meta.thr]; [meta, files]
    }

    mixed_fasta = skesa_fasta.combine(skesa_fasta, by: 0) | MIX_FASTA

    pangenome_gfa = mixed_fasta | MAKE_PANGENOME

    panassembly_gfa = MAKE_PANASSEMBLY(pangenome_gfa.join(skesa_gfa).join(unicycler_gfa)) 

    emit:
    pangenome = pangenome_gfa
    panassembly = panassembly_gfa
    pan_fasta = mixed_fasta
    ske_fasta = skesa_fasta
    uni_fasta = unicycler_fasta
    ske_gfa = skesa_gfa
    uni_gfa = unicycler_gfa
    all = pangenome_gfa.mix(panassembly_gfa).mix(mixed_fasta).mix(skesa_fasta).mix(unicycler_fasta).mix(skesa_gfa).mix(unicycler_gfa)
}

workflow {
    input_ch = Channel.fromFilePairs("${params.input}/*-S*-{s,u}.gfa.gz", type: 'file', checkIfExists: true)

    input_ch = input_ch.map{ meta, filepair -> def fmeta = [:]; fmeta.id = meta[5..-1]; fmeta.species = meta[0..3]
        [fmeta, filepair]
    }
    meta_ch = input_ch.map{meta, files -> meta}

    unicycler_ch = input_ch | map { meta, files -> meta += [asm: "uni"]; uni = files.findAll { it.toString().contains("-u.gfa.gz")}; [meta, uni[0]]}
    skesa_ch = input_ch | map { meta, files -> meta += [asm: "ske"]; ske = files.findAll { it.toString().contains("-s.gfa.gz")}; [meta, ske[0]]}

    PREPROCESS(unicycler_ch.mix(skesa_ch))
    
    PUBLISH(PREPROCESS.out.all)

}