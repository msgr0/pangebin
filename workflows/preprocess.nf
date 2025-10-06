#!/usr/bin/env nextflow

include { PANGENOME                } from "$projectDir/modules/nf-core/pangenome/workflows/pangenome"

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
def species_map = ["efea": "Enterococcus faecium", "kpne": "Klebsiella pneumoniae", "abau": "Acinetobacter baumannii", "ecol": "Escherichia coli", "oth": "other"]
def pctid = [0.95, 0.92, 0.98]
def cutlen = [0, 1, 50, 100]
def thr = [250, 500, 1000, 2000, 5000]

process extractGfa {
    cache 'lenient'

    input:
    tuple val(meta), path(gfa)
    
    output:
    tuple val(meta), path(gfao)

    script:
    if (meta.asm == 's') { // transform skesa graph into unicycler
        param = '-c'
        asm = "ske"
    }
    else if (meta.asm == 'u') {
        param = ''
        asm = "uni" 
    }

    basen = metaname(meta)
    gfao = basen + ".gfa"

    """
    #!/usr/bin/env bash
    bgzip -d -c ${gfa} > ${basen}.tmp.gfa
    python $projectDir/bin/rename_gfa.py -i ${basen}.tmp.gfa -o ${gfao} -p ${asm} ${param}
    """
}

process removeNodes {
    cache 'lenient'
    
    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path(gfat), emit: gfa

    script:
    gfat = metaname(meta) + ".gfa"

    if (meta.cutlen == 0 || meta.cutlen == null)
    """
    #!/usr/bin/env bash
    mv ${gfa} ${gfat}
    """
    else
    """
    #!/usr/bin/env bash
    python $projectDir/bin/remove_nodes.py -v 2 -i ${gfa} -o ${gfat} -t ${meta.cutlen}
    """
}

process extractFasta {
    cache 'lenient'

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val (meta), path(fasta), emit: fasta

    script:
    fasta = metaname(meta) + ".fa"

    """
    #!/usr/bin/env bash
    awk '/^S/{print ">"\$2; print ""\$3}' ${gfa} > ${fasta}    
    """
}

process mixFasta {
    cache 'lenient'

    input:
    tuple val(meta), path(fasta_a), path(fasta_b)
    output:
    tuple val(meta), path(fasta)

    script:
    fasta = metaname(meta) + ".fa"
    """
    cat ${fasta_a} ${fasta_b} > ${fasta}
    """
}

process makePangenome {
    executor 'slurm'
    memory '16 GB'
    cpus 16
    time '12h'
    cache 'lenient'

    input:
    tuple val(meta), path(mixed_fasta)
    
    output:
    tuple val(meta), path(pangenome), emit: pangenome

    script:

    PANGENOME()

    pangenome = metaname(meta) + "pan.gfa"
    out_core = "${meta.id}_nfcore"
    haplos = 2
    paramfile = "${meta.id}.params.json"
    release = params.pan_release
    profile = params.panexec
    """
    #!/usr/bin/env bash
    export NXF_OPTS="-XX:ActiveProcessorCount=16 -Xmx16g"
    module restore pbf
    source $projectDir/.venv/bin/activate
    touch ${paramfile}
    echo -e '{\n\t"max_memory": "16.GB",\n\t"wfmash_segment_length": ${meta.thr},\n\t"seqwish_min_match_length": 0,\n\t"smoothxg_poa_params": "asm5",\n\t"wfmash_map_pct_id": ${meta.pctid},\n\t"max_cpus": ${task.cpus},\n\t"wfmash_merge_segments": true\n}' > ${paramfile}

    bgzip ${mixed_fasta}
    nextflow run nf-core/pangenome -r ${release} -profile $profile -resume --input ${mixed_fasta}.gz --n_haplotypes ${haplos} --outdir ${out_core} -params-file ${paramfile} --wfmash_segment_length ${meta.thr} --wfmash_map_pct_id ${meta.pctid}
    cp ${out_core}/FINAL_GFA/${mixed_fasta}.gz.gfaffix.unchop.Ygs.view.gfa ${pangenome}
    """
}

process makePanassembly {
    executor 'slurm'
    memory '16 GB'
    cpus 8
    time '6h'
    cache 'lenient'

    input:
    tuple val(meta), path(pangenome), path(ske), path(uni)

    output:
    tuple val(meta), path(panassembly), emit: panassembly
    // tuple val(meta), path(panassemblycut), emit: panassembly

    script:
    panassembly = metaname(meta) + "asm.gfa"
    pangenomecl = metaname(meta) + "cl.gfa"
    """
    #!/usr/bin/env bash

    module restore pbf
    source $projectDir/.venv/bin/activate
    python $projectDir/bin/gfa_cleaner.py --input ${pangenome} --output ${pangenomecl}
    python $projectDir/bin/pangenome_enancher.py --input ${pangenomecl} --output ${panassembly} --assembly ${uni} ${ske}
    """
}


process computeScores {
    storeDir "${params.output}/${meta.id}/"
    cache 'lenient'

    input:
    tuple val(meta), path(skesa), path(unicycler)

    output:
    tuple val(meta), path(gc_ske), emit: gc_ske
    tuple val(meta), path(gd_ske), emit: gd_ske
    
    tuple val(meta), path(gc_uni), emit: gc_uni
    tuple val(meta), path(gd_uni), emit: gd_uni
    

    script:
    base_ske = metaname(meta) + ".s"
    base_uni = metaname(meta) + ".u"

    gd_ske = base_ske + ".gd.tsv"
    gd_uni = base_uni + ".gd.tsv"

    gc_ske = base_ske + ".gc.tsv"
    gc_uni = base_uni + ".gc.tsv"

    input_csv = metaname(meta) + ".in.csv"

    putils = "$projectDir/PlasBin-flow-pangenome/code/plasbin_utils.py"
    pls_db_file = "$projectDir/PlasBin-flow-pangenome/database/genes.fasta"

    """
    #!/usr/bin/env bash

    bgzip -k ${skesa}
    bgzip -k ${unicycler}

    echo "sample,gfa" > "${input_csv}"
    echo "${base_ske},${skesa}.gz" >> "${input_csv}"
    echo "${base_uni},${unicycler}.gz" >> "${input_csv}"

    python ${putils} preprocessing --input_file ${input_csv} --out_dir . --tmp_dir ./tmp --out_file out_file.csv --db_file ${pls_db_file}
  
    """
}

process computeAllScores {
    storeDir "${params.input}/"
    cache 'lenient'

    input:
    tuple val(meta), path(pangenome), path(skesa), path(unicycler)

    output:
    tuple val(meta), path(gc_pan), emit: gc_pan
    tuple val(meta), path(gd_pan), emit: gd_pan
    tuple val(meta), path(gc_ske), emit: gc_ske
    tuple val(meta), path(gd_ske), emit: gd_ske
    tuple val(meta), path(gc_uni), emit: gc_uni
    tuple val(meta), path(gd_uni), emit: gd_uni
    

    script:
    base_pan = metaname(meta) + ".p"
    base_ske = metaname(meta) + ".s"
    base_uni = metaname(meta) + ".u"
    base_mix = metaname(meta) + ".m"

    gd_pan = base_pan + ".gd.tsv"
    gd_ske = base_ske + ".gd.tsv"
    gd_uni = base_uni + ".gd.tsv"

    gc_pan = base_pan + ".gc.tsv"
    gc_ske = base_ske + ".gc.tsv"
    gc_uni = base_uni + ".gc.tsv"

    input_csv = metaname(meta) + ".in.csv"

    putils = "$projectDir/PlasBin-flow-pangenome/code/plasbin_utils.py"
    pls_db_file = "$projectDir/PlasBin-flow-pangenome/database/genes.fasta"

    """
    #!/usr/bin/env bash

    bgzip -k ${skesa}
    bgzip -k ${unicycler}

    echo "sample,gfa" > "${input_csv}"
    echo "${base_ske},${skesa}.gz" >> "${input_csv}"
    echo "${base_uni},${unicycler}.gz" >> "${input_csv}"

    python ${putils} preprocessing --input_file ${input_csv} --out_dir . --tmp_dir ./tmp --out_file out_file.csv --db_file ${pls_db_file}



    cat ${gd_ske} ${gd_uni} > ${base_mix}.gd.tsv
    cat ${gc_ske} ${gc_uni} > ${base_mix}.gc.tsv

    python3 $projectDir/bin/pbf_preprocess.py --input ${pangenome} --gd  ${base_mix}.gd.tsv --gc ${base_mix}.gc.tsv --output ${base_pan}
  
    """
}
workflow PANGENOME_WRAP {
    take:
    mixFasta_ch

    main:

    meta = mixFasta_ch.map { meta, _files -> meta }
    input_ch = mixFasta_ch.map { meta, files -> [[meta.id], files] }
    PANGENOME(input_ch)

    out_ch = Channel.empty()
    out_ch = meta.combine(PANGENOME.out.gfa)

    emit:
    pangenome = out_ch
}

workflow PREPROCESS {
    take:
    input_ch // short reads assembled by unicycler and skesa in 2 .gfa.gz

    main:
    extracted_gfa = extractGfa(input_ch)

    extracted_gfa.combine( Channel.fromList(cutlen) ).map{meta, files, _cutlen -> meta += [cutlen: _cutlen]; [meta, files]}
    | removeNodes
    | extractFasta

    // extracted_gfa.map{meta, files -> meta += [cutlen:params.cutlen]; [meta, files]}
    // | removeNodes
    // | extractFasta

    skesaGfa_ch = removeNodes.out.gfa.filter { meta, file -> meta['asm'] == 's' }

    uniGfa_ch = removeNodes.out.gfa.filter { meta, file -> meta['asm'] == 'u' }

    skesaFasta_ch = extractFasta.out.fasta.filter { meta, file -> meta['asm'] == 's' }

    uniFasta_ch = extractFasta.out.fasta.filter { meta, file -> meta['asm'] == 'u' }
    mixFasta_ch = Channel.empty() 
    pangeGfa_ch = Channel.empty()
    panasmGfa_ch = Channel.empty()

    gcPan_ch = Channel.empty()
    gdPan_ch = Channel.empty()
    mlPan_ch = Channel.empty()
    
    gcSke_ch = Channel.empty()
    gdSke_ch = Channel.empty()
    mlSke_ch = Channel.empty()

    gcUni_ch = Channel.empty()
    gdUni_ch = Channel.empty()
    mlUni_ch = Channel.empty()

    result_ch = Channel.empty()


    if (params.pangenome) {
        mixFasta_ch = skesaFasta_ch.map{meta, files -> meta_t = [:]; meta_t['id'] = meta['id']; meta_t['cutlen'] = meta['cutlen']; [meta_t, files]
            }.combine(uniFasta_ch.map{meta, files -> meta_t = [:]; meta_t['id'] = meta['id']; meta_t['cutlen'] = meta['cutlen']; [meta_t, files]}, by: 0).view()
        | mixFasta
        mixFasta_ch = mixFasta_ch.combine( Channel.fromList(pctid) ).combine( Channel.fromList(thr) ).map{meta, files, _pctid, _thr -> meta += [pctid: _pctid, thr: _thr]; [meta, files]}


        
        // pangenome_ch = mixFasta_ch | makePangenome
        
        PANGENOME_WRAP(mixFasta_ch)
        pangenome_ch = PANGENOME_WRAP.out.pangenome

        skesa_ch = skesaGfa_ch.combine( Channel.fromList(pctid) ).combine( Channel.fromList(thr) ).map{meta, files, _pctid, _thr -> meta += [pctid: _pctid, thr: _thr]; [meta, files]}
        uni_ch = uniGfa_ch.combine( Channel.fromList(pctid) ).combine( Channel.fromList(thr) ).map{meta, files, _pctid, _thr -> meta += [pctid: _pctid, thr: _thr]; [meta, files]}
        result = Channel.empty()
        panasmGfa_ch = makePanassembly(pangenome_ch.join(skesa_ch).join(uni_ch)) 
        result = computeAllScores(panasmGfa_ch.join(skesa_ch).join(uni_ch))
        gcPan_ch = gcPan_ch.mix(result.gc_pan)
        gdPan_ch = gdPan_ch.mix(result.gd_pan)

    }   
    if (params.assemblers && ! params.pangenome) {         
        result = computeScores(skesaGfa_ch.join(uniGfa_ch))
        // if (params.ml) {
        //     mlres = computeMLscores(skesaGfa_ch.join(uniGfa_ch))
        //     mlSke_ch = mlSke_ch.mix(mlres.ske)
        //     mlUni_ch = mlUni_ch.mix(mlres.uni)
        // }
    }

    gcSke_ch = result.gc_ske
    gdSke_ch = result.gd_ske
    gcUni_ch = result.gc_uni
    gdUni_ch = result.gd_uni

    
    emit:

    mixFasta = mixFasta_ch
    skeFasta = skesaFasta_ch
    uniFasta = uniFasta_ch

    pangeGfa = pangeGfa_ch
    panasmGfa = panasmGfa_ch
    skesaGfa = skesaGfa_ch
    uniGfa = uniGfa_ch

    gcPan = gcPan_ch
    gdPan = gdPan_ch
    mlPan = mlPan_ch

    gcSke = gcSke_ch
    gdSke = gdSke_ch
    mlSke = mlSke_ch

    gcUni = gcUni_ch
    gdUni = gdUni_ch
    mlUni = mlUni_ch

}
