def metaname(meta) {
    def name = meta.id
    if (meta.containsKey('spc')) {
        name += "-${meta.spc}"
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

def cutl = [0, 1, 50, 100]
def species_map = ["efea": "Enterococcus faecium", "kpne": "Klebsiella pneumoniae", "abau": "Acinetobacter baumannii", "ecol": "Escherichia coli", "oth": "other"]
def pctid = [95, 92, 98]
def thr = [250, 500, 1000, 2000, 5000]
// Gfa to Panassembly
process bgzip_d {
    input:
    tuple val(meta), path(gz)

    output:
    tuple val(meta), path(out)

    script:
    out = gz.baseName 

    """
    #!/usr/bin/env bash
    bgzip -d ${gz}
    """
}

process bgzip  {
    input:
    tuple val(meta), path(in)

    output:
    tuple val(meta), path(gz)

    script:
    gz = in + ".gz" 

    """
    #!/usr/bin/env bash
    bgzip -k ${in}
    """
}

process remove_nodes {
    cache 'lenient'
    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path(gfat)


    script:
    gfat = metaname(meta) + ".gfa"

    if (meta.cut == 0 || meta.cut == null)
    """
    #!/usr/bin/env bash
    cp ${gfa} ${gfat}
    """
    else
    """
    #!/usr/bin/env bash
    python $projectDir/bin/remove_nodes.py -v 2 -i ${gfa} -o ${gfat} -t ${meta.cut}
    """

}
//

process rename_gfa {
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

    gfao = metaname(meta) + ".gfa"
    """
    python $projectDir/bin/rename_gfa.py -i ${gfa} -o ${gfao} -p ${asm} ${param}
    """



}

process extract_fasta {
    cache 'lenient'

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val (meta), path(fasta)

    script:
    fasta = metaname(meta) + ".fa"

    """
    #!/usr/bin/env bash
    awk '/^S/{print ">"\$2; print ""\$3}' ${gfa} > ${fasta}    
    """

}

process mix_fasta {
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

// GROUND TRUTH
process refcheck{
    errorStrategy 'ignore'
    input:

    tuple val(meta), path(refin)

    output:
    tuple val(meta), path(ref)

    script:
    ref = metaname(meta) + ".fna"
    """
    python $projectDir/bin/evaluation/strip_plasmid_fasta.py --input ${refin} --output ${ref}
    
    """

}

process gt {

    maxForks 6
    errorStrategy { sleep(Math.pow(2, task.attempt) * 5 as long); return 'retry' }
    maxRetries 2

    input:
    val (meta)

    output:
    tuple val(meta), path(reference)

    script:
    name = metaname(meta)
    referencegz = "${name}.fna.gz"
    reference = "${name}.tmp.fna"

    """
    #!/usr/bin/env bash

    python $projectDir/bin/evaluation/ncbi_link.py --input ${meta.id} --output ${referencegz}
    bgzip -d -c ${referencegz} > ${reference}
    if [ -f ${reference} ] && [ ! -s ${reference} ]; then
        echo "File exists and is empty" # we need to retry - maybe a neterror
        exit 127
    fi

    """

}
//................//
process make_pangenome {
    memory '16 GB'
    executor 'slurm'
    cpus 8 
    time '12h'
    cache 'lenient'

    input:
    tuple val(meta), path(mixed_fasta)
    
    output:
    tuple val(meta), path(pangenome)

    script:
    pangenome = metaname(meta) + ".pan.gfa"
    out_core = metaname(meta) + "_nfcore"
    haplos = 2
    paramfile = metaname(meta) + ".params.json"
    release = params.pan_release
    profile = params.panexec
    """
    #!/usr/bin/env bash
    
    touch ${paramfile}
    echo -e '{\n\t"max_memory": "16.GB",\n\t"wfmash_segment_length": ${meta.thr},\n\t"seqwish_min_match_length": 0,\n\t"smoothxg_poa_params": "asm5",\n\t"wfmash_map_pct_id": ${meta.pctid},\n\t"max_cpus": ${task.cpus},\n\t"wfmash_merge_segments": true\n}' > ${paramfile}

    bgzip ${mixed_fasta}
    nextflow run nf-core/pangenome -r ${release} -profile $profile -resume --input ${mixed_fasta}.gz --n_haplotypes ${haplos} --outdir ${out_core} -params-file ${paramfile} --wfmash_segment_length ${meta.thr} --wfmash_map_pct_id ${meta.pctid}
    cp ${out_core}/FINAL_GFA/${mixed_fasta}.gz.gfaffix.unchop.Ygs.view.gfa ${pangenome}
    """

    stub:
    pangenome = metaname(meta) + ".pan.gfa"
    """
    touch ${pangenome}
    """
}

process make_panassembly {
    executor 'slurm'

    memory '16 GB'
    cpus 1
    time '15m'
    cache 'lenient'

    input:
    tuple val(meta), path(pangenome), path(ske), path(uni)

    output:
    tuple val(meta), path(panassembly)
    // tuple val(meta), path(panassemblycut), emit: panassembly

    script:
    panassembly = metaname(meta) + ".pasm.gfa"
    pangenomecl = metaname(meta) + ".cl.gfa"
    """
    #!/usr/bin/env bash

    # module restore pbf
    # source $projectDir/.venv/bin/activate
    python $projectDir/bin/gfa_cleaner.py --input ${pangenome} --output ${pangenomecl}
    python $projectDir/bin/pangenome_enancher.py --input ${pangenomecl} --output ${panassembly} --assembly ${uni} ${ske}
    """
    stub:
    panassembly = metaname(meta) + ".pasm.gfa"
    """
    touch ${panassembly}
    """
}
process compute_scores {
    cache 'lenient'
    errorStrategy 'terminate'

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
    stub:
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

"""
touch ${gd_pan}
touch ${gd_ske}
touch ${gd_uni}
touch ${gc_pan}
touch ${gc_ske}
touch ${gc_uni}
"""

}

process model {
    executor 'slurm'
    cpus 8
    memory '16 GB'    
    cache 'lenient'

    input:
    tuple val(meta), path(gfa), path(gc), path(plscore)

    
    output:
    tuple val(meta), path(bins)

    script:

    if (meta.t == "pan") {
        assembler = "pangenome"
        // asm = "pan"
        if (meta.pty == 0) {
            args = "--nopenalty"
        }
        else if (meta.pty == 1) {
            args = ""
        }
        seed_len = params.pangenome_seedlen
        seed_score = params.pangenome_seedscore
        min_pls_len = params.pangenome_minplaslen
    }

    else if (meta.t == "asm") {
        assembler = meta.asm == "s" ? "skesa" : "unicycler"
        // asm = meta.asm == "s" ? "ske" : "uni"
        args = "--nopenalty"
        seed_len = params.assembly_seedlen
        seed_score = params.assembly_seedscore
        min_pls_len = params.assembly_minplaslen
    }

    base = metaname(meta)
    bins = "${base}.pbf.bins.tsv"
    pflow = "$projectDir/PlasBin-flow-pangenome/code/plasbin_flow.py"

    """
    bgzip -k ${gfa}
    python ${pflow} ${args} -alpha4 ${meta.pty} -ag ${gfa}.gz -gc ${gc} -out_dir . -out_file ${bins}  -score ${plscore} -assembler ${assembler} -seed_len ${seed_len}  -seed_score ${seed_score} -min_pls_len ${min_pls_len}
    """

    stub:
    base = metaname(meta)
    bins = "${base}.pbf.bins.tsv"
    """
    touch ${bins}
    """

    
}

process transform {
    cache false
    input:
    tuple val(meta), path(bins), path(gfa)

    output:
    tuple val(meta), path(res)
    
    script:

    res = metaname(meta) + ".pred.tab"
    
    """
    python $projectDir/bin/evaluation/transform_pbf_pred.py --input ${bins} --gfa ${gfa} --output ${res} 
    """
} 

process transform_nve {
    cache false

    input:
    tuple val(meta), path(pred), path(graph)

    output:
    tuple val(meta), path(modded)
    
    script:
        modded = metaname(meta) + ".nve.tab"
        """
        #!/usr/bin/env bash
        python $projectDir/bin/extend_bins.py --pred ${pred} --out ${modded} --naive --n 1
        """
}
process transform_ovl {
    cache false

    input:
    tuple val(meta), path(pred), path(graph)

    output:
    tuple val(meta), path(modded)
    
    script:
        modded = metaname(meta) + ".ovl.tab"
        """
        #!/usr/bin/env bash
        python $projectDir/bin/extend_bins.py --pred ${pred} --out ${modded} --super --graph ${graph}
        """
}

process labeling {
    cache false

    input:
    tuple val(meta), path(prediction), path(gt)

    output:
    tuple val(meta), path(stats)

    script:

    output = prediction.getBaseName()
    id = output.split("\\.")[0]
    asm = output.split("\\.")[1]
    thr = output.split("\\.")[2]
    tool = output.split("\\.")[3]

    reference = gt.getBaseName().split("\\.")[2]

    description = "sample ${id}, ${asm} graph (cut ${thr}, ref ${reference}) with ${tool}"
    stats = "${output}.${reference}.lab.txt"

    """
    python $projectDir/bin/evaluate_bins.py --bin ${prediction} --csv ${gt} --sample ${output} --output ${stats} --description '${description}'
    """
}

process binning { 
    cache false

    input:
    tuple val(meta), path(prediction), path(gt)

    output:
    tuple val(meta), path(stats)
    // tuple val(meta), path(plots), emit: plots

    script:
    output = prediction.getBaseName()
    reference = gt.getBaseName().split("\\.")[2]
    stats = "${output}.${reference}.bin.txt"
    // awk -i inplace '{\$0=gensub(/\s*\S+/,\\"\\",3)}1' ${gt} 
    // awk -i inplace '{\$0=gensub(/\s*\S+/,\\"\\",4)}1' ${gt} 
    """
    #!/bin/bash
    echo 'plasmid\tcontig\tnan\tnan\tcontig_len' | cat - ${gt} > temp && mv temp ${gt}
    python $projectDir/PlasEval/src/plaseval.py eval --pred ${prediction} --gt ${gt} --out ${stats} 
    """

}
workflow { // single input workflow, for dataset input use pipeline.nf

    // params input is the dataset folder
    input_ch = Channel.fromFilePairs("${params.dataset}/*-S*-{s,u}.gfa.gz", type: 'file', checkIfExists: true)

    input_ch = input_ch.map{
        m, filepair ->
        def fmeta = [:]
        fmeta.id = m.split("-")[1]
        fmeta.spc = m.split("-")[0]
        [fmeta, filepair]
    }

    meta_ch = input_ch.map { meta, _files -> meta }
    ref_ch = Channel.empty()
    meta_filt_ch = meta_ch | gt | refcheck.tap{ref_ch} | map {meta, fna -> [meta]}

    input_filt_ch = meta_filt_ch.join(input_ch)


    ske_ch = input_filt_ch.map { meta, files -> meta += [asm: "s"]; def f = files[0] ; [meta, f]}
    // | view
    uni_ch = input_filt_ch.map { meta, files -> meta += [asm: "u"]; def f = files[1] ; [meta, f]}
    // | view

    // |view
    asm_ch = ske_ch.concat(uni_ch)


// from meta retrive GT


    gfa_ch = asm_ch
    | bgzip_d // decompress gfa.gz
    | rename_gfa 
    | combine(Channel.fromList(cutl)) // add cut parameter
    | map {_m, _f, _c -> _m['cut'] = _c; [_m, _f]} // map to meta
    | remove_nodes

    fasta_ch = gfa_ch
    | extract_fasta

    fasta_ch_s = fasta_ch
    | filter {meta, f -> meta['asm'] == 's'}
    | map {meta, f -> def m = [:]; m['id'] = meta['id']; m['spc'] = meta['spc']; m['cut'] = meta['cut']; [m, f]}

    fasta_ch_u = fasta_ch
    | filter {meta, f -> meta['asm'] == 'u'}
    | map {meta, f -> def m = [:]; m['id'] = meta['id']; m['spc'] = meta['spc']; m['cut'] = meta['cut']; [m, f]}

    pangenome_ch = Channel.empty()
    panassembly_ch = Channel.empty()
    
    gfa_ch_s = gfa_ch
    | filter {meta, f -> meta['asm'] == 's' }
    | map {meta, f -> def m = [:]; m['id'] = meta['id']; m['spc'] = meta['spc']; m['cut'] = meta['cut']; [m, f]}
    | combine(Channel.fromList(pctid))
    | combine(Channel.fromList(thr))
    | map {meta, f, _p, _t -> meta += ['pctid': _p, 'thr': _t] ;[meta, f] }
    
    gfa_ch_u = gfa_ch
    | filter {meta, f -> meta['asm'] == 'u' }
    | map {meta, f -> def m = [:]; m['id'] = meta['id']; m['spc'] = meta['spc']; m['cut'] = meta['cut']; [m, f]}
    | combine(Channel.fromList(pctid))
    | combine(Channel.fromList(thr))
    | map {meta, f, _p, _t -> meta += ['pctid': _p, 'thr': _t] ;[meta, f] }
    
    fastamix_ch = fasta_ch_s.join(fasta_ch_u)

    pangenome_ch = fastamix_ch
    | mix_fasta
    | combine(Channel.fromList(pctid))
    | combine(Channel.fromList(thr))
    | map {meta, f, _p, _t -> meta += [pctid: _p, thr: _t] ;[meta, f] } 
    | make_pangenome
    
    panassembly_ch = pangenome_ch
    | join(gfa_ch_s)
    | join(gfa_ch_u)
    | make_panassembly


    scores_ch = panassembly_ch 
    | join(gfa_ch_s)
    | join(gfa_ch_u)
    | compute_scores


    pangebin_bin_ch = panassembly_ch
    | join(compute_scores.out.gc_pan)
    | join(compute_scores.out.gd_pan)
    | model

    pangebin_pred_ch = pangebin_bin_ch
    | join(panassembly_ch)
    | transform

    pangebin_naive_pred_ch = pangebin_pred_ch
    | join(panassembly_ch)
    | transform_nve
    
    pangebin_ovl_pred_ch = pangebin_pred_ch
    | join(panassembly_ch)
    | transform_ovl


    pred_ch = pangebin_pred_ch
    | concat(pangebin_naive_pred_ch)
    | concat(pangebin_ovl_pred_ch)
    // | combine()

return
    pred_ch
    | labeling

    pred_ch
    | binning


    // blast_ch_s
    // blast_ch_u


// make pangenome


}
