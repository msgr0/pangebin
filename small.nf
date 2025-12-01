def metaname(meta) {
    def name = meta.id
    if (meta.containsKey('spc')) {
        name += "-${meta.spc}"
    }
    if (meta.containsKey('asm')) {
        name += "-${meta.asm}"
    }
    if (meta.containsKey('cut')) {
        name += "-c${meta.cut}"
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
    if (meta.containsKey('bintype')) {
        name += "-${meta.bintype}"
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
def pty = [0, 0.5, 1, 2]
def asmlist = ["u", "s"]

// def bintypeL = ["bins", "nve", "ovl"]
def bintypeL = ["bins"]
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
    bgzip -k -d ${gz}
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
    storeDir null
    errorStrategy 'ignore'
    input:

    tuple val(meta), path(refin)

    output:
    tuple val(meta), path(ref), optional: true

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

process blast {

    errorStrategy 'terminate'
    cache false
    maxForks 32

    input:
    tuple val(meta), path(graph), path(mix), path(ske), path(uni), path(reference)

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

    stub:

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
    touch ${pan_mix_gt}
    touch ${mix_gt}
    touch ${pan_ske_gt}
    touch ${ske_gt}
    touch ${pan_uni_gt}
    touch ${uni_gt}
    """


}
//................//
process make_pangenome {
    if (params.executor == 'slurm') {
        executor 'slurm'
        memory '16 GB'
        array 60
	maxForks 180
        cpus 8
        time { task.attempt > 1 ? 8.hour : 2.hour}
    } else {
        executor 'local'
        memory '16 GB'
        cpus 8
    }
    errorStrategy 'retry'
    maxRetries 2

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
    
    bgzip -k ${mixed_fasta}
    mkdir -p ./work-tmp
    nextflow run nf-core/pangenome -r ${release} -profile $profile -w ./work-tmp  --input ${mixed_fasta}.gz --n_haplotypes ${haplos} --outdir ${out_core} -params-file ${paramfile} --wfmash_segment_length ${meta.thr} --wfmash_map_pct_id ${meta.pctid}
    cp ${out_core}/FINAL_GFA/${mixed_fasta}.gz.gfaffix.unchop.Ygs.view.gfa ${pangenome} 
    rm -rf ${out_core}
    rm -rf ./work-tmp

    """

    stub:
    pangenome = metaname(meta) + ".pan.gfa"
    """
    touch ${pangenome}
    """
}

process make_panassembly {

      maxForks 8
//     executor 'slurm'
//     memory '8 GB'
//     array 60
//     cpus 1
//     time '15m'

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
    if (params.executor == 'slurm') {   
        executor 'slurm'
        array 60
        maxForks 120
        cpus 16
        time { task.attempt > 1 ? 4.hour * task.attempt : 4.hour}
        memory { task.attempt > 1 ? (task.previousTrace.memory > (16.GB) ? (task.previousTrace.memory * 2) : (32.GB)) : (32.GB) }
    }
    else {
        executor 'local'
        maxForks 1
        cpus 8
        memory '16 GB'
    }

    errorStrategy 'ignore'

    input:
    tuple val(meta), path(gfa), path(gc), path(plscore)

    
    output:
    tuple val(meta), path(bins)

    script:

    assembler = "pangenome"
    args = ""
    seed_len = params.pangenome_seedlen
    seed_score = params.pangenome_seedscore
    min_pls_len = params.pangenome_minplaslen
    pty_local = meta.pty

    if (pty_local == 0) {
      args = "--nopenalty"
    }
    
    if (meta.asm != null && meta.asm != "pan" && meta.asm != "p") { 
        assembler = "unicycler"
        args = "--nopenalty"
        pty_local = 0
        seed_len = params.assembly_seedlen
        seed_score = params.assembly_seedscore
        min_pls_len = params.assembly_minplaslen
    }

    base = metaname(meta)
    bins = "${base}.pbf.bins.tsv"
    pflow = "$projectDir/PlasBin-flow-pangenome/code/plasbin_flow.py"

    if (task.attempt > 3) {
    """
    touch ${bins}
    """
    }

    """
    bgzip -k ${gfa}
    python ${pflow} ${args} -alpha4 ${pty_local} -ag ${gfa}.gz -gc ${gc} -out_dir . -out_file ${bins}  -score ${plscore} -assembler ${assembler} -seed_len ${seed_len} -seed_score ${seed_score} -min_pls_len ${min_pls_len} -log_file loglog.log
    
    rm -rf loglog.log
    rm -rf m.log
    
    """

    stub:
    base = metaname(meta)
    bins = "${base}.pbf.bins.tsv"
    """
    touch ${bins}
    """

    
}

process check_model {

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path(bins2), optional: true

    script:
    bins2 = metaname(meta) + ".pbf.binsx.tsv"
    """
    #!/usr/bin/env bash

    if [[ \$(wc -l <${bins}) -ge 2 ]]
    then
        cp ${bins} ${bins2}
    else
        echo "No bins predicted, outputting no file"
    fi

    """

}

process transform {
    maxForks 24
    input:
    tuple val(meta), path(bins), path(gfa)

    output:
    tuple val(meta), path(res)
    
    script:

    meta += ['bintype': 'bins']
    res = metaname(meta) + ".tab"
    
    """
    python $projectDir/bin/evaluation/transform_pbf_pred.py --input ${bins} --gfa ${gfa} --output ${res} 
    
    """

    stub:
    meta += ['bintype': 'bins']
    res = metaname(meta) + ".tab"
    """
    touch ${res}
    """
} 

process transform_nve {

    input:
    tuple val(meta), path(pred), path(graph)

    output:
    tuple val(meta), path(modded), optional: true
    
    script:
    meta += ['bintype': 'nve']
        modded = metaname(meta) + ".tab"
        """
        #!/usr/bin/env bash
        if [[ \$(wc -l <${pred}) -ge 2 ]];
        then
            python $projectDir/bin/extend_bins.py --pred ${pred} --out ${modded} --naive --n 1
        fi
        """

    stub:
    meta += ['bintype': 'nve']
        modded = metaname(meta) + ".tab"
        """
        touch ${modded}
        """


}
process transform_ovl {

    input:
    tuple val(meta), path(pred), path(graph)

    output:
    tuple val(meta), path(modded), optional: true
    
    script:
    meta += ['bintype': 'ovl']
        modded = metaname(meta) + ".tab"
        """
        #!/usr/bin/env bash
        if [[ \$(wc -l <${pred}) -ge 2 ]];
        then
            python $projectDir/bin/extend_bins.py --pred ${pred} --out ${modded} --super --graph ${graph}
        fi
        """
        stub:
    meta += ['bintype': 'ovl']
        modded = metaname(meta) + ".tab"
        """
       touch ${modded}
       """

}

process labeling {

    input:
    tuple val(meta), path(prediction), path(gt)

    output:
    tuple val(meta), path(stats)

    script:

    meta += ['t': 'l']
    name = metaname(meta)

    description = "sample ${name}"
    stats = metaname(meta) + ".txt"

    """
    python $projectDir/bin/evaluate_bins.py --bin ${prediction} --csv ${gt} --sample ${name} --output ${stats} --description '${description}'
    """
    stub:
    meta += ['t': 'l']
    stats = metaname(meta) + ".txt"
    """
    touch ${stats}
    """
}

process binning { 

    input:
    tuple val(meta), path(prediction), path(gt)

    output:
    tuple val(meta), path(stats)

    script:
    meta += ['t': 'b']

    stats = metaname(meta) + ".txt"
    // awk -i inplace '{\$0=gensub(/\s*\S+/,\\"\\",3)}1' ${gt} 
    // awk -i inplace '{\$0=gensub(/\s*\S+/,\\"\\",4)}1' ${gt} 
    """
    #!/bin/bash
    echo 'plasmid\tcontig\tnan\tnan\tcontig_len' | cat - ${gt} > temp && mv temp ${gt}
    python $projectDir/PlasEval/src/plaseval.py eval --pred ${prediction} --gt ${gt} --out ${stats} 
    """
    stub:
    meta += ['t': 'b']
    stats = metaname(meta) + ".txt"

"""
touch ${stats}
"""
}

// process extract_results {
//     storeDir null
//   cache false

//   input:
//   tuple val(meta), path(res)

//   output:
//   tuple val(meta),  env(prec), env(rec), env(f1)//   val(prec), val(rec), val(f1)


//   script:
//   def prec = 0
//   def rec = 0
//   def f1 = 0

//   if (meta.t == 'l') {

//     """
// if [ \$(wc -l < ${res}) -gt 1 ]; then
//     tail -n 4 "${res}" | head -n 3 | awk 'NR==1{prec=\$NF} NR==2{rec=\$NF} NR==3{f1=\$NF} END{print prec, rec, f1}'
//     export prec
//     export rec
//     export f1
// else
//     echo "File has one line or less"
//     export prec=0
//     export rec=0
//     export f1=0
// fi

//     """
// }
// else if (meta.t == 'b') {

//   """
// if [ \$(wc -l < ${res}) -gt 1 ]; then
//     tail -n 4 "${res}" | head -n 3 | awk 'NR==1{prec=\$NF} NR==2{rec=\$NF} NR==3{f1=\$NF} END{print prec, rec, f1}'
//     export prec
//     export rec
//     export f1
// else
//     echo "File has one line or less"
//     export prec=0
//     export rec=0
//     export f1=0
// fi
//     """
// }
// else {
//     """
//         echo "Unknown type"
//         export prec=0
//         export rec=0
//         export f1=0
//     """

// }


//   stub:
//   """
//   export prec = 0.9
//     export rec = 0.4
//     export f1 = 0.4
//   """

// }



workflow { // single input workflow, for dataset input use pipeline.nf
    main:
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

    ref_ch = meta_ch | gt | refcheck
      
    meta_filt_ch = ref_ch
    | map {meta, fna -> [meta]}

    input_filt_ch = meta_filt_ch.join(input_ch)


    ske_ch = input_filt_ch.map { meta, files -> meta += [asm: "s"]; def f = files[0] ; [meta, f]}
    // | view
    uni_ch = input_filt_ch.map { meta, files -> meta += [asm: "u"]; def f = files[1] ; [meta, f]}
    // | view

    // |view
    asm_ch = ske_ch.concat(uni_ch)


// from meta retrive GT


    gfa_ch = asm_ch
    | bgzip_d
    | rename_gfa
    | combine(Channel.fromList(cutl))
    | map { meta, file, cutvalue -> meta += ["cut": cutvalue]; [meta, file] } 
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


    gfa_ch_s_pan = gfa_ch_s 
    | combine(Channel.fromList(pctid))
    | combine(Channel.fromList(thr))
    | map {meta, f, _p, _t -> meta += ['pctid': _p, 'thr': _t] ;[meta, f] }
    
    gfa_ch_u = gfa_ch
    | filter {meta, f -> meta['asm'] == 'u' }
    | map {meta, f -> def m = [:]; m['id'] = meta['id']; m['spc'] = meta['spc']; m['cut'] = meta['cut']; [m, f]}

    gfa_ch_u_pan = gfa_ch_u
    | combine(Channel.fromList(pctid))
    | combine(Channel.fromList(thr))
    | map {meta, f, _p, _t -> meta += ['pctid': _p, 'thr': _t] ;[meta, f] }
    
    fastamix_ch = fasta_ch_s.join(fasta_ch_u)
    | mix_fasta

    pangenome_ch = fastamix_ch
    | combine(Channel.fromList(pctid))
    | combine(Channel.fromList(thr))
    | map {meta, f, _p, _t -> meta += ['pctid': _p, 'thr': _t] ;[meta, f] }
    | make_pangenome
    
    panassembly_ch = pangenome_ch
    | join(gfa_ch_s_pan)
    | join(gfa_ch_u_pan)
    | make_panassembly


    scores_ch = panassembly_ch 
    | join(gfa_ch_s_pan)
    | join(gfa_ch_u_pan)
    | compute_scores


    panassembly_model_ch = panassembly_ch
    | join(compute_scores.out.gc_pan)
    | join(compute_scores.out.gd_pan)
    | combine(Channel.fromList(pty))
    | map { meta, file, gc, pls, _pty -> meta += ['pty': _pty]; [meta, file, gc, pls]}
    
    // panassembly_gt_ch = panassembly_ch

    // panassembly_ch = panassembly_model_ch 
    // | map { meta, file, gc, pls -> [meta, file]}
    

    skesa_model_ch = gfa_ch 
    | filter {meta, f -> meta['asm'] == 's'}
    | join( (compute_scores.out.gc_ske).filter{meta, f -> meta['thr'] == thr[0] && meta['pctid'] == pctid[0]}.map {meta, f -> def m = [:]; m['id'] = meta['id']; m['spc'] = meta['spc']; m['asm'] = 's'; m['cut'] = meta['cut']; [m, f]})
    | join( (compute_scores.out.gd_ske).filter{meta, f -> meta['thr'] == thr[0] && meta['pctid'] == pctid[0]}.map {meta, f -> def m = [:]; m['id'] = meta['id']; m['spc'] = meta['spc']; m['asm'] = 's'; m['cut'] = meta['cut']; [m, f]})

    uni_model_ch = gfa_ch
    | filter {meta, f -> meta['asm'] == 'u' }
    | join( (compute_scores.out.gc_uni).filter{meta, f -> meta['thr'] == thr[0] && meta['pctid'] == pctid[0]}.map {meta, f -> def m = [:]; m['id'] = meta['id']; m['spc'] = meta['spc']; m['asm'] = 'u'; m['cut'] = meta['cut']; [m, f]})
    | join( (compute_scores.out.gd_uni).filter{meta, f -> meta['thr'] == thr[0] && meta['pctid'] == pctid[0]}.map {meta, f -> def m = [:]; m['id'] = meta['id']; m['spc'] = meta['spc']; m['asm'] = 'u'; m['cut'] = meta['cut']; [m, f]})



    bin_ch = panassembly_model_ch
    | mix (skesa_model_ch)
    | mix (uni_model_ch)
    

    model_ch = bin_ch
    | model
    // check for empty predictions and remove them
    // model_ok_ch = model_ch
    // | check_model
    // | view

    model_pan_ch = model_ch
    | join (panassembly_ch.combine(Channel.fromList(pty)).map{meta, file, _pty -> meta += ['pty': _pty]; [meta, file]})

    model_asm_ch = model_ch
    | join (gfa_ch)


// gfa_ch = (meta(asm, cutln, ...), gfa)

    // transform predictions to required format
    transform_ch = model_pan_ch
    | mix (model_asm_ch)
    | transform


    // transform_nve_ch = model_pan_ch
    // | transform_nve

    // transform_ovl_ch = model_pan_ch
    // | transform_ovl

    prediction_ch = transform_ch
    // | mix(transform_ovl_ch)
    // | mix(transform_nve_ch)


    fasta_join_ch = fastamix_ch
    | join (fasta_ch_s)
    | join (fasta_ch_u)
    | map {m, f1, f2, f3 -> [m['id'], m, f1, f2, f3]}
    | combine (ref_ch.map{m, f -> [m['id'], f]}, by:0 )
    | map {mid, m, f1, f2, f3, f4 -> [mid, m['cut'], f1, f2, f3, f4]}


    blast_ch = panassembly_ch.map {meta, f -> [meta['id'], meta['cut'], meta, f]}
    | combine(fasta_join_ch, by:[0,1])
    | map {mid, mcut, meta, f1, f2, f3, f4, f5 -> [meta, f1, f2, f3, f4, f5]}
    | blast


    prediction_pan_ch = prediction_ch.filter { meta, pred -> meta['pty'] != null }
    prediction_asm_ch = prediction_ch.filter { meta, pred -> meta['pty'] == null }
    // (meta, bins.tab)
    // labeling()
    // tuple val(meta), path(prediction), path(gt)
    
    // binning()
    // tuple val(meta), path(prediction), path(gt)

    // gt()
    // tuple val(meta), path(pan_mix_gt), path(mix_gt), emit: pangt
    // tuple val(meta), path(pan_ske_gt), path(ske_gt), emit: skegt
    // tuple val(meta), path(pan_uni_gt), path(uni_gt), emit: unigt

    // blast.out.skegt | view
    // meta[id, spc, cut, pctid, thr], pangt, gt

    blast_asm_ch = blast.out.skegt.filter{meta, pan, gt -> meta['thr']==1000 && meta['pctid']==95}.map{meta, pan, gt-> meta+=['asm':'s']; [meta, gt]}
    | mix (
       blast.out.unigt.filter{meta, pan, gt -> meta['thr']==1000 && meta['pctid']==95}.map{meta, pan, gt-> meta+=['asm':'u']; [meta, gt]} 
    )
    blast_pan_ch = blast.out.skegt.map{meta, pan, gt-> meta+=['asm':'s']; [meta, pan]}
    | mix (
       blast.out.unigt.map{meta, pan, gt-> meta+=['asm':'u']; [meta, pan]} 
    )


    scores_asm_ch = 
    prediction_asm_ch.map {meta, bins -> [meta['id'], meta['asm'], meta['cut'], meta, bins]}
    | join(blast_asm_ch.map{meta, gt -> [meta['id'], meta['asm'], meta['cut'], gt]}, by: [0,1,2])
    | map {_id, _asm, _cut, meta, bins, gt -> [meta, bins, gt]}

    scores_pan_ch = 
    prediction_pan_ch.combine(Channel.fromList(asmlist)).map{meta, bins, _asm -> [meta['id'], meta['cut'], meta['pctid'], meta['thr'], _asm, meta, bins]}
    | combine(blast_pan_ch.map{meta, gt -> [meta['id'], meta['cut'], meta['pctid'], meta['thr'], meta['asm'], gt]}, by: [0,1,2,3,4])
    | map {_id, _cut, _pctid, _thr, _asm, meta, bins, gt -> meta+=['asm':_asm]; [meta, bins, gt]}

    results_ch = 
    scores_asm_ch
    | mix (scores_pan_ch)
    | (labeling & binning)


}

