#!/usr/bin/env nextflow

// def convert_species(spec) {
//   switch(spec) {
//     case 'efea':
//       species = 'Enterococcus faecium'
//       break
//     case 'kpne':
//       species = 'Klebsiella pneumoniae'
//       break
//     case 'abau':
//       species = 'Acinetobacter baumannii'
//       break
//     case 'ecol':
//       species = 'Escherichia coli'
//       break
//     default:
//       species = 'Other species'
//   }
//   return species
// }

process extractGfa {
    storeDir "${params.input}/"
    cache 'lenient'

    input:
    tuple val(meta), path(gfa)
    
    output:
    tuple val(meta), path(gfa_std)

    script:
    if (meta.asm == 's') { // transform skesa graph into unicycler
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
    #!/usr/bin/env bash
    bgzip -d ${gfa}
    python $projectDir/bin/rename_gfa.py -i ${gfa_xtracted} -o ${gfa_std} -p ${asm} ${param}
    """
}

process removeNodes {
    storeDir "${params.input}/"
    cache 'lenient'
    
    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path(trimmed), emit: gfa

    script:
    trimmed = gfa.getBaseName() + ".${meta.thr}.gfa"

    if (meta.thr == 0 )
    """
    #!/usr/bin/env bash
    mv ${gfa} ${trimmed}
    """
    else
    """
    #!/usr/bin/env bash
    python $projectDir/bin/remove_nodes.py -v 2 -i ${gfa} -o ${trimmed} -t ${meta.thr}
    """
}

process extractFasta {
    storeDir "${params.input}/"
    cache 'lenient'

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val (meta), path(fasta), emit: fasta

    script:
    fasta = gfa.getBaseName() + ".fasta"

    """
    #!/usr/bin/env bash
    awk '/^S/{print ">"\$2; print ""\$3}' ${gfa} > ${fasta}    
    """
}

process mixFasta {
    storeDir "${params.input}/"
    cache 'lenient'

    input:
    tuple val(meta), path(fasta_a), path(fasta_b)
    output:
    tuple val(meta), path(fasta)

    script:
    fasta = "${meta.id}.${meta.thr}.fa"
    """
    cat ${fasta_a} ${fasta_b} > ${fasta}
    """
}

process makePangenome {
    storeDir "${params.input}/"
    maxForks 4
    time '60m'
    cache 'lenient'

    input:
    tuple val(meta), path(mixed_fasta)
    
    output:
    tuple val(meta), path(pangenome), emit: pangenome

    script:
    pangenome = "${meta.id}.pan.gfa"
    out_core = "${meta.id}_nfcore"
    haplos = 2
    paramfile = "${meta.id}.params.json"
    release = '1.1.2'
    profile = 'apptainer'

    """
    #!/usr/bin/env bash
    touch ${paramfile}
    echo -e '{\n\t"max_memory": "16.GB",\n\t"wfmash_segment_length": $params.minlen,\n\t"seqwish_min_match_length": 0,\n\t"smoothxg_poa_params": "asm5",\n\t"wfmash_map_pct_id": $params.pctid,\n\t"max_cpus": 16,\n\t"wfmash_merge_segments": true\n}' > ${paramfile}


    bgzip ${mixed_fasta}
    nextflow run nf-core/pangenome -r ${release} -profile $profile -resume --input ${mixed_fasta}.gz --n_haplotypes ${haplos} --outdir ${out_core} -params-file ${paramfile}
    cp ${out_core}/FINAL_GFA/${mixed_fasta}.gz.gfaffix.unchop.Ygs.view.gfa ${pangenome}
    """
    // bgzip /data/proj/pangebin/data-test/SAMN02786856/SAMN02786856.1.fa 
    // nextflow run nf-core/pangenome -r 1.1.2 -profile podman -resume --input /data/proj/pangebin/data-test/SAMN02786856/SAMN02786856.1.fa.gz --n_haplotypes 2 --outdir test-nfcore.txt -params-file pangenome-params.json
}

process makePanassembly {
    storeDir "${params.input}/"
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
    #!/usr/bin/env bash

    python $projectDir/bin/gfa_cleaner.py --input ${pangenome} --output ${pangenomecl}
    python $projectDir/bin/pangenome_enancher.py --input ${pangenomecl} --output ${panassembly} --assembly ${uni} ${ske}
    """
}


process computeScores {
    storeDir "${params.input}/"
    cache 'lenient'

    input:
    tuple val(meta), path(skesa), path(unicycler)

    output:
    tuple val(meta), path(gc_ske), emit: gc_ske
    tuple val(meta), path(gd_ske), emit: gd_ske
    
    tuple val(meta), path(gc_uni), emit: gc_uni
    tuple val(meta), path(gd_uni), emit: gd_uni
    

    script:
    base_ske = "${meta.id}.ske.${meta.thr}"
    base_uni = "${meta.id}.uni.${meta.thr}"

    gd_ske = base_ske + ".gd.tsv"
    gd_uni = base_uni + ".gd.tsv"

    gc_ske = base_ske + ".gc.tsv"
    gc_uni = base_uni + ".gc.tsv"

    input_csv = "${meta.id}.${meta.thr}.input.csv"

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
    base_pan = "${meta.id}.pan.${meta.thr}"
    base_ske = "${meta.id}.ske.${meta.thr}"
    base_uni = "${meta.id}.uni.${meta.thr}"

    gd_pan = base_pan + ".gd.tsv"
    gd_ske = base_ske + ".gd.tsv"
    gd_uni = base_uni + ".gd.tsv"

    gc_pan = base_pan + ".gc.tsv"
    gc_ske = base_ske + ".gc.tsv"
    gc_uni = base_uni + ".gc.tsv"

    input_csv = "${meta.id}.${meta.thr}.input.csv"

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



    cat ${gd_ske} ${gd_uni} > ${meta.id}.${meta.thr}.mix.gd.tsv
    cat ${gc_ske} ${gc_uni} > ${meta.id}.${meta.thr}.mix.gc.tsv

    python3 $projectDir/bin/pbf_preprocess.py --input ${pangenome} --gd  ${meta.id}.${meta.thr}.mix.gd.tsv --gc ${meta.id}.${meta.thr}.mix.gc.tsv --output ${meta.id}.pan.${meta.thr}
  
    """
}

// process computeMLscores {

// }

// process computeAllMLscores {
//     cache 'lenient'

//     input:
//     tuple val(meta), path(gpan), path(gske), path(guni)
//     tuple val(meta2), path(fpan), path(fske), path(funi)
//     val(species)

//     output:
//     tuple val(meta), path(pred_pan), emit: pan
//     tuple val(meta), path(pred_uni), emit: uni
//     tuple val(meta), path(pred_ske), emit: ske
    
//     script:
//     mlplas_threshold = '0.5'
//     species = convert_species("${species}")
//     pred_pan = fpan.getBaseName() + ".mlp.tab"
//     pred_uni = funi.getBaseName() + ".mlp.tab"
//     pred_ske = fske.getBaseName() + ".mlp.tab"
//     """
//     #!/bin/bash
//     Rscript $projectDir/bin/run_mlplasmids.R ${fpan} ${pred_pan} ${mlplas_threshold} '${species}' TRUE
//     Rscript $projectDir/bin/run_mlplasmids.R ${funi} ${pred_uni} ${mlplas_threshold} '${species}' TRUE
//     Rscript $projectDir/bin/run_mlplasmids.R ${fske} ${pred_ske} ${mlplas_threshold} '${species}' TRUE

//     python $projectDir/bin/mlpl.asmtopan.py --pred ${pred_pan} --graph ${gfa}  --pbf ${pred_pbf}
//     python $projectDir/bin/mlpl.asm.py --pred ${pred_uni} --graph ${gfa}  --pbf ${pred_ske}
//     python $projectDir/bin/mlpl.asm.py --pred ${pred_ske} --graph ${gfa}  --pbf ${pred_uni}

//     """
//     // python $projectDir/bin/mlpl.asm.py --pred ${pred} --graph ${gfa}  --output ${pred_gplas}

//     input:
//     tuple val(meta), path(pangenome), path(skesa), path(unicycler)

//     output:
//     tuple val(meta), path(gc_pan), emit: gc_pan
//     tuple val(meta), path(gd_pan), emit: gd_pan
//     tuple val(meta), path(gc_ske), emit: gc_ske
//     tuple val(meta), path(gd_ske), emit: gd_ske
//     tuple val(meta), path(gc_uni), emit: gc_uni
//     tuple val(meta), path(gd_uni), emit: gd_uni
    

//     script:
//     base_pan = "${meta.id}.pan.${meta.thr}"
//     base_ske = "${meta.id}.ske.${meta.thr}"
//     base_uni = "${meta.id}.uni.${meta.thr}"

//     gd_pan = base_pan + ".gd.tsv"
//     gd_ske = base_ske + ".gd.tsv"
//     gd_uni = base_uni + ".gd.tsv"

//     gc_pan = base_pan + ".gc.tsv"
//     gc_ske = base_ske + ".gc.tsv"
//     gc_uni = base_uni + ".gc.tsv"

//     input_csv = "${meta.id}.${meta.thr}.input.csv"

//     putils = "$projectDir/PlasBin-flow-pangenome/code/plasbin_utils.py"
//     pls_db_file = "$projectDir/PlasBin-flow-pangenome/database/genes.fasta"

//     """
//     #!/usr/bin/env bash

//     bgzip -k ${skesa}
//     bgzip -k ${unicycler}

//     echo "sample,gfa" > "${input_csv}"
//     echo "${base_ske},${skesa}.gz" >> "${input_csv}"
//     echo "${base_uni},${unicycler}.gz" >> "${input_csv}"

//     python ${putils} preprocessing --input_file ${input_csv} --out_dir . --tmp_dir ./tmp --out_file out_file.csv --db_file ${pls_db_file}



//     cat ${gd_ske} ${gd_uni} > ${meta.id}.${meta.thr}.mix.gd.tsv
//     cat ${gc_ske} ${gc_uni} > ${meta.id}.${meta.thr}.mix.gc.tsv

//     python3 $projectDir/bin/pbf_preprocess.py --input ${pangenome} --gd  ${meta.id}.${meta.thr}.mix.gd.tsv --gc ${meta.id}.${meta.thr}.mix.gc.tsv --output ${meta.id}.pan.${meta.thr}
  
//     """
// }


workflow PREPROCESS {
    take:
    input_ch // short reads assembled by unicycler and skesa in 2 .gfa.gz

    main:
    extracted_gfa = extractGfa(input_ch)
    extracted_gfa.map{meta, files -> meta += [thr:params.cutlen]; [meta, files]}
    | removeNodes
    | extractFasta

    skesaGfa_ch = removeNodes.out.gfa.filter { meta, file -> meta['asm'] == 's' }.map{
        meta, files -> meta = [id: meta.id, thr:meta.thr]; [meta, files]
    }

    uniGfa_ch = removeNodes.out.gfa.filter { meta, file -> meta['asm'] == 'u' }.map{
        meta, files -> meta = [id: meta.id, thr:meta.thr]; [meta, files]
    }

    skesaFasta_ch = extractFasta.out.fasta.filter { meta, file -> meta['asm'] == 's' }.map{
        meta, files -> meta = [id: meta.id, thr:meta.thr]; [meta, files]
    } | view

    uniFasta_ch = extractFasta.out.fasta.filter { meta, file -> meta['asm'] == 'u' }.map{
        meta, files -> meta = [id: meta.id, thr:meta.thr]; [meta, files]
    } | view

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
        mixFasta_ch = skesaFasta_ch.combine(uniFasta_ch, by: 0)
        | mixFasta
     
        pangenome_ch = mixFasta_ch | makePangenome
        
        panasmGfa_ch = makePanassembly(pangenome_ch.join(skesaGfa_ch).join(uniGfa_ch)) 
        result = computeAllScores(panasmGfa_ch.join(skesaGfa_ch).join(uniGfa_ch))
        gcPan_ch = gcPan_ch.mix(result.gc_pan)
        gdPan_ch = gdPan_ch.mix(result.gd_pan)
        // if (params.ml) {
        //     mlres = computeAllMLscores(panasmGfa_ch.join(skesaGfa_ch).join(uniGfa_ch))
        //     mlPan_ch = mlPan_ch.mix(mlres.pan)
        //     mlSke_ch = mlSke_ch.mix(mlres.ske)
        //     mlUni_ch = mlUni_ch.mix(mlres.uni)


        // }

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
