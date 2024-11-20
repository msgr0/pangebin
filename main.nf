#!/usr/bin/env nextflow

/*  PANGEBIN
----------------------------------------------
gtihub: https://github.com/msgr0/pangebin.git
----------------------------------------------

1. preprocessing
2. pangenome
3. pbf pangneome model 
optional. a. ground-truth
          b. evaluation */

include { GT                   } from './workflows/gt'
include { PREPROCESS           } from './workflows/preprocess'
include { PANGENOME            } from 'nf-core/pangenome/workflows/pangenome'
include { PAN                  } from './workflows/pangenome'
include { AUGMENT              } from './workflows/augment'
include { PBF as MODEL         } from './workflows/pbf'
include { PBF as MODELPAN      } from './workflows/pbf'
include { EVAL                 } from './workflows/eval'
include { PLASEVAL             } from './workflows/eval'


workflow {

    input_ch = Channel.fromFilePairs("${params.input}/*-S*-{s,u}.gfa.gz", type: 'file', checkIfExists: true)

    input_ch = input_ch.map{ meta, filepair -> def fmeta = [:]; fmeta.id = meta[5..-1]; fmeta.species = meta[0..3]
        [fmeta, filepair]
    }
    meta_ch = input_ch.map{meta, files -> meta}

    unicycler_ch = input_ch | map { meta, files -> meta += [asm: "uni"]; uni = files.findAll { it.toString().contains("-u.gfa.gz")}; [meta, uni[0]]}
    skesa_ch = input_ch | map { meta, files -> meta += [asm: "ske"]; ske = files.findAll { it.toString().contains("-s.gfa.gz")}; [meta, ske[0]]}

    // skesa_ch | view
    PREPROCESS(unicycler_ch.mix(skesa_ch))
    PUBLISH(PREPROCESS.out.all)



    ch_input = Channel.empty()
    ch_input = Channel.fromFilePairs("${params.input}/*-S*-{s,u}.gfa.gz", type:'file', checkIfExists: true)
    if (ch_input.isEmpty())

    if (params.strict) { // run only on the samples that give a result in the GT
        GT()
        PREPROCESS()
    }

    else if(params.gt) {
        GT()
        PREPROCESS()
    }
    else {
        PREPROCESS()
    }
    


    if (params.pan) {
        PREPREOCESS.out.fastagz
        PANGENOME()
        AUGMENT()
        ch_gfa = PANGENOME.out.gfa
        MODELPAN()
    }
    else if (params.asm) {
        PREPROCESS.out.gfau
        PREPROCESS.out.gfas
        ch_gfa = mix() 
        MODEL()
    }
    

    if (params.gt) {
        EVAL()
        PLASEVAL()
    }

}



process MAKE_PANGENOME {
    cache 'lenient'
    errorStrategy 'ignore'

    maxForks 4
    time '60m'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(pangenome), emit: pangenome

    script:
    haplos = 2
    pangenome = fasta.baseName + ".pan.gfa"

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
    panassembly = pangenome.baseName + "asm.gfa"
    pangenomecl = pangenome.baseName + "cl.gfa"
    """
    python $projectDir/bin/gfa_cleaner.py --input ${pangenome} --output ${pangenomecl}
    python $projectDir/bin/pangenome_enancher.py --input ${pangenomecl} --output ${panassembly} --assembly ${uni} ${ske}
    """
}
