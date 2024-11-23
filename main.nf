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
include { DOWNLOAD_GT          } from './workflows/gt'
include { MODEL as PBFp        } from './workflows/plasbin'
include { MODEL as PBFa        } from './workflows/plasbin'
// include { EVAL                 } from './workflows/eval'
// include { PLASEVAL             } from './workflows/eval'


workflow {

    // input_ch = Channel.fromFilePairs("${params.input}/*-S*-{s,u}.gfa.gz", type: 'file', checkIfExists: true)
    
    input_ch = Channel.fromFilePairs("${params.input}/{short,skesa}.gfa.gz", type: 'file', checkIfExists: true)

    metaid = file(params.input).toString().split("-")[1]
    metaspecies = file(params.input).toString().split("/")[-1].split("-")[0]

    input_ch = input_ch.map { _m, filepair ->
        def fmeta = [:];
        fmeta.id = metaid;
        fmeta.species = metaspecies;
        [fmeta, filepair]
    } | view

    meta_ch = input_ch.map { meta, _files -> meta }

    uni_ch = input_ch.map { meta, files -> meta += [asm: "u"]; f1 = files[0]  ; [meta, f1]}

    ske_ch = input_ch.map { meta, files -> meta += [asm: "s"]; f2 = files[1] ; [meta, f2]}


    PREPROCESS(uni_ch.mix(ske_ch))

    gt_ch = PREPROCESS.out.gfa_psm.join(PREPROCESS.out.fasta_mix).join( PREPROCESS.out.fasta_uni ).join ( PREPROCESS.out.fasta_ske)

    GT(meta_ch, gt_ch)

    PBFp('pan', PREPROCESS.out.gfa_psm, PREPROCESS.out.fasta_mix)
    // PBFp('asm', PREPROCESS.out.gfa_uni.mix(PREPROCESS.out.gfa_ske), PREPROCESS.out.fasta_uni.mix(PREPROCESS.out.fasta_ske))

    // EVAL(PBFp.out.bins.mix(PBFa.out.bins), graphs, fastas, gt_everything)
    

}
