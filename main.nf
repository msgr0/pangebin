#!/usr/bin/env nextflow

/*  PANGEBIN
----------------------------------------------
gtihub: https://github.com/msgr0/pangebin.git
----------------------------------------------
*/

include { GT                   } from './workflows/gt'
include { PREPROCESS           } from './workflows/preprocess'
// include { DOWNLOAD_GT          } from './workflows/gt'
include { MODEL as PBFp        } from './workflows/plasbin'
include { MODEL as PBFu        } from './workflows/plasbin'
include { MODEL as PBFs        } from './workflows/plasbin'

// include { EVAL                 } from './workflows/eval'
// include { PLASEVAL             } from './workflows/eval'


workflow {

    // input_ch = Channel.fromFilePairs("${params.input}/*-S*-{s,u}.gfa.gz", type: 'file', checkIfExists: true)
    
    input_ch = Channel.fromFilePairs("${params.input}/{unicycler,skesa}.gfa.gz", type: 'file', checkIfExists: true)

    metaid = file(params.input).toString().split("/")[-1]
    // if ... metaspecies ... ok 
    // metaspecies = file(params.input).toString().split("/")[-1].split("-")[0]
    // metaspecies = "generic"

    input_ch = input_ch.map { _m, filepair ->
        def fmeta = [:];
        fmeta.id = metaid;
        // fmeta.species = metaspecies;
        [fmeta, filepair]
    } | view

    meta_ch = input_ch.map { meta, _files -> meta }

    uni_ch = input_ch.map { meta, files -> meta += [asm: "u"]; def f = files[1] ; [meta, f]}

    ske_ch = input_ch.map { meta, files -> meta += [asm: "s"]; def f = files[0] ; [meta, f]}


    PREPROCESS(uni_ch.mix(ske_ch))

    gt_ch = PREPROCESS.out.panasmGfa.join(PREPROCESS.out.mixFasta).join( PREPROCESS.out.uniFasta ).join ( PREPROCESS.out.skeFasta)

    GT(meta_ch, gt_ch)

    bin_ch = Channel.empty()
    res_ch = Channel.empty()

    if (params.pangenome) {
        PBFp('pan', PREPROCESS.out.panasmGfa, PREPROCESS.out.gcPan, PREPROCESS.out.gdPan)

        bin_ch = bin_ch.mix(PBFp.out.bins)
        res_ch = res_ch.mix(PBFp.out.res)  
        // if (params.ml) {
            
        // }
    } 
    if (params.assembly) {
        PBFu("uni", PREPROCESS.out.uniGfa, PREPROCESS.out.gcUni, PREPROCESS.out.gdUni)
        bin_ch = bin_ch.mix(PBFu.out.bins)
        res_ch = res_ch.mix(PBFu.out.res) 
        
        PBFs("ske", PREPROCESS.out.skesaGfa, PREPROCESS.out.gcSke, PREPROCESS.out.gdSke)
        bin_ch = bin_ch.mix(PBFs.out.bins)
        res_ch = res_ch.mix(PBFs.out.res)  
    }


    
// output every bin modification for the pangenome
// run PLASEVAL
// run normal_EVAL
// include in the ID the bin result?


    // EVAL(PBFp.out.bins.mix(PBFa.out.bins), graphs, fastas, gt_everything)
    

}
