#!/usr/bin/env nextflow

/*  PANGEBIN
----------------------------------------------
gtihub: https://github.com/msgr0/pangebin.git
----------------------------------------------
*/

include { GT                   } from './workflows/gt'
include { PREPROCESS           } from './workflows/preprocess'
include { MODEL as PBFp        } from './workflows/plasbin'
include { MODEL as PBFu        } from './workflows/plasbin'
include { MODEL as PBFs        } from './workflows/plasbin'
include { EVALUATION           } from './workflows/evaluation'


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
    }

    meta_ch = input_ch.map { meta, _files -> meta }

    uni_ch = input_ch.map { meta, files -> meta += [asm: "u"]; def f = files[1] ; [meta, f]}

    ske_ch = input_ch.map { meta, files -> meta += [asm: "s"]; def f = files[0] ; [meta, f]}


    PREPROCESS(uni_ch.mix(ske_ch))



    gt_ch = PREPROCESS.out.panasmGfa.join(PREPROCESS.out.mixFasta).join( PREPROCESS.out.uniFasta ).join ( PREPROCESS.out.skeFasta)

    GT(meta_ch, gt_ch)

    bin_ch = Channel.empty()
    res_ch = Channel.empty()
    naive_ch = Channel.empty()
    overlap_ch = Channel.empty()

    evalp_ch = Channel.empty()
    evals_ch = Channel.empty()
    evalu_ch = Channel.empty()

    if (params.pangenome) {
        PBFp('pan', PREPROCESS.out.panasmGfa, PREPROCESS.out.gcPan, PREPROCESS.out.gdPan)

        bin_ch = bin_ch.mix(PBFp.out.bins)
        res_ch = res_ch.mix(PBFp.out.res)
        naive_ch = naive_ch.mix(PBFp.out.naive)
        overlap_ch = overlap_ch.mix(PBFp.out.overlap)
        
        
        evalpu_ch = res_ch.mix(naive_ch).mix(overlap_ch) | view
        evalpu_ch = evalpu_ch.join(GT.out.uniReference.map{ meta, fragments, contigs -> [meta, fragments]})
        
        evalps_ch = res_ch.mix(naive_ch).mix(overlap_ch) | view
        evalps_ch = evalps_ch.join(GT.out.skeReference.map { meta, fragments, contigs -> [meta, fragments]})        
        
        evalp_ch = evalpu_ch.mix(evalps_ch)
    }

    if (params.assembly) {
        PBFu("uni", PREPROCESS.out.uniGfa, PREPROCESS.out.gcUni, PREPROCESS.out.gdUni)
        bin_ch = bin_ch.mix(PBFu.out.bins)
        res_ch = res_ch.mix(PBFu.out.res) 
        evalu_ch = res_ch.join(GT.out.uniReference.map{ meta, fragments, contigs -> [meta, contigs]})

        PBFs("ske", PREPROCESS.out.skesaGfa, PREPROCESS.out.gcSke, PREPROCESS.out.gdSke)
        bin_ch = bin_ch.mix(PBFs.out.bins)
        res_ch = res_ch.mix(PBFs.out.res)  
        evals_ch = res_ch.join(GT.out.skeReference.map{ meta, fragments, contigs -> [meta, contigs]})

    }

    pan_ch = PREPROCESS.out.panasmGfa.map{meta, files -> meta += [asm: "p"]; [meta, files]}

    
    EVALUATION(evalp_ch.mix(evalu_ch).mix(evals_ch))
    

}
