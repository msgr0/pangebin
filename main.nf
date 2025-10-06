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




workflow MAIN {
    take:
    input_ch

    main:

    meta_ch = input_ch.map { meta, _files -> meta }

    uni_ch = input_ch.map { meta, files -> meta += [asm: "u"]; def f = files[1] ; [meta, f]}

    ske_ch = input_ch.map { meta, files -> meta += [asm: "s"]; def f = files[0] ; [meta, f]}


    PREPROCESS(uni_ch.mix(ske_ch))



    gt_ch = PREPROCESS.out.panasmGfa.join(PREPROCESS.out.mixFasta).join( PREPROCESS.out.skeFasta ).join ( PREPROCESS.out.uniFasta)

    GT(meta_ch, gt_ch)
    
    GT.out.mixReference
    GT.out.skeReference
    GT.out.uniReference

    bin_ch = Channel.empty()
    res_ch = Channel.empty()
    naive_ch = Channel.empty()
    overlap_ch = Channel.empty()

    evalp_ch = Channel.empty()
    evals_ch = Channel.empty()
    evalu_ch = Channel.empty()

    if (params.pangenome) {
        PBFp(PREPROCESS.out.panasmGfa, PREPROCESS.out.gcPan, PREPROCESS.out.gdPan)


        bin_ch = bin_ch.mix(PBFp.out.bins)
        res_ch = res_ch.mix(PBFp.out.res)
        naive_ch = naive_ch.mix(PBFp.out.naive)
        overlap_ch = overlap_ch.mix(PBFp.out.overlap)
        
        
        evalpu_ch = res_ch.mix(naive_ch).mix(overlap_ch)
        evalpu_ch = evalpu_ch.combine(GT.out.uniReference.map{ meta, fragments, contigs -> [meta, fragments]}, by: 0)
        
        evalps_ch = res_ch.mix(naive_ch).mix(overlap_ch)
        evalps_ch = evalps_ch.combine(GT.out.skeReference.map { meta, fragments, contigs -> [meta, fragments]}, by: 0)        
        
        evalp_ch = evalpu_ch.mix(evalps_ch)
    }

    if (params.assemblers) {
            PBFu(PREPROCESS.out.uniGfa, PREPROCESS.out.gcUni, PREPROCESS.out.gdUni)
            PBFs(PREPROCESS.out.skesaGfa, PREPROCESS.out.gcSke, PREPROCESS.out.gdSke)
        
        // bin_ch = bin_ch.mix(PBFu.out.bins)
        // res_ch = res_ch.mix(PBFu.out.res) 
        evalu_ch = PBFu.out.res.combine(GT.out.uniReference.map{ meta, fragments, contigs -> [meta, contigs]}, by: 0)

        evalu_ch | view

        // bin_ch = bin_ch.mix(PBFs.out.bins)
        // res_ch = res_ch.mix(PBFs.out.res)  
        evals_ch = PBFs.out.res.combine(GT.out.skeReference.map{ meta, fragments, contigs -> [meta, contigs]}, by: 0)
        evals_ch | view

    }

    // pan_ch = PREPROCESS.out.panasmGfa.map{meta, files -> meta += [asm: "p"]; [meta, files]}

    // evalp_ch | view
    // evalu_ch | view
    // evals_ch | view 
    EVALUATION(evalp_ch.mix(evals_ch).mix(evalu_ch))
    

}


workflow { // single input workflow, for dataset input use pipeline.nf

    // input as in plasgraph
    // input_ch = Channel.fromFilePairs("${params.input}/*-S*-{s,u}.gfa.gz", type: 'file', checkIfExists: true)

    // standard input type for a sample folder SAMPLE/{skesa,unicycler}.gfa.gz
    input_ch = Channel.fromFilePairs("${params.input}/{skesa,unicycler}.gfa.gz", type: 'file', checkIfExists: true)


    metaid = file(params.input).toString().split("/")[-1]


    input_ch = input_ch.map { _m, filepair ->
        def fmeta = [:];
        fmeta.id = metaid;
        [fmeta, filepair]
    } | view


    MAIN(input_ch) 
}