
process NCBI_FETCH {
    beforeScript 'sleep 2'
    maxForks 1 
    errorStrategy 'ignore'

    input:
    val(meta)
    
    output:
    tuple val(meta), path(genomic_fna_gz), emit: db

    script:
    genomic_fna_gz = "${meta.id}.genomic.fna.gz"
    """
    #!/bin/bash

    link=\$(python $projectDir/bin/get_assembly_link.py --input ${meta.id})
    rsync -vPp --chmod=u+w \$link ${genomic_fna_gz}
    """
}

process BLAST {
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(fasta_AB), path(genomic_fna_gz)

    output:
    tuple val(meta), path(ground_truth_csv), emit: gt
    tuple val(meta), path(blast_tsv)

    script:
    fna = "${meta.id}.fna"
    ground_truth_csv = "${meta.id}.gt.csv"
    blast_tsv = "${meta.id}.blast.tsv"

    """
    bgzip -d -c ${genomic_fna_gz} > ${fna}

    makeblastdb -in ${fna} -dbtype nucl
    blastn -db ${fna} -query ${fasta_AB} -out ${blast_tsv} -outfmt "6 qseqid sseqid pident length mismatch qstart qend bitscore"
    python3 $projectDir/bin/blast_utils.py --assembly ${fna} --csv ${blast_tsv} --fasta ${fasta_AB} --output ${ground_truth_csv}
    """

}

process EVAL {

    input:
    tuple val(meta), path(graph), path(gt_csv), path(bin_tsv)
    val type

    output:
    tuple val(meta), path("*.pdf"), path("*.stats.txt"), path("*.mix.csv"), emit: mix, optional: true

    script:
    """
    python3 $projectDir/bin/evaluate_bins.py --sample ${meta.id} --gfa ${graph} --csv ${gt_csv} --bin ${bin_tsv} --type ${type} --output .
    
    """
}

process MODE0 {
    input:
    tuple val(meta), path(mix_csv)
    val type

    output:
    tuple val(meta), path(eval_out), optional: true, emit: res
    tuple val(meta), path(eval_gt), path(eval_bin), optional: true, emit: mode1

    script:
    eval_out = "${meta.id}.${type}.eval.out"
    eval_gt = "${meta.id}.${type}.eval.gt.tsv"
    eval_bin = "${meta.id}.${type}.eval.bin.tsv"
    """
    python3 $projectDir/bin/convert_pred.py --input ${mix_csv} --output-gt ${eval_gt} --output-bin ${eval_bin}

    python3 $projectDir/PlasEval/src/plaseval.py eval --pred ${eval_bin} --gt ${eval_gt} --out_file ${eval_out}
    """

    
}

process MODE1 {
    input:
    tuple val(meta), path(gt), path(bins)
    val type

    output:
    tuple val(meta), path(eval_out), optional: true, emit: res
    
    script:
    eval_out = "${meta.id}.${type}.mod1.eval.out"
    eval_bins = "${meta.id}.${type}.mod1.bins.tsv"
    """
    python3 $projectDir/bin/modify_bins.py --mod 1 --bins ${bins} --gtruth ${gt} --output ${eval_bins}
    python3 $projectDir/PlasEval/src/plaseval.py eval --pred ${eval_bins} --gt ${gt} --out_file ${eval_out}

    """


}


workflow PLASEVAL {
    take:
    bins
    graph
    fasta_AB
    type
    

    main:
    eval_ch = Channel.empty()

    NCBI_FETCH(bins.map{meta, bins -> meta})
    BLAST(fasta_AB.join(NCBI_FETCH.out.db))
    EVAL(graph.join(BLAST.out.gt).join(bins), type)

    MODE0(EVAL.out.mix.map{meta, pdf, stats, mix -> tuple(meta, mix)}, type)
    if (type.matches("pangenome(.*)")) {
        MODE1(MODE0.out.mode1, type)
        out_mode1 = MODE1.out.res
    }
    else {
        out_mode1 = Channel.empty()
    }
    
    // MODE2(EVAL.out.mix.map{meta, pdf, stats, mix -> tuple(meta, mix)}, type)

    emit:
    mode0 = MODE0.out.res
    mode1 = out_mode1

}