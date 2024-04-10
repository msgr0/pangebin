ass_A = "short"
ass_B = "skesa"
ass_A_type = "uni"
ass_B_type = "ske"
src_dir = "$projectDir/bin"


process FASTA {
    input:
    tuple val(meta), path(A), path(B)

    output:
    tuple val(meta), path(ren_A), path(ren_B), emit: ren_asm
    tuple val(meta), path(fa_A), path(fa_B), emit: fa_asm
    tuple val(meta), path(fasta_AB), path(fasta_AB_gz), emit: pang

    script:
    dec_A = "${meta.id}.${ass_A}.gfa"
    dec_B = "${meta.id}.${ass_B}.gfa"

    ren_A = "${meta.id}.${ass_A}.ren.gfa"
    ren_B = "${meta.id}.${ass_B}.ren.gfa"
    
    fa_A = "${meta.id}.${ass_A}.fa"
    fa_B = "${meta.id}.${ass_B}.fa"

    fasta_AB = "${meta.id}.fasta"
    fasta_AB_gz = "${meta.id}.fasta.gz"


    """
    #!/bin/bash
    bgzip -d -c ${A} > ${dec_A}
    bgzip -d -c ${B} > ${dec_B}

    python $src_dir/rename_gfa.py --input ${dec_A} --output ${ren_A} --type ${ass_A_type}
    python $src_dir/rename_gfa.py --input ${dec_B} --output ${ren_B} --type ${ass_B_type}

    awk '/^S/{print ">"\$2; print ""\$3}' ${ren_A} > ${fa_A}
    awk '/^S/{print ">"\$2; print ""\$3}' ${ren_B} > ${fa_B}

    cat ${fa_A} ${fa_B} > ${fasta_AB}

    bgzip -c ${fasta_AB} > ${fasta_AB_gz}
    """
}

workflow PREPROCESS {

    take:
    input

    main:
    FASTA(input)

    emit:
    assemblies = FASTA.out.ren_asm
    fastas = FASTA.out.fa_asm
    pangenome = FASTA.out.pang
}
