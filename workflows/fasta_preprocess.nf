process fastas_preprocess {
    input:
    tuple val(id), path(A), path(B)

    output:
    tuple val(id), path(ren_A), path(ren_B), path(fasta_AB_gz), path(fasta_AB)

    script:
    dec_A = "${id}.${ass_A}.gfa"
    dec_B = "${id}.${ass_B}.gfa"

    ren_A = "${id}.${ass_A}.ren.gfa"
    ren_B = "${id}.${ass_B}.ren.gfa"
    
    fa_A = "${id}.${ass_A}.fa"
    fa_B = "${id}.${ass_B}.fa"

    fasta_AB = "${id}.fasta"
    fasta_AB_gz = "${id}.fasta.gz"


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


workflow FASTA_PREPROCESS {

}
