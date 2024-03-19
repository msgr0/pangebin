workflows = "$projectDir/workflows"
src_dir = "$projectDir/bin"
test_sample = "ecol-SAMN04014855"
test_dataset = "$projectDir/test_data"
ass_A = "short"
ass_B = "skesa"
ass_A_type = "uni"
ass_B_type = "ske"

process ncbi_retrive {
    beforeScript 'sleep 2'
    maxForks 1 
    errorStrategy 'ignore'

    input:
    val(id)
    
    output:
    tuple val(id), path(genomic_fna_gz), optional: true

    script:
    genomic_fna_gz = "${id}.genomic.fna.gz"
    """
    #!/bin/bash

    link=\$(python ${src_dir}/get_assembly_link.py --input ${id})
    rsync -vPp --chmod=u+w \$link ${genomic_fna_gz}
    """
}

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

process blast {
    errorStrategy 'ignore'

    input:
    tuple val(id), path(fasta_AB), path(genomic_fna_gz)

    output:
    tuple val(id), path(ground_truth_csv), path(blast_tsv)

    script:
    fna = "${id}.fna"
    ground_truth_csv = "${id}.gt.csv"
    blast_tsv = "${id}.blast.tsv"

    """
    bgzip -d -c ${genomic_fna_gz} > ${fna}

    makeblastdb -in ${fna} -dbtype nucl
    blastn -db ${fna} -query ${fasta_AB} -out ${blast_tsv} -outfmt "6 qseqid sseqid pident length mismatch qstart qend bitscore"
    python3 $src_dir/blast_utils.py --assembly ${fna} --csv ${blast_tsv} --fasta ${fasta_AB} --output ${ground_truth_csv}
    """

}

process make_pangenome {
    maxForks 1

    input:
    tuple val(id), path(fasta_AB_gz)

    output:
    tuple val(id), path(pangenome)

    script:
    pangenome = "${id}.pangenome.gfa"
    haplos = 2
    paramfile = "${src_dir}/param-pangenome.yaml"
    """
    nextflow run nf-core/pangenome -r 1.1.0 -profile conda --input ${fasta_AB_gz} --n_haplotypes ${haplos} --outdir ${id}-nfcore -params-file ${paramfile}
    mv ${id}-nfcore/FINAL_GFA/${fasta_AB_gz}.gfaffix.unchop.Ygs.view.gfa ${id}.pangenome.gfa
    """
    
}

workflow {
    // sample_ids = Channel.of(test_sample)
    sample_folders = Channel.fromPath("${test_dataset}/*", type: 'dir')
    sample_ids = sample_folders.map{it-> it.getName()}
    genomic_files = ncbi_retrive(sample_ids)

    fastas = sample_folders.map{dir ->  [dir.getName(), file("$dir/${ass_A}.gfa.gz"), file("$dir/${ass_B}.gfa.gz")]}
    fastas_prep = fastas_preprocess(fastas)
    blast_info = blast(fastas_prep.map{it -> [it[0], it[-1]]}.join(genomic_files))
    pangenome = make_pangenome(fastas_prep.map{it->[it[0], it[-2]]})

    
     
}