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



workflow BLAST {

}