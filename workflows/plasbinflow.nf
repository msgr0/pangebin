
putils_executable="${projectDir}/src/plasbin-flow/bin/plasbin_utils.py"
pflow_executable="${baseDir}/src/plasbin-flow/bin/plasbin_flow.py"
database="${baseDir}/src/plasbin-flow/database/plas_db.gfa"


process compress {

}


process make_input {

    input:
    tuple val(id), path(pangenome)

    output:
    tuple val(id), path(pangenome), path(input_csv)
    
    script:
    input_csv = "${id}-input.csv"
    """
    echo "sample,gfa" > "${input_csv}"
    echo "${id},${pangenome}" >> "${input_csv}"
    """

}

process gc_probs {
    input:
    tuple val(id), path(pangenome), path(input_csv)

    output:
    tuple val(id), path(pangenome), path(input_csv), path(gc_content)

    script:
    gc_content = "${id}.gc.tsv"
    """
    #!/bin/bash
    python ${putils_executable} gc_probabilities --input_file ${input_csv} --out_dir . --tmp_dir ./tmp
    """


}

process map_genes_to_ctgs {

}

process gene_density {

}

process model {

}

