
putils_executable="${projectDir}/PlasBin-flow-pangenome/code/plasbin_utils.py"
pflow_executable="${projectDir}/PlasBin-flow-pangenome/code/plasbin_flow.py"
database="${projectDir}/PlasBin-flow-pangenome/database/plas_db.gfa"

include { PACK_GZ as COMPRESS } from "${projectDir}/workflows/other/gfautils"


process MAKE_INPUT {

    input:
    tuple val(meta), path(graph)

    output:
    tuple val(meta), path(graph), path(input_csv)
    
    script:
    input_csv = "${meta.id}-input.csv"
    """
    #!/bin/bash
    echo "sample,gfa" > "${input_csv}"
    echo "${meta.id},${graph}" >> "${input_csv}"
    """

}

process GC_PROBS {
    conda 'matplotlib pandas seaborn biopython scipy numpy'

    input:
    tuple val(meta), path(graph), path(input_csv)

    output:
    tuple val(meta), path(graph), path(input_csv), path(gc_content)

    script:
    gc_content = "${meta.id}.gc.tsv"
    """
    #!/bin/bash
    python ${putils_executable} gc_probabilities --input_file ${input_csv} --out_dir . --tmp_dir ./tmp
    """

}

process GENES2CTG {
    conda 'matplotlib pandas seaborn biopython scipy numpy blast'

    input:
    tuple val(meta), path(graph), path(input_csv)

    output:
    tuple val(meta), path(graph), path(input_csv), path(mapping), path("*.genes_mappings.txt")

    script:
    mapping = "${meta.id}.mapping.csv"
    """
    #!/bin/bash
    python ${putils_executable} map_genes_to_ctgs --input_file ${input_csv} --out_dir . --tmp_dir ./tmp --db_file ${database} --out_file ${mapping}

    """

}

process GD_SCORE {
    conda 'matplotlib pandas seaborn biopython scipy numpy'

    input:
    tuple val(meta), path(graph), path(input_csv), path(mapping_csv), path(mapping_txt)

    output:
    tuple val(meta), path("*.gd.tsv")

    script:
    """
    python ${putils_executable} gene_density --input_file ${mapping_csv}  --out_dir . --tmp_dir ./tmp
    """
}

process MILP {
    conda 'matplotlib pandas seaborn biopython scipy numpy gurobi::gurobi pip:networkx'

    maxForks 1

    input:
    tuple val(meta), path(graph), path(input_csv), path(gc_tsv), path(gd_tsv)
    val asm
    val seed_len
    val seed_score
    val min_pls_len
    output:
    tuple val(meta), path(bins), optional: true

    script:
    bins = "${meta.id}.${asm}.bins.tsv"
    if (asm == "pangenome_no_penalty") {
        asm = "pangenome"
    }
    if (asm != "pangenome") {
        penalty = "--nopenalty"
    } else {
        penalty = ""
    }

    """
    python ${pflow_executable} ${penalty} -ag ${graph} -gc ${gc_tsv} -out_dir . -out_file ${bins} -alpha4 1 -score ${gd_tsv} -assembler ${asm} -seed_len ${seed_len}  -seed_score ${seed_score} -min_pls_len ${min_pls_len}
    """

}

workflow PLASBIN {

    take:
    graph
    asm
    seed_len
    seed_score
    min_pls_len

    main:
    input_csv = graph | COMPRESS | MAKE_INPUT
    gc_csv = input_csv | GC_PROBS
    mapping_csv = input_csv | GENES2CTG
    gd_csv = mapping_csv | GD_SCORE
    model_results = MILP(gc_csv.join(gd_csv), asm, seed_len, seed_score, min_pls_len)


    emit:
    bins = model_results

}
