
putils_executable="${projectDir}/PlasBin-flow-pangenome/code/plasbin_utils.py"
pflow_executable="${projectDir}/PlasBin-flow-pangenome/code/plasbin_flow.py"
database="${projectDir}/PlasBin-flow-pangenome/database/plas_db.gfa"

include { COMPRESS } from "${projectDir}/workflows/utils"
include { MAKE_INPUT } from "${projectDir}/workflows/plasbinflow"
include { GENES2CTG } from "${projectDir}/workflows/plasbinflow"
include { GD_SCORE } from "${projectDir}/workflows/plasbinflow"


// process PLASMIDEC {

//     input:
//     tuple val(meta), path(graph)
    
//     output:
//     tuple val(meta), path(out_prediction)

//     script:
//     name = graph.getName().substring(0, graph.getName().length() - 4)
//     out_prediction = "${meta.id}/gplas2_format/${name}_plasmid_prediction.tab"
    
//     """
//     $projectDir/plasmidEC/plasmidEC.sh -i ${graph} -o ${meta.id} -g
//     """
// }


process GD_PREDICTION {
    input:
    tuple val(meta), path(graph), path(gd)

    output:
    tuple val(meta), path(out_prediction)

    script:
    out_prediction = "${meta.id}.gd.tab"
    """
    python $projectDir/bin/transform_gd_gplas.py --input ${gd} --output ${out_prediction} --gfa ${graph}
    """
}

process MODEL {
    errorStrategy 'ignore'
    
    conda "/home/sgro/micromamba/envs/gplas2"

    input:
    tuple val(meta), path(graph), path(prediction)

    output:
    tuple val(meta), path(res), emit: bins, optional: true

    script:
    res = "./results/${meta.id}_bins.tab"
    """
    gplas -i ${graph} -c extract -n '${meta.id}'
    gplas -c predict -i ${graph} -P ${prediction} -n '${meta.id}' -l 0
    """
}

workflow GPLAS {
    take:
    graph

    main:
    gd_csv = graph | COMPRESS | MAKE_INPUT | GENES2CTG | GD_SCORE 
    
    gd_pred = graph.join(gd_csv) | GD_PREDICTION

    MODEL(graph.join(gd_pred))

    emit:
    bins = MODEL.out.bins

}
