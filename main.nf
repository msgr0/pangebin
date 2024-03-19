workflows = "$projectDir/workflows"
src_dir = "$projectDir/src"
test_sample = "ecol-SAMN04014855"

process ncbi_retrive {
    beforeScript 'sleep 2'
    maxForks 1 
    errorStrategy 'ignore'

    input:
    val(id)
    
    output:
    tuple val(id), path("${id}.genomic.fna.gz"), optional: true

    shell:
    """
    link=\$(python ${src_dir}/get_assembly_link.py --input ${id})
    rsync -vPp --chmod=u+w \$link ${id}.genomic.fna.gz
    """
}

workflow {
    sample_ids = Channel.of(test_sample)    
    sample_ids | view
    genomic_files = ncbi_retrive(sample_ids) | view
}