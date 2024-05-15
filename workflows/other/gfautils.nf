nameA = "short"
nameB = "skesa"
typeA = "uni"
typeB = "ske"
src_dir = "$projectDir/bin"

process PROCESS_INPUT_GFA_GPLAS {
    conda "/opt/mambaforge/envs/gplas2"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path(gfa), emit: gfa
    tuple val(meta), path(fa), emit: fasta

    script:
    fa = "${meta.id}.${typeA}.gplas.fa"
    """
    gplas -i ${gfa} -c extract -n '${fa}'
    """

}

process PROCESS_INPUT_GFA {
    conda "bioconda::htslib"

    input:
    tuple val(meta), path(gfaA), path(gfaB) 

    output:
    tuple val(meta), path(renA), path(renB), emit: gfa
    tuple val(meta), path(faA), path(faB), emit: fasta
    tuple val(meta), path(faAB), path(faABgz), emit: fastapangenome

    script:
    renA = "${gfaA.getName()[0..-3]}.ren.gfa"
    renB = "${gfaB.getName()[0..-3]}.ren.gfa"
    
    renA= "${meta.id}.${typeA}.ren.gfa"
    renB = "${meta.id}.${typeB}.ren.gfa"
    faA = "${meta.id}.${typeA}.fa"
    faB = "${meta.id}.${typeB}.fa"

    faAB = "${meta.id}.mixed.fasta"
    faABgz = "${meta.id}.mixed.fasta.gz"

    """
    #!/bin/bash
    python $src_dir/rename_gfa.py --input ${gfaA} --output ${renA} --type ${typeA}
    python $src_dir/rename_gfa.py --input ${gfaB} --output ${renB} --type ${typeB}

    awk '/^S/{print ">"\$2; print ""\$3}' ${renA} > ${faA}
    awk '/^S/{print ">"\$2; print ""\$3}' ${renB} > ${faB}

    cat ${faA} ${faB} > ${faAB}

    bgzip -c ${faAB} > ${faABgz}
    """
}

process PACK_GZ {
    cache 'lenient'
    conda 'bioconda::htslib'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(compressed)

    script:
    compressed = "${file}.gz"
    """
    #!/bin/bash
    bgzip -c ${file} > ${compressed}
    """
}

process UNPACK_GZ {
    cache 'lenient'
    conda 'bioconda::htslib'

    input:
    path(compressed)

    output:
    path(file_)

    script:
    file_ = "${compressed.getName()[0..-2]}"
    """
    #!/bin/bash
    bgzip -d -k ${compressed}
    """
}


process GFA_TO_FASTA {
    /// extract fasta using GPLAS TOOL for compatibility acroos diverse tools
    cache 'lenient'
    conda "/opt/mambaforge/envs/gplas2"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path(fasta)

    script:
    """
    gplas -i ${gfa} -c extract -n '${meta.id}'
    """
}