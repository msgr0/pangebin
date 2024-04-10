process COMPRESS {
    input:
    tuple val(id), path(file)

    output:
    tuple val(id), path (compressed)

    script:
    compressed = "${file}.gz"
    """
    #!/bin/bash
    bgzip -c ${file} > ${compressed}
    """
}
