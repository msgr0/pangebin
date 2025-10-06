//
// Check input FASTA and prepare indices
//

include { TABIX_BGZIP                 } from '../../modules/nf-core/tabix/bgzip/main.nf'
include { SAMTOOLS_FAIDX              } from '../../modules/nf-core/samtools/faidx/main.nf'

workflow INPUT_CHECK {
    take:
    input_ch // file: /path/to/sequences.fasta

    main:

    ch_versions = Channel.empty() // we collect all versions here
    ch_fasta = Channel.empty() // final output channel [ val(meta) , [ fasta ] ]
    fai = Channel.empty() // we store the .fai index here [ fai ]
    gzi = Channel.empty() // we store the .gzi index here [ gzi ]


    SAMTOOLS_FAIDX(input_ch, [[],[]])
    fai = SAMTOOLS_FAIDX.out.fai
    gzi = SAMTOOLS_FAIDX.out.gzi
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_fasta = input_ch


    emit:
    fasta = ch_fasta         // channel: [ val(meta), [ fasta ] ]
    fai = fai                // channel: [ val(meta), fasta.fai ]
    gzi = gzi                // channel: [ val(meta), fasta.gzi ]
    versions = ch_versions   // channel: [ versions.yml ]
}
