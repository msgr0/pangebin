#!/bin/bash
# $1 is the sample folder with (uni.gfa.gz)(ske.gfa.gz)

nextflow workflows/preprocess.nf --input $1 --output $1/prep

# nextflow workflows/plasbin-flow-pangenome.nf --input $1/prep/panassembly.gfa --output $1/res

# nextflow workflows/gt.nf --input $1 --output $1/gt

# nextflow workflows/eval.nf --input $1(prep,res,gt) --output $1/out



