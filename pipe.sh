#!/bin/bash

# preprocess folder/(uni.gfa.gz)(ske.gfa.gz) -> folder/out/()
nextflow workflows/preprocess.nf --input $1 --output $1/out
nextflow workflos/plasbin-flow-pangenome.nf --input $1/out/panassembly.gfa

