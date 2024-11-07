# Pangebin Pipeline (WIP)

Pangebin is a collection of methods aimed at performing the task of "plasmid binning" on a bacterial sample,
exploiting a pangenome construction on a pair of assembly graph, running a modified version of @chauvec 's Plasbin-Flow
and, lastly, 

**Requirements**
- nextflow
- PBF requirements:

using conda
using docker/podman, a single docker image prebuilt for this task
using singularity, a single singularity image prebuilt for this task
 
using it locally (we assume that every requirement is satisfied locally)


**Running the pipeline**
1. pangenome construction   `--profile step1`\
     `input`: `u.gfa.gz`, `s.gfa.gz`\
     `output`: `psm.gfa.gz`\
     
2. model + binning          `--profile step2`\
    `input`:  `psm.gfa.gz`\
    `output`:  `bins.tsv`\

3. evaluation               `--profile step3`\
     `input`:  `psm.gfa.gz`, `bins.tsv`, ID\
     `output`:  `eval.labeling.txt`, `eval.binning.txt`\
     
- whole pipeline             `--profile all` (or no profile)\
     `input`:  `u.gfa.gz`, `s.gfa.gz`, ID\
     `output`:  `pangenome`, `bins`, `eval.labeling.txt`, `eval.binning.txt`\

## Assembly graphs and Pangenome construction

The input is composed by:
 - unicycler and skesa graphs in .gfa or .gfagz format
 - pangenome graph built with nf-core Pangenome (PGGB).


## Plasbin-Flow + pangenome
tbd

## Output Bins
tbd


