# Pangebin (WIP)
Pangebin is a nextflow pipeline that permites to run plasmid binning software (Plasbinflow, Gplas for example) on multiple assembly graphs at once, exploiting pangenome graphs capabilities.

Pangebin pipeline takes as input two assembly graphs of a bacterial sample,
built from Skesa and Unicycler, builds the augmented pangenome that could be then used as infput for PlasBin-Flow or Gplas.

## Requirements
Software required for running the pipeline
- conda/mamba
- nextflow (can be installed via conda -c bioconda)

## Input
The input graphs should be named `short.gfa.gz` for unicycler graphs and `skesa.gfa.gz` for skesa graphs. Both should be placed into a folder named after the `sample_id`. 

## Running the pipeline
`nextflow run main.nf -profile mamba --db folder_sample_id --out output`
will run the pipeline, placing the output gfa into the `output` folder.

# still work-in-progress

