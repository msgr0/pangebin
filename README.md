# Pangebin Pipeline (WIP)

Pangebin is a pipeline aimed at performing the task of "plasmid binning" on a bacterial sample
exploiting a pangenome construction on a pair of assembly graph using `nf-core/pangenome`,
running a modified version of @chauvec's `Plasbin-Flow`.

**Requirements**
- nextflow + java 17
- python 3.10
- gurobi license
- conda/mamba

**Prior to run the pipeline**
Clone the repoistory, cd into it and checkout the `ismb2025` branch:
```
git clone ...
cd pangebin
git checkout ismb2025
```

Prepare your samples: provide for each sample a folder named after its biosample_id,
inside put the gfa(s) built for your short read sample using Unicycler and Skesa, named
respectively as `unicycler.gfa.gz` and `skesa.gfa.gz`.
Note that the `gfa` files are compressed using `bgzip` (from `htslib`).

The dataset used for the experiment used for ISMB2025 submission can be found in the `dataset.tar.gz` file.
De-compress it, the command will create a `database` folder.
```
tar xavf dataset.tar.gz
```
<!-- Moreover, you can find the whole dataset (complete of intermediate and output files) at ...
```
tar xvf dataset-whole.tar.gz --directory=./dataset-processed
``` -->
**Running the pipeline**


Run the pipeline on a sample of choice in the dataset (replace xxxxxx with an actual folder name found in `dataset`).
```
export SAMPLE_ID="xxxxxxxx"
```

Chose either `conda` or `mamba` as running backend, to run the main `pangebin` pipeline:
```
nextflow run . --input dataset/$SAMPLE_ID -profile conda
```
or 
```
nextflow run . --input dataset/$SAMPLE_ID -profile mamba
```

To run also original `plasbin-flow` on skesa and unicycler assemblies separately run:
```
nextflow run . --input dataset/$SAMPLE_ID -profile conda --assemblers
```
or
```
nextflow run . --input dataset/$SAMPLE_ID -profile mamba --assemblers
```
You can find the results in the same input folder:
- binning output
     - plasbin-flow `$SAMPLE_ID.{ske,uni}.1.pbf.bins.tsv`
     - pangebin `$SAMPLE_ID.pan.1.pbf.bins.tsv`
- binning results 
     - plasbin-flow `$SAMPLE_ID.{ske,uni}.1.pbf.pred.{ske,uni}.bin.txt`
     - pangebin: `$SAMPLE_ID.pan.1.pbf.pred.{uni,ske}.bin.txt`
- labeling results
     - plasbin-flow `$SAMPLE_ID.{ske,uni}.1.pbf.pred.{ske,uni}.lab.txt`
     - pangebin `$SAMPLE_ID.pan.1.pbf.pred.{uni,ske}.lab.txt`
