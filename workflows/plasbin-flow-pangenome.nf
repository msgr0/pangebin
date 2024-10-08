// input: pan-assembly graph  + mlplasmid prediction/plasmid database

include { MLPLASMIDS } from "./mlplasmids"

binfolder = "~/bin"

putils = "$projectDir/PlasBin-flow-pangenome/code/plasbin_utils.py"
pflow = "$projectDir/PlasBin-flow-pangenome/code/plasbin_flow.py"
pdatabase = "$projectDir/Plasbin-flow-pangenome/database/genes.fasta"


#SBATCH --time=08:00:00
#SBATCH --account=def-chauvec
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=mattia.sgro@unimib.it
#SBATCH --mail-type=ALL
#SBATCH --mem=64G

echo 'Starting GPLAS PANGENOME pipeline'
module load nextflow/23.10.0
module load apptainer

for INPUT_FASTA in mixed_fastas_ecolxx/*.gz; do
	nextflow run nf-core/pangenome -r '1.1.2' -w $SLURM_TMPDIR/work -profile apptainer -resume --input $INPUT_FASTA --n_haplotypes 2 --outdir nf_output -params-file pangenome-params.json
done