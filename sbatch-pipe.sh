#SBATCH --time=01:00:00
#SBATCH --account=def-chauvec
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=mattia.sgro@unimib.it
#SBATCH --mail-type=ALL
#SBATCH --mem=64G


module 

bash ./pipe.sh ~/scratch/pangebin/test_data/ecol-SAMN04014855