#==================================================================================================
#   File: GTscore_runHPCC.sh
#   Date: July 8, 2025
#   Description: Run GTscore on the MSU HPCC
#--------------------------------------------------------------------------------------------------
#	Authors: Jared Homola
#==================================================================================================

#Define alias for project root directory
PROJ_ID=2523_waeBroodstock20242025
MEER_DIR=/mnt/research/meerLab/2025/$PROJ_ID 

mkdir $MEER_DIR/shell
cd $MEER_DIR/shell

# Load necessary modules
module purge
module load R/4.4.1-gfbf-2023b

echo '#!/bin/sh 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 3:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH -J GTscore
#SBATCH -o '$MEER_DIR'/stdout/GTscore.o

# Set paths
PROJECT_DIR="/mnt/research/meerLab/2025/2523_waeBroodstock20242025/genotyping"
GTSCORE_PATH="/mnt/research/meerLab/software/GTscore/analysis/GTscore.R"
OUTPUT_PREFIX="2523_"

# Run the R script
Rscript /mnt/research/meerLab/software/GTscore/GTscore_HPCC.R "$PROJECT_DIR" "$GTSCORE_PATH" "$OUTPUT_PREFIX"

scontrol show job ${SLURM_JOB_ID} $' > GTscore_runHPCC.sh

sbatch GTscore_runHPCC.sh
