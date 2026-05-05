#==================================================================================================
#   File: 07_populations.sh
#   Date: December 15, 2025
#   Description: Run Stacks' 'Populations' to generate vcf file
#--------------------------------------------------------------------------------------------------
#	Authors: Jared Homola
#==================================================================================================

#Define alias for project root directory
PROJ_ID=2505_mnBurbotStructure
YEAR=2025

MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID 
SCR=/mnt/gs21/scratch/homolaj1/$PROJ_ID 

module load Stacks/2.68-foss-2023a

cd $MEER_DIR/shell

echo '#!/bin/sh 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 3:59:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=24G
#SBATCH -J populations
#SBATCH -o '$MEER_DIR'/qstat/populations.o

populations -P '$MEER_DIR'/out -t 16 --vcf 

scontrol show job ${SLURM_JOB_ID} $' > populations.sh

sbatch populations.sh
