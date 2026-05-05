#==================================================================================================
#   File: 03_removeClones.sh
#   Date: December 15, 2025
#   Description: Filter clonal reads from data
#--------------------------------------------------------------------------------------------------
#	Authors: Jared Homola
#==================================================================================================

#Define alias for project root directory
PROJ_ID=2505_mnBurbotStructure
YEAR=2025

MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID 
SCR=/mnt/gs21/scratch/homolaj1/$PROJ_ID 

mkdir $MEER_DIR/stdout/cloneFilter
mkdir $MEER_DIR/shell/cloneFilter
mkdir $SCR/cloneFilter

module load Stacks/2.64-foss-2023a

cd $MEER_DIR/shell/cloneFilter
ls $SCR/demultiplex/ | grep "fq" | grep -v ".2.fq.gz" | grep -v "rem" | grep -v "L2" | sed 's/\.1\.fq\.gz//g' | sort | uniq | while read -r LINE
do

echo '#!/bin/sh 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH -J cloneFilter.'$LINE'
#SBATCH -o '$MEER_DIR'/stdout/cloneFilter/cloneFilter.'$LINE'.o

cd '$SCR'/demultiplex
clone_filter -1 ./'$LINE'.1.fq.gz -2 ./'$LINE'.2.fq.gz -i gzfastq -o '$SCR'/cloneFilter/

scontrol show job ${SLURM_JOB_ID}' > ./cloneFilter."$LINE".sh

sbatch ./cloneFilter."$LINE".sh

done