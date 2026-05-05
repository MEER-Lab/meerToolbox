#==================================================================================================
#   File: 06_genotype.sh
#   Date: December 15, 2025
#   Description: Genotype using GStacks
#--------------------------------------------------------------------------------------------------
#	Authors: Jared Homola
#==================================================================================================

#Define alias for project root directory
PROJ_ID=2505_mnBurbotStructure
YEAR=2025

MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID 
SCR=/mnt/gs21/scratch/homolaj1/$PROJ_ID 

module load Stacks/2.64-foss-2023a

mkdir $MEER_DIR/out

cd $MEER_DIR/shell

echo '#!/bin/sh 
#SBATCH --nodes=1-8
#SBATCH --ntasks=1
#SBATCH -t 2-0:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH -J genotype
#SBATCH -o '$MEER_DIR'/stdout/genotype.o

gstacks -I '$SCR'/mapped/filtered/sorted \
    -S .pe.filtered.sorted.bam \
    --threads 32 \
    -M '$MEER_DIR'/dependencies/pop.map \
    --rm-unpaired-reads \
	--details \
    -O '$MEER_DIR'/out

scontrol show job ${SLURM_JOB_ID} $' > genotype.sh

sbatch genotype.sh


