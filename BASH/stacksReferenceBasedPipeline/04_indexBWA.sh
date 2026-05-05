#==================================================================================================
#   File: 04_indexBWA.sh
#   Date: December 15, 2025
#   Description: Perform BWA indexing
#--------------------------------------------------------------------------------------------------
#	Authors: Jared Homola
#==================================================================================================

#Define alias for project root directory
PROJ_ID=2505_mnBurbotStructure
YEAR=2025

MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID 
SCR=/mnt/gs21/scratch/homolaj1/$PROJ_ID 

cd $MEER_DIR/referenceGenome 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/801/175/GCA_016801175.1_ASM1680117v1/GCA_016801175.1_ASM1680117v1_genomic.fna.gz
gunzip GCA_016801175.1_ASM1680117v1_genomic.fna.gz

module load BWA/0.7.17-GCCcore-12.3.0

cd $MEER_DIR/shell

echo '#!/bin/sh 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 3:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=72G 
#SBATCH -J bwaIndex
#SBATCH -o '$MEER_DIR'/stdout/bwaIndex.o

bwa index '$MEER_DIR'/referenceGenome/GCA_026745465.1_ASM2674546v1_genomic.fna' > bwaIndex.sh

sbatch bwaIndex.sh


