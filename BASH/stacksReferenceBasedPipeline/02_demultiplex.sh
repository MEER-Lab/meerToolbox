#==================================================================================================
#   File: 02_demultiplex.sh
#   Date: December 15, 2025
#   Description: Demultiplex libraries using process_radtags
#--------------------------------------------------------------------------------------------------
#   Author: Jared Homola
#==================================================================================================

#Define alias for project root directory
PROJ_ID=2505_mnBurbotStructure
YEAR=2025

MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID 
SCR=/mnt/gs21/scratch/homolaj1/$PROJ_ID 

mkdir $MEER_DIR/stdout/demultiplex
mkdir $MEER_DIR/shell/demultiplex
mkdir $SCR/demultiplex

cd $MEER_DIR/shell/demultiplex

module load Stacks/2.64-foss-2023a

ls $SCR/raw | grep -v ".2.fq" | awk -F ".1.fq" '{print $1}' | uniq | sort | while read -r LINE
do
echo '#!/bin/sh 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G
#SBATCH -J demultiplex_'$LINE'
#SBATCH -o '$MEER_DIR'/stdout/demultiplex/demultiplex_'$LINE'.o

process_radtags \
    -1 '$SCR'/raw/'$LINE'.1.fq.gz \
    -2 '$SCR'/raw/'$LINE'.2.fq.gz \
    -i gzfastq \
    -y gzfastq \
    -b '$MEER_DIR'/dependencies/'$LINE'_barcodes.txt \
    -e sbfI \
    -c \
    -q \
    -r \
    --bestrad \
    -t 140 \
    --threads 8 \
    -o '$SCR'/demultiplex/ 

scontrol show job ${SLURM_JOB_ID}' > demultiplex.$LINE.sh

sbatch demultiplex.$LINE.sh
done
