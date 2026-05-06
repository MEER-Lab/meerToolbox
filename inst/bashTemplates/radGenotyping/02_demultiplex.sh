#==================================================================================================
#   File: 02_demultiplex.sh
#   Date: May 5, 2026
#   Description: Demultiplex libraries using process_radtags
#--------------------------------------------------------------------------------------------------
#   Author: Jared Homola
#==================================================================================================

# --- USER DEFINED VARIABLES ---
PROJ_ID=2505_mnBurbotStructure          ## Update to your MEER Lab project ID
YEAR=2025                               ## Update for the year
NETID=homolaj1                          ## Update to your MSU NetID

ENZYME=sbfI
TRIM_LEN=140
THREADS=8

# --- PATH DEFINITIONS ---
MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID
SCR=/mnt/gs21/scratch/$NETID/$PROJ_ID

# --- SETUP ---
mkdir -p $MEER_DIR/stdout/demultiplex
mkdir -p $MEER_DIR/shell/demultiplex
mkdir -p $MEER_DIR/dependencies

cd $MEER_DIR/shell/demultiplex

module load Stacks/2.64-foss-2023a

# --- GENERATE AND SUBMIT DEMULTIPLEX JOBS ---
# Libraries are inferred from *.1.fq.gz files in scratch/raw/
# Barcode files are expected as <library>_barcodes.txt

ls '$SCR'/raw | grep "\.1\.fq\.gz" | awk -F ".1.fq.gz" '{print $1}' | uniq | sort | while read -r LINE
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
  -e '$ENZYME' \
  -c \
  -q \
  -r \
  --bestrad \
  -t '$TRIM_LEN' \
  --threads '$THREADS' \
  -o '$SCR'/demultiplex
' > demultiplex.$LINE.sh

  sbatch demultiplex.$LINE.sh

done

module reset
