#==================================================================================================
#   File: 00_setupDirectories.sh
#   Date: December 15, 2025
#   Description: Set up directories, rename files to unifing, and move to scratch
#--------------------------------------------------------------------------------------------------
#   Author: Jared Homola
#==================================================================================================
#Define alias for project root directory
PROJ_ID=2505_mnBurbotStructure
YEAR=2025

MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID 
SCR=/mnt/gs21/scratch/homolaj1/$PROJ_ID 

mkdir $SCR
mkdir $SCR/raw
mkdir $MEER_DIR/raw  #### REMOVE??
mkdir $MEER_DIR/stdout
mkdir $MEER_DIR/shell
mkdir $MEER_DIR/out
mkdir $MEER_DIR/dependencies
mkdir $MEER_DIR/referenceGenome 
