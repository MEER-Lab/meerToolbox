#==================================================================================================
#   File: 01_rename.sh
#   Date: May 5, 2026
#   Description: Standardize raw file names and move to scratch
#--------------------------------------------------------------------------------------------------
#   Author: Jared Homola
#==================================================================================================


##### Untested as of 5/6/2026 #####


# --- USER DEFINED VARIABLES ---
PROJ_ID=2505_mnBurbotStructure          ## Update to your MEER Lab project ID
YEAR=2025                               ## Update for the year
NETID=homolaj1                          ## Update to your MSU NetID
RAW_DATA_SUBDIR="20250603_DNASeq_PE150" ## The folder name provided by the sequencing center

# --- PATH DEFINITIONS ---
MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID
SCR=/mnt/gs21/scratch/$NETID/$PROJ_ID
INPUT_DIR=$MEER_DIR/$RAW_DATA_SUBDIR

# --- AUTOMATED RENAMING LOOP ---
echo "Standardizing file names for project $PROJ_ID..."

# This looks for any R1 fastq file in your input directory
# It assumes a pattern of: projectID-LibraryName_R1.fastq.gz
ls $INPUT_DIR/*_R1.fastq.gz | while read -r R1_PATH; do
    
    # Extract the library number (e.g., P1, P2) from the filename
    LIB_NUM=$(basename "$R1_PATH" | cut -d'-' -f4 | cut -d'_' -f1 | sed 's/P//')
    
    # Define the R2 path by swapping R1 for R2
    R2_PATH=${R1_PATH/_R1/_R2}
    
    echo "Processing Library $LIB_NUM..."
    
    # Copy and rename to your scheme: project_library.1or2.fq.gz
    cp "$R1_PATH" "$SCR/raw/${PROJ_ID:0:4}_${LIB_NUM}.1.fq.gz"
    cp "$R2_PATH" "$SCR/raw/${PROJ_ID:0:4}_${LIB_NUM}.2.fq.gz"
done

echo "Renaming complete. Files are now in $SCR/raw/"