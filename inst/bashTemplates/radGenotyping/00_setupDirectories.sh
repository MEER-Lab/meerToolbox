#==================================================================================================
#   File: 00_setupDirectories.sh
#   Date: May 5, 2026
#   Description: Set up directories, rename files to unifing, and move to scratch
#--------------------------------------------------------------------------------------------------
#   Author: Jared Homola
#==================================================================================================

# --- USER DEFINED VARIABLES ---
PROJ_ID=2505_mnBurbotStructure  ## Update to your MEER Lab project ID
YEAR=2025                       ## Update for the year
NETID=homolaj1                  ## Update to your MSU NetID

# --- PATH DEFINITIONS ---
MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID
SCR=/mnt/gs21/scratch/$NETID/$PROJ_ID         ## May need to update "gs21"

# --- DIRECTORY CREATION ---
echo "Initializing directories for project: $PROJ_ID"

mkdir -p $SCR/raw
mkdir -p $MEER_DIR/dependencies
mkdir -p $MEER_DIR/out
mkdir -p $MEER_DIR/stdout/{demultiplex,cloneFilter,bwaMemMap}
mkdir -p $MEER_DIR/shell/{demultiplex,cloneFilter,bwaMemMap}

echo "Setup complete. Paths established at:"
echo "MEER: $MEER_DIR"
echo "SCRATCH: $SCR"
