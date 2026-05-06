#==================================================================================================
#   File: 02_demultiplex.sh
#   Date: May 5, 2026
#   Description: Demultiplex libraries using process_radtags via SLURM
#--------------------------------------------------------------------------------------------------
#   Author: Jared Homola
#==================================================================================================

# Note: barcode files are currently made manually by saving each tab of the Excel sheet as a .txt file. Future iterations may automate this step if the format is consistent.

# --- USER DEFINED VARIABLES ---
PROJ_ID=2528_MnRedhorse                 ## Update to your MEER Lab project ID
YEAR=2025                               ## Update for the year
NETID=homolaj1                          ## Update to your MSU NetID

ENZYME=sbfI
TRIM_LEN=140
THREADS=16

# --- PATH DEFINITIONS ---
MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID
SCR=/mnt/gs21/scratch/$NETID/$PROJ_ID

# --- SETUP ---
mkdir -p "$MEER_DIR/stdout/demultiplex" "$MEER_DIR/shell/demultiplex" "$SCR/demultiplex"
cd "$MEER_DIR/shell/demultiplex"

module purge
module load Stacks/2.64-foss-2023a

# --- GENERATE AND SUBMIT JOBS ---
ls "$SCR/raw" | grep "\.1\.fq\.gz" | sed 's/\.1\.fq\.gz//' | sort | uniq | while read -r LINE
do
  BARCODE_FILE="$MEER_DIR/dependencies/${LINE}_barcodes.txt"

  # Check for barcodes before submitting to avoid immediate job failure
  if [ ! -f "$BARCODE_FILE" ]; then
    echo "Warning: Barcode file missing for $LINE, skipping..."
    continue
  fi

  cat <<EOF > demultiplex."$LINE".sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem-per-cpu=12G
#SBATCH -J demultiplex_$LINE
#SBATCH -o $MEER_DIR/stdout/demultiplex/demultiplex_$LINE.o

module load Stacks/2.64-foss-2023a

process_radtags \\
  -1 $SCR/raw/$LINE.1.fq.gz \\
  -2 $SCR/raw/$LINE.2.fq.gz \\
  -i gzfastq \\
  -y gzfastq \\
  -b $BARCODE_FILE \\
  -e $ENZYME \\
  -c -q -r \\
  --bestrad \\
  -t $TRIM_LEN \\
  --threads $THREADS \\
  -o $SCR/demultiplex

scontrol show job \${SLURM_JOB_ID}
EOF

  sbatch demultiplex."$LINE".sh
done

module reset