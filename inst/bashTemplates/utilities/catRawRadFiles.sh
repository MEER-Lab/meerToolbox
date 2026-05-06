#!/bin/bash
#==================================================================================================
#   File: 01b_concatenateReads.sh
#   Date: May 6, 2026
#   Description: SLURM submission to merge reads from multiple sequencing instances
#--------------------------------------------------------------------------------------------------
#   Author: Jared Homola
#==================================================================================================

# --- USER DEFINED VARIABLES ---
PROJ_ID=2528_MnRedhorse         ## Update to your MEER Lab project ID
YEAR=2025                       ## Update for the year
NETID=homolaj1                  ## Update to your MSU NetID

# Define source directories. Add or remove additional raw data locations as needed)
LOC1=/mnt/research/meerLab/rawData/2528_MnRedhorse-runB/20260429_KAN17711_RADSeq_PE150
LOC2=/mnt/research/meerLab/rawData/2528_MnRedhorse-runA/20260408_KAN17711_RADSeq_PE150

# --- PATH DEFINITIONS ---
MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID
SCR=/mnt/gs21/scratch/$NETID/$PROJ_ID 
OUT_DIR=$SCR/raw

# --- SETUP ---
mkdir -p "$SCR"
mkdir -p "$OUT_DIR"
mkdir -p "$MEER_DIR"
mkdir -p "$MEER_DIR/stdout"
mkdir -p "$MEER_DIR/shell"
cd "$MEER_DIR/shell"


# --- GUIDANCE FOR HPCC RESOURCE REQUESTS --- 
# Estimating time (-t): (Total GB / 10 GB per minute) x 1.5 (buffer).
## Calculate total GB as (add LOCs as needed): 
du -shc $LOC1 $LOC2 $LOC3

# Estimating memory (--mem): Under 100GB: 8G; 100–500GB: 16G; Over 500GB: 32G

# Estimating CPUs (--cpus-per-task): cat is a single thread processes, so only request 1 CPU


# --- GENERATE SLURM BATCH SCRIPT ---
cat <<EOF > run_concat_$PROJ_ID.sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -t 01:30:00
#SBATCH -J concat_$PROJ_ID
#SBATCH -o $MEER_DIR/stdout/concat_$PROJ_ID.o

# Improved sample discovery: Specifically look for 2528- prefix and strip extension
SAMPLES=\$(ls $LOC1/2528-*.1.fastq.gz $LOC1/2528-*_R1.fastq.gz 2>/dev/null | xargs -n 1 basename | sed -E 's/(_R1\.fastq\.gz|\.1\.fastq\.gz)//' | sort | uniq)

echo "Starting concatenation of samples..."

for SAMPLE in \$SAMPLES; do
    echo "Processing: \$SAMPLE"

    # Use wildcards to catch the files regardless of minor naming differences
    # This will merge matching files from LOC1 and LOC2 into your clean .1.fq.gz format
    cat $LOC1/\${SAMPLE}*R1.fastq.gz $LOC2/\${SAMPLE}*R1.fastq.gz > $OUT_DIR/\${SAMPLE}.1.fq.gz 2>/dev/null || true
    cat $LOC1/\${SAMPLE}*R2.fastq.gz $LOC2/\${SAMPLE}*R2.fastq.gz > $OUT_DIR/\${SAMPLE}.2.fq.gz 2>/dev/null || true
done

echo "Concatenation complete. Files located in: $OUT_DIR"
scontrol show job \${SLURM_JOB_ID}
EOF

# --- SUBMIT TO QUEUE ---
sbatch run_concat_$PROJ_ID.sh