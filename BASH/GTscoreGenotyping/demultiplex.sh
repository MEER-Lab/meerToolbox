#==================================================================================================
#   File: demultiplex.sh
#   Date: July 6, 2025
#   Description: Demultiplex using fastq-multx
#--------------------------------------------------------------------------------------------------
#	Authors: Jared Homola
#==================================================================================================

salloc -t 48:00:00 --mem=64G


#Define alias for project root directory
PROJ_ID=/mnt/research/meerLab/2025/2523_waeBroodstock20242025

cd $PROJ_ID

# Load conda module
module purge
module load Miniforge3

# Activate your conda environment
conda activate ea_utils_env

# Change to your working directory
cd $PROJ_ID
mkdir demux

# Run your command
fastq-multx \
-B barcodes25092523.txt \
/mnt/research/meerLab/rawData/NEW_2523-2509_Download/20251024_DNASeq_PE150/WAE_2523-2509_I1.fastq.gz \
/mnt/research/meerLab/rawData/NEW_2523-2509_Download/20251024_DNASeq_PE150/WAE_2523-2509_I2.fastq.gz \
/mnt/research/meerLab/rawData/NEW_2523-2509_Download/20251024_DNASeq_PE150/WAE_2523-2509_R1.fastq.gz \
/mnt/research/meerLab/rawData/NEW_2523-2509_Download/20251024_DNASeq_PE150/WAE_2523-2509_R2.fastq.gz \
-o demux/%_I1.fastq.gz \
-o demux/%_I2.fastq.gz \
-o demux/%_R1.fastq.gz \
-o demux/%_R2.fastq.gz \
-m 0 \
-d 1 

zcat /mnt/research/meerLab/rawData/2523_waeBroodstock20242025_2509_upWaePbt/20251024_DNASeq_PE150/WAE_2523-2509_R1.fastq.gz | wc -l
zcat /mnt/research/meerLab/rawData/2523_waeBroodstock20242025_2509_upWaePbt/20251024_DNASeq_PE150/WAE_2523-2509_R2.fastq.gz | wc -l