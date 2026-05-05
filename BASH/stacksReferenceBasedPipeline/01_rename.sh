#==================================================================================================
#   File: 01_rename.sh
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


### Name files manually to a scheme of "projectNumber_libraryNumber.1or2.fq.gz"
cp $MEER_DIR/20250603_DNASeq_PE150/2505-Ll-RAD-P1_R1.fastq.gz $SCR/raw/2505_1.1.fq.gz
cp $MEER_DIR/20250603_DNASeq_PE150/2505-Ll-RAD-P1_R2.fastq.gz $SCR/raw/2505_1.2.fq.gz

cp $MEER_DIR/20250603_DNASeq_PE150/2505-Ll-RAD-P2_R1.fastq.gz $SCR/raw/2505_2.1.fq.gz
cp $MEER_DIR/20250603_DNASeq_PE150/2505-Ll-RAD-P2_R2.fastq.gz $SCR/raw/2505_2.2.fq.gz

cp $MEER_DIR/20250603_DNASeq_PE150/2505-Ll-RAD-P3_R1.fastq.gz $SCR/raw/2505_3.1.fq.gz
cp $MEER_DIR/20250603_DNASeq_PE150/2505-Ll-RAD-P3_R2.fastq.gz $SCR/raw/2505_3.2.fq.gz

cp $MEER_DIR/20250603_DNASeq_PE150/2505-Ll-RAD-P4_R1.fastq.gz $SCR/raw/2505_4.1.fq.gz
cp $MEER_DIR/20250603_DNASeq_PE150/2505-Ll-RAD-P4_R2.fastq.gz $SCR/raw/2505_4.2.fq.gz

cp $MEER_DIR/20250603_DNASeq_PE150/2505-Ll-RAD-P5_R1.fastq.gz $SCR/raw/2505_5.1.fq.gz
cp $MEER_DIR/20250603_DNASeq_PE150/2505-Ll-RAD-P5_R2.fastq.gz $SCR/raw/2505_5.2.fq.gz
