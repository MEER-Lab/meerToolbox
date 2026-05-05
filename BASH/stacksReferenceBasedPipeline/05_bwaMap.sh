#==================================================================================================
#   File: 05_bwaMap.sh
#   Date: December 15, 2025
#   Description: Map reads to reference genome using bwa-mem
#--------------------------------------------------------------------------------------------------
#	Authors: Jared Homola
#==================================================================================================

#Define alias for project root directory
PROJ_ID=2505_mnBurbotStructure
YEAR=2025

MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID 
SCR=/mnt/gs21/scratch/homolaj1/$PROJ_ID 

mkdir $SCR/mapped
mkdir $SCR/mapped/filtered
mkdir $SCR/mapped/filtered/sorted
mkdir $MEER_DIR/stdout/bwaMemMap
mkdir $MEER_DIR/shell/bwaMemMap

module load BWA/0.7.17-GCCcore-12.3.0 SAMtools/1.18-GCC-12.3.0

cd $MEER_DIR/shell/bwaMemMap
ls $SCR/cloneFilter | grep ".1.fq.gz" | sed 's/\.1\.1\.fq\.gz//g' | sort | uniq | while read -r LINE
do
echo '#!/bin/sh 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 2:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH -J '$LINE'.bwaMemMap
#SBATCH -o '$MEER_DIR'/stdout/bwaMemMap/bwaMemMap.'$LINE'.o

cd '$SCR'/cloneFilter

bwa mem -R "@RG\tID:'$LINE'\tSM:'$LINE'\tPL:ILLUMINA\tLB:LB1" '$MEER_DIR'/referenceGenome/GCA_026745465.1_ASM2674546v1_genomic.fna ./'$LINE'.1.1.fq.gz ./'$LINE'.2.2.fq.gz | samtools view -Sb - > ../mapped/'$LINE'.pe.sam 

samtools view -q 20 -h -f2 -F2308 ../mapped/'$LINE'.pe.sam  | samtools view -Sb > ../mapped/filtered/'$LINE'.pe.filtered.bam
samtools sort --threads 32 ../mapped/filtered/'$LINE'.pe.filtered.bam -o ../mapped/filtered/sorted/'$LINE'.pe.filtered.sorted.bam

scontrol show job ${SLURM_JOB_ID}' > ./bwaMemMap."$LINE".sh

sbatch ./bwaMemMap."$LINE".sh

done


### Log results of mapping

#Define alias for project root directory
RUN_PATH=/mnt/home/homolaj1/gobies/lib18

cd $RUN_PATH/SHELL

echo '#!/bin/sh 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G 
#SBATCH -J bwaFlagstats
#SBATCH -o '$RUN_PATH'/bwaFlagstats.o

module purge
module load GCC/8.3.0 SAMtools/1.10

cd '$SCR'/gobies/lib18/mapped

ls | grep ".pe.bam" | sort | uniq | while read -r LINE
do
   echo "$LINE" >> '$RUN_PATH'/OUT/mappingResults.txt
   samtools flagstat "$LINE" >> '$RUN_PATH'/OUT/mappingResults.txt
done' > bwaFlagstats.sh

sbatch bwaFlagstats.sh
