#==================================================================================================
#   File: GTscore-setup.sh
#   Date: July 7, 2025
#   Description: Set everything up to genotype via GT-score
#--------------------------------------------------------------------------------------------------
#	Authors: Jared Homola
#==================================================================================================

#Define alias for project root directory
PROJ_ID=2523_waeBroodstock20242025
MEER_DIR=/mnt/research/meerLab/2025/$PROJ_ID 

# Create necessary directories
mkdir $MEER_DIR/genotyping

# Remove unnecessary files from demux directory
cd $MEER_DIR/demux
rm *_I* *_R2.fastq.gz Undetermined* unmatched*

# Delete empty files if needed
find ./ -type f -size -100c -delete


# Drop the _R1
for file in *_R1*; do mv -- "$file" "${file//_R1/}"; done

# Unzip everything
for f in *.gz; do
  gunzip "$f"
done

# Set up sampleFiles.txt for GTscore
ls $MEER_DIR/demux >> $MEER_DIR/demux/sampleFiles.txt



#### To run this perl script for the first time, you must...
#cpanm --local-lib=~/perl5 --notest Algorithm::Combinatorics
#export PERL5LIB=~/perl5/lib/perl5:$PERL5LIB
#export PATH=~/perl5/bin:$PATH
#source ~/.bashrc

module purge
module load Perl-bundle-CPAN/5.36.1-GCCcore-12.3.0

cd $MEER_DIR/demux

# Create primer and probe files for GTscore
perl /mnt/research/meerLab/software/GTscore/AmpliconReadCounter.pl \
    -p /mnt/research/meerLab/software/GTscore/GTscore_primer_probe_files/Svit_n400_primerProbe_file.txt \
    --files $MEER_DIR/demux/sampleFiles.txt 

# Move files to genotyping directory
cp $MEER_DIR/demux/* $MEER_DIR/genotyping/






