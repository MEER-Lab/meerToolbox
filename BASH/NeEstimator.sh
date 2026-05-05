
#Define alias for project root directory (e.g., 2505_mnBurbotStructure)
PROJ_ID=2514_lakeSuperiorBnt
YEAR=2025

MEER_DIR=/mnt/research/meerLab/$YEAR/$PROJ_ID 
SCR=/mnt/gs21/scratch/homolaj1/$PROJ_ID 

cd $MEER_DIR/out

dos2unix clusters.2514.gen
sed -i 's/\t/ /g' clusters.2514.gen
LC_ALL=C tr -cd '\11\12\15\40-\176' < clusters.2514.gen > clusters.2514.clean.gen
mv clusters.2514.clean.gen clusters.2514.gen
/mnt/research/meerLab/software/NeEstimator2.X/Ne2x ./clusters.2514.gen