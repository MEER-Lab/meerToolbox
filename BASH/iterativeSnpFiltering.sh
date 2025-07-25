
#Define alias for project root directory (e.g., 2505_mnBurbotStructure)
PROJ_ID=

MEER_DIR=/mnt/research/meerLab/2025/$PROJ_ID 
SCR=/mnt/gs21/scratch/homolaj1/$PROJ_ID 

cd $MEER_DIR/out

module load VCFtools/0.1.16-GCC-12.3.0 BCFtools/1.18-GCC-12.3.0

# === User-defined thresholds ===
STEPS=10
FINAL_LOCUS_MISSING=0.1     # Allow up to X% missing per locus
FINAL_INDV_MISSING=0.3      # Allow up to X% missing per individual
STATS_FILE="missingness_summary.txt"

# Convert to data present thresholds for vcftools
FINAL_LOCUS_MIN_DATA=$(echo "scale=4; 1 - $FINAL_LOCUS_MISSING" | bc)
FINAL_INDV_MIN_DATA=$(echo "scale=4; 1 - $FINAL_INDV_MISSING" | bc)

# Create working copy
WORKING_VCF="working_temp.vcf"
cp "${INPUT_VCF}" "${WORKING_VCF}"

# Step size for each threshold
LOCUS_STEP_SIZE=$(echo "scale=4; $FINAL_LOCUS_MIN_DATA / ($STEPS - 1)" | bc)
INDV_STEP_SIZE=$(echo "scale=4; $FINAL_INDV_MIN_DATA / ($STEPS - 1)" | bc)

# Initialize output
echo -e "Step\tLocusThreshold\tIndvThreshold\tNumIndividuals\tNumLoci" > "${STATS_FILE}"

# Function to log VCF stats
count_vcf_stats() {
    local vcf_file="$1"
    local step="$2"
    local locus_data="$3"
    local indv_data="$4"

    local num_indv
    num_indv=$(bcftools query -l "$vcf_file" | wc -l)

    local num_loci
    num_loci=$(grep -vc "^#" "$vcf_file")

    echo -e "${step}\t${locus_data}\t${indv_data}\t${num_indv}\t${num_loci}" >> "${STATS_FILE}"
}

# Step 0: stats on unfiltered VCF
count_vcf_stats "${WORKING_VCF}" "0" "0.0" "0.0"

# Iterative filtering loop
for ((i=1; i<=STEPS; i++)); do
    LOCUS_MIN_DATA=$(echo "scale=4; $LOCUS_STEP_SIZE * ($i - 1)" | bc)
    INDV_MIN_DATA=$(echo "scale=4; $INDV_STEP_SIZE * ($i - 1)" | bc)
    echo "Step $i: Locus threshold = ${LOCUS_MIN_DATA}, Individual threshold = ${INDV_MIN_DATA}"

    # 1. Filter loci first
    vcftools --vcf "${WORKING_VCF}" \
             --max-missing "$LOCUS_MIN_DATA" \
             --recode --recode-INFO-all \
             --out temp_locus_filtered > /dev/null 2>&1

    if [[ ! -f temp_locus_filtered.recode.vcf ]]; then
        echo "Step $i: Locus filtering failed. Exiting."
        exit 1
    fi

    # 2. Compute individual missingness based on filtered loci
    vcftools --vcf temp_locus_filtered.recode.vcf --missing-indv --out temp_missing > /dev/null 2>&1
    awk -v min_data="$INDV_MIN_DATA" 'NR>1 && (1 - $5) < min_data {print $1}' temp_missing.imiss > temp_remove.txt

    # 3. Remove individuals
    vcftools --vcf temp_locus_filtered.recode.vcf \
             --remove temp_remove.txt \
             --recode --recode-INFO-all \
             --out temp_step_filtered > /dev/null 2>&1

    if [[ ! -f temp_step_filtered.recode.vcf ]]; then
        echo "Step $i: Individual filtering failed. Exiting."
        exit 1
    fi

    # 4. Replace working file and save stats
    mv temp_step_filtered.recode.vcf "${WORKING_VCF}"
    count_vcf_stats "${WORKING_VCF}" "$i" "$LOCUS_MIN_DATA" "$INDV_MIN_DATA"
done

# Final cleanup and rename
mv "${WORKING_VCF}" "missingnessFiltered.vcf"
rm -f temp_missing* temp_step_* temp_locus_filtered.* temp_remove.txt

echo "Filtering complete. Final VCF: missingnessFiltered.vcf"
echo "Summary saved to: ${STATS_FILE}"

vcftools --vcf missingnessFiltered.vcf --min-meanDP 10 --recode --recode-INFO-all --out missingness_depthFiltered
grep -v "^#" missingness_depthFiltered.vcf | wc -l