cval <- 2.2
genomeSize <- cval * (0.978*10^9)  ## Genome size in bp (e.g., reference genome size)
cutSites <- genomeSize/4^8         ## Number of enzyme cut sites (If PstI use 4^6, if SbfI use 4^8)
radTags <- cutSites*2              ## Number of RAD tags created from enzyme cuts
readLength <- 150*2                ## Read length (e.g., 150 for 150bp). Include *2 if paired reads
basesCovered <- radTags*readLength ## Number of bases to sequence based on read length (if 2x150 use 300)
laneReads <- 1000000000            ## Reads per lane minus PhiX spike-in (i.e., 450000000 for HiSeq X, minus 50000000 for PhiX)
qualityFilter <- 0.6               ## Reads passing process_radtags
dupRemoval <- 0.75                 ## Proportion of reads retained after removing PCR duplicates
desiredDepth <- 40                 ## Desired sequencing depth
nSamples <- 446

### Calculate proportion of genome to be sequenced
basesCovered/genomeSize

### Number of raw reads to sequence one individual to desiredDepth at one locus
rawReadsNeeded.locus <- (desiredDepth / dupRemoval) / qualityFilter

### Number of raw reads to sequence per individual to achieve desiredDepth
rawReadsNeeded.ind <- rawReadsNeeded.locus * radTags

### Number of individuals that can fit on one lane to achieve desired depth at all loci
indPerLane <- laneReads / rawReadsNeeded.ind

### Lanes needed for samples
nSamples/indPerLane


rawReadsNeeded.ind * nSamples
## Reads needed = 2,000,000,000