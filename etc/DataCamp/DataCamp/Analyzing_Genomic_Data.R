# Load the BiocInstaller package
library(BiocInstaller)

# Check R version
version

# Explicit syntax to check the Bioconductor version
biocVersion()

# Load the BSgenome package
library(BSgenome)


# Check the version of the BSgenome package
packageVersion("BSgenome")

# Investigate the a_genome using show()
show(a_genome)

# Investigate some other accesors
organism(a_genome)
provider(a_genome)
seqinfo(a_genome)

# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Get the head of seqnames and tail of seqlengths for yeastGenome
head(seqnames(yeastGenome))
tail(seqlengths(yeastGenome))

# Print chromosome M, alias chrM
yeastGenome$chrM

# Count characters of the chrM sequence
nchar(yeastGenome$chrM)

# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Get the first 30 bases of each chromosome
getSeq(yeastGenome, end = 30)

# Load packages
library(Biostrings)

# Check the alphabet of the zikaVirus
alphabet(zikaVirus)

# Check the alphabetFrequency of the zikaVirus
alphabetFrequency(zikaVirus)

# Check alphabet of the zikaVirus using baseOnly = TRUE
alphabet(zikaVirus, baseOnly = TRUE)


# read in with
readDNAStringSet()


# Unlist the set, select the first 21 letters, and assign to dna_seq
dna_seq <- subseq(unlist(zikaVirus), end = 21)
dna_seq

# Transcribe dna_seq into an RNAString object and print it
rna_seq <- RNAString(dna_seq) 
rna_seq

# Translate rna_seq into an AAString object and print it
aa_seq <- translate(rna_seq)
aa_seq


# Unlist the set, select the first 21 letters, and assign to dna_seq
dna_seq <- subseq(unlist(zikaVirus), end = 21)
dna_seq

# Transcribe and translate dna_seq into an AAString object and print it
aa_seq <- translate(dna_seq)
aa_seq


#### SEQUENCE HANDLING ########################################################
# Create zikv with one collated sequence using zikaVirus
zikv <- unlist(zikaVirus)

# Check the length of zikaVirus and zikv
length(zikaVirus)
length(zikv)

# Check the width of zikaVirus
width(zikaVirus)

# Subset zikv to only the first 30 bases
subZikv <- subseq(zikv, end = 30)
subZikv


