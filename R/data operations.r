# source: https://microbiome.github.io/tutorials/Preprocessing.html
# using practice data from http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html

library(phyloseq)
library(microbiome)
library(knitr)
data(atlas1006)   
# Rename the example data (which is a phyloseq object)
pseq <- atlas1006


##### SUMMARIZING THE DATA OF ATLAS1006 #####
# provides a summary of min, max, total, average, and median number of reads.
summarize_phyloseq(pseq)

##### RETRIEVE DATA ELEMENTS FROM THE DATA #####
# A phyloseq object contains OTU table (taxa abundances), sample metadata, taxonomy table (mapping between OTUs and higher-level taxonomic classifications), and phylogenetic tree (relations between the taxa).


# Pick metadata as data.frame:
meta <- meta(pseq)

# Taxonomy table:
taxonomy <- tax_table(pseq)

# Abundances for taxonomic groups (‘OTU table’) as a TaxaxSamples matrix:
# Absolute abundances
otu.absolute <- abundances(pseq)

# Relative abundances
otu.relative <- abundances(pseq, "compositional")

# Total read counts:
reads_sample <- readcount(pseq)

# check for first 5 samples
reads_sample[1:5]

# Add read per sample to phyloseq object metadata.
sample_data(pseq)$reads_sample <- reads_sample

# reads_sample is add to the last column in sample_data of pseq object.
head(meta(pseq)[,c("sample", "reads_sample")])

# Melting phyloseq data for easier plotting:
df <- psmelt(pseq)
kable(head(df))

##### SAMPLE OPERATIONS #####
# Sample names and variables
head(sample_names(pseq))

# Total OTU abundance in each sample
s <- sample_sums(pseq)

# Abundance of a given species in each sample
head(abundances(pseq)["Akkermansia",])

# Select a subset by metadata fields:
pseq.subset <- subset_samples(pseq, nationality == "AFR")

# Select a subset by providing sample names:
# Check sample names for African Females in this phyloseq object
s <- rownames(subset(meta(pseq), nationality == "AFR" & sex == "Female"))

# Pick the phyloseq subset with these sample names
pseq.subset2 <- prune_samples(s, pseq)

# Pick samples at the baseline time points only:
pseq0 <- baseline(pseq)

##### DATA TRANSFORMATIONS #####
# The microbiome package provides a wrapper for standard sample/OTU transforms.Relative abundances (note that the input data needs to be in absolute scale, not logarithmic!):
pseq.compositional <- microbiome::transform(pseq, "compositional")

# The CLR (“clr”) transformation is also available, and comes with a pseudocount to avoid zeroes. 
data(dietswap)
x <- dietswap

# Compositional data 
x2 <- microbiome::transform(x, "compositional")


##### VARIABLE OPERATIONS #####
# Sample variable names
sample_variables(pseq)

# Pick values for a given variable
head(get_variable(pseq, sample_variables(pseq)[1]))



# Assign new fields to metadata
# Calculate diversity for samples
div <- microbiome::alpha(pseq, index = "shannon")

# Assign the estimated diversity to sample metadata
sample_data(pseq)$diversity <- div

##### TAXA OPERATIONS #####
# Number of taxa
n <- ntaxa(pseq)

# Most abundant taxa
topx <- top_taxa(pseq, n = 10)

# Names
ranks <- rank_names(pseq)  # Taxonomic levels
taxa  <- taxa(pseq)        # Taxa names at the analysed level

# Subset taxa
pseq.bac <- subset_taxa(pseq, Phylum == "Bacteroidetes")

#Prune (select) taxa:
# List of Genera in the Bacteroideted Phylum
taxa <- map_levels(NULL, "Phylum", "Genus", pseq)$Bacteroidetes

# With given taxon names
ex2 <- prune_taxa(taxa, pseq)

# Taxa with positive sum across samples
ex3 <- prune_taxa(taxa_sums(pseq) > 0, pseq)

# Filter by user-specified function values (here variance):
f <- filter_taxa(pseq, function(x) var(x) > 1e-05, TRUE)

# List unique phylum-level groups:
head(get_taxa_unique(pseq, "Phylum"))


# Pick the taxa abundances for a given sample:
samplename <- sample_names(pseq)[[1]]

# Pick abundances for a particular taxon
tax.abundances <- abundances(pseq)[, samplename]

##### MERGING OPERATIONS #####
# Aggregate taxa to higher taxonomic levels. This is particularly useful if the phylogenetic tree is missing. When it is available, see merge_samples, merge_taxa and tax_glom).
pseq2 <- aggregate_taxa(pseq, "Phylum") 

# Merge the desired taxa into “Other” category. Here, we merge all Bacteroides groups into a single group named Bacteroides.
pseq2 <- merge_taxa2(pseq, pattern = "^Bacteroides", name = "Bacteroides") 

# Merging phyloseq objects with the phyloseq package:
merge_phyloseq(pseqA, pseqB)

##### Joining otu/asv table and taxonomy in one data frame #####
library(dplyr) 
library(microbiome)
data("atlas1006") # example data from microbiome pkg
x <-atlas1006
asv_tab <- as.data.frame(abundances(x)) # get asvs/otus
asv_tab$asv_id <- rownames(asv_tab) # add a new column for ids

#tax_tab <- as.data.frame(tax_table(x)) # get taxonomy note: can be slow
tax_tab <- as(x@tax_table,"matrix") # get taxonomy note as matrix
tax_tab <- as.data.frame(tax_tab) # convert to data frame
tax_tab$asv_id <- rownames(tax_tab) # add a new column for ids
asv_tax_tab <- tax_tab %>% 
  left_join(asv_tab, by="asv_id") # join to get taxonomy and asv table

head(asv_tax_tab)[,1:8]

##### RAREFACTION #####
pseq.rarified <- rarefy_even_depth(pseq)

##### TAXONOMY #####
# Convert between taxonomic levels (here from Genus (Akkermansia) to Phylum (Verrucomicrobia):
m <- map_levels("Akkermansia", "Genus", "Phylum", tax_table(pseq))
print(m)

##### METADATA #####
# Visualize frequencies of given factor (sex) levels within the indicated groups (group):
p <- plot_frequencies(sample_data(pseq), "bmi_group", "sex")
print(p)

# Custom functions are provided to cut age or BMI information into discrete classes.
group_bmi(c(22, 28, 31), "standard")
group_age(c(17, 41, 102), "decades")


# Add metadata to a phyloseq object. For reproducibility, we just use the existing metadata in this example, but this can be replaced by another data.frame (samples x fields).

# Example data
data(dietswap)
pseq <- dietswap

# Pick the existing metadata from a phyloseq object
# (or retrieve this from another source)
df <- meta(pseq)

# Merge the metadata back in the phyloseq object
pseq2 <- merge_phyloseq(pseq, sample_data(df))