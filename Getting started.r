# installing and get packages from Bioconductor (open source software for bioinformatics)
# source: https://bioconductor.org/packages/devel/bioc/html/microbiome.html

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version='devel')

# make sure this package installs through the bioconductor package
BiocManager::install('microbiome')

# Loading the package in R (after installation from Bioconductor)
library(microbiome)

# This microbiome package relies on external phyloseq datasets, which should already have been installed through BioManager::install('microbiome'), which you can use as tutorials (optional)

# data from http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html
data(atlas1006) 


# or you can use http://joey711.github.io/phyloseq/import-data
data(GlobalPatterns)


