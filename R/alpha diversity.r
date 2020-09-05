# source: https://microbiome.github.io/tutorials/Diversity.html
# load example data : global ecosystem state variables

library(microbiome)
library(knitr)
data(dietswap)
pseq <- dietswap

# global indicators
# A comprehensive list of global indicators of the ecosystem state can be obtained as follows. This includes various measures of richness, evenness, diversity, dominance, and rarity with default parameters. See the individual functions for more options regarding parameter tuning.
tab <-microbiome::alpha(pseq, index = "all")
kable(head(tab))


# alpha diversity 
# This returns a table with selected diversity indicators.
tab <-microbiome::alpha(pseq, index = "all")
kable(head(tab))

# richness
# This returns observed richness with given detection threshold(s).
tab <- richness(pseq)
kable(head(tab))

# dominance
# The dominance index refers to the abundance of the most abundant species. Various dominance indices are available (see the function help for a list of options).
# Absolute abundances for the single most abundant taxa in each sample
tab <- dominance(pseq, index = "all")
kable(head(tab))

# We also have a function to list the dominating (most abundant) taxa in each sample.
dominant(pseq)

# rarity and low abundance
# The rarity indices quantify the concentration of rare or low abundance taxa. Various rarity indices are available (see the function help for a list of options).
tab <- rarity(pseq, index = "all")
kable(head(tab))

# coverage
# The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).
tab <- coverage(pseq, threshold = 0.5)
kable(head(tab))

# core abundance
# The core_abundance function refers to the relative proportion of the core species. Non-core abundance provides the complement (1-x; see rare_abundance).
tab <- core_abundance(pseq, detection = .1/100, prevalence = 50/100)

# Gini index
# Gini index is a common measure for inequality in economical income. The inverse gini index (1/x) can also be used as a community diversity measure.
tab <- inequality(pseq)

# evenness
# Various evenness measures are also available.
tab <- evenness(pseq, "all")
kable(head(tab))




