# source: https://microbiome.github.io/tutorials/Composition.html

# Read example data from a diet swap study:
library(microbiome)
library(dplyr)
data(dietswap)

# Make sure we use functions from correct package
transform <- microbiome::transform

# Merge rare taxa to speed up examples
pseq <- transform(dietswap, "compositional")
pseq <- aggregate_rare(pseq, level = "Genus", detection = 1/100, prevalence = 50/100)

# Pick sample subset
library(phyloseq)
pseq2 <- subset_samples(pseq, group == "DI" & nationality == "AFR" & timepoint.within.group == 1)

# Normal western adults
data(atlas1006)
pseq3 <- atlas1006 %>%
          subset_samples(DNA_extraction_method == "r") %>%
          aggregate_taxa(level = "Phylum") %>%  
          microbiome::transform(transform = "compositional")

# Composition barplots
# Same with compositional (relative) abundances; for each sample (left), or averafged by group (right).

# Try another theme
# from https://github.com/hrbrmstr/hrbrthemes
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
theme_set(theme_bw(21))
p <- pseq3 %>%
    plot_composition(sample.sort = "Firmicutes", otu.sort = "abundance") +
         # Set custom colors
          scale_fill_manual(values = default_colors("Phylum")[taxa(pseq3)]) +
      scale_y_continuous(label = scales::percent)

print(p)


# Limit the analysis on core taxa and specific sample group
p <- plot_composition(pseq2,
              taxonomic.level = "Genus",
                      sample.sort = "nationality",
                      x.label = "nationality") +
     guides(fill = guide_legend(ncol = 1)) +
     scale_y_percent() +
     labs(x = "Samples", y = "Relative abundance (%)",
                                   title = "Relative abundance data",
                                   subtitle = "Subtitle",
                                   caption = "Caption text.") + 
     theme_ipsum(grid="Y")
print(p)  


# Averaged by group
p <- plot_composition(pseq2,
                      average_by = "bmi_group", transform = "compositional")
print(p)


p <- NULL
Composition heatmaps
Heatmap for CLR-transformed abundances, with samples and OTUs sorted with the neatmap method:

p <- plot_composition(microbiome::transform(pseq, "compositional"),
                    plot.type = "heatmap",
                        sample.sort = "neatmap", otu.sort = "neatmap")
print(p)


# Plot taxa prevalence
# This function allows you to have an overview of OTU prevalences alongwith their taxonomic affiliations. This will aid in checking if you filter OTUs based on prevalence, then what taxonomic affliations will be lost.

data(atlas1006)

# Use sample and taxa subset to speed up example
p0 <- subset_samples(atlas1006, DNA_extraction_method == "r")

# Define detection and prevalence thresholds to filter out rare taxa
p0 <- core(p0, detection = 0.1/100, prevalence = 1/100)

# For the available taxonomic levels
plot_taxa_prevalence(p0, "Phylum", detection = 0.1/100)


# Amplicon data
# Also see phyloseq barplot examples.

# Here the data from Thompson, Luke R., et al. “A communal catalogue reveals Earth’s multiscale microbial diversity.” Nature 551.7681 (2017): 457-463. will be used which is stored as example in jeevanuDB

# Check the core microbiome page which shows how to read the your files into R and make a phyloseq object.

# Example data
library(microbiome)

# Try another theme
# from https://github.com/hrbrmstr/hrbrthemes
# you can install these if you don't have it already.
# devtools::install_github("hrbrmstr/hrbrthemes")
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(dplyr)
library(jeevanuDB)
ps1 <- emp_human
colnames(tax_table(ps1))


# As you can see the taxonomic classification is just lablled as “Rank1” … “Rank7”. We need to change this to proper designation and also do some formatting of the data. This can be a useful example for understanding simple file processing in R.

# In case you see the taxonomic classification is just lablled as “Rank1” … “Rank7” we can change it as follows

# First change the column names of the taxonomy table in phyloseq to following:
colnames(tax_table(ps1)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species" )
tax_table(ps1)[tax_table(ps1)[,"Domain"]== "NA", "Domain" ] <- "Unidentified_Domain"
tax_table(ps1)[tax_table(ps1)[,"Phylum"]== "p__", "Phylum" ] <- "p__Unidentified_Phylum"



# Composition barplots
# The compositon plots can be shown either as barplots or heatmaps. Both examples are show below.

# Plot counts abundance
# Now we can improve the plot further.
# Let’s try at Family level.
library(phyloseq)

# merge at family level.
# check how many samples are there
# Use only saliva samples 
ps1.saliva <- subset_samples(ps1, env_material == "saliva")
total_samples <- nsamples(ps1.saliva)
ps1.saliva.pruned <- prune_taxa(taxa_sums(ps1.saliva) >0, ps1.saliva)

# merge all taxa that are detected rare
pseq.fam <- aggregate_rare(ps1.saliva.pruned, level="Family", detection = 50, prevalence = 25/total_samples)

p.fam <- plot_composition(pseq.fam, sample.sort = NULL, 
                          otu.sort = NULL,
                          x.label = "empo_3", # sample type
                          plot.type = "barplot", 
                          verbose = FALSE) + 
  theme_bw() + scale_fill_brewer("Family", palette = "Paired")
  
# we can rotate x axis labels 
print(p.fam + theme(axis.text.x = element_text(angle = 90)))


# Plot relative abundance
pseq.famrel <- microbiome::transform(pseq.fam, "compositional")

p.famrel <- plot_composition(pseq.famrel, sample.sort = NULL, otu.sort = NULL,
  x.label = "empo_3", plot.type = "barplot", verbose = FALSE)

print(p.famrel)


# further improvements can be done as follows  

p.famrel <- plot_composition(pseq.famrel, 
                             sample.sort = NULL, 
                             otu.sort = NULL, 
                             x.label = "empo_3", 
                             plot.type = "barplot", 
                             verbose = FALSE) + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") + 
     scale_fill_brewer("Family", palette = "Paired")

print(p.famrel)




# Averaged by group

# Use all samples 
ps1 <- emp_human

# get relative abudance
ps1.rel <- microbiome::transform(ps1, "compositional")
ps1.fam.rel <-aggregate_rare(ps1.rel, level = "Family", detection = 0.005, prevalence = 0.5)

p <- plot_composition(ps1.fam.rel,
                      average_by = "empo_3") + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") 
print(p + scale_fill_brewer("Family", palette = "Paired") + theme_bw())


# Heatmap composition
# Use all samples 
ps1 <- emp_human
ps1.rel <-aggregate_rare(ps1, level = "Family", detection = 10, prevalence = 0.5)

pseq.famlog <- microbiome::transform(ps1.rel, "log10")

p.famrel.heatmap <- plot_composition(pseq.famlog, 
                             sample.sort = NULL, 
                             otu.sort = NULL, 
                             x.label = "empo_3", 
                             plot.type = "heatmap", 
                             verbose = FALSE)

print(p.famrel.heatmap)


Plot core taxa time trajectory
library(dplyr)

# select core
ps <- moving_pictures
table(meta(ps)$sample_type, meta(ps)$host_subject_id)


taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
# Filter the data to include only gut samples from M3 subject
ps.m3 <- subset_samples(ps, sample_type == "stool" & host_subject_id == "M3") 
#print(ps.m3)

ps.m3.rel <- microbiome::transform(ps.m3, "compositional")
pseq.core <- core(ps.m3.rel, detection = 0.001, prevalence = .95)

ps.stool.df <- psmelt(pseq.core)
#head(ps.stool.df)

# add genus name to ASVid
ps.stool.df <- ps.stool.df %>% 
  mutate(asv_gen= paste0(OTU, "-",Genus))

ps.stool.rel.plot <- ggplot(ps.stool.df) + 
  geom_line(aes(days_since_experiment_start, 
                Abundance, color = asv_gen)) +
  theme_bw() + 
  theme(legend.position="top") + 
  xlab("Days since experiment start") + 
  ylab("Relative abundance") + 
  scale_color_brewer("Core ASVs",palette = "Paired") +
  guides(col = guide_legend(ncol = 3, nrow = 3))

ps.stool.rel.plot 


# Highlight only one ASVs of interest.

ps.highlight.plot <- ggplot(ps.stool.df) + 
  geom_line(aes(days_since_experiment_start, 
                Abundance), color="grey80") 
				
# pick only data for ASV996-g__Faecalibacterium

asv996 <- subset(ps.stool.df, asv_gen =="ASV996-g__Faecalibacterium")

ps.highlight.plot <- ps.highlight.plot + 
  geom_line(data= asv996,aes(x=days_since_experiment_start, 
                y=Abundance, color=asv_gen)) +
  theme_bw() + 
  theme(legend.position="top") + 
  xlab("Days since experiment start") + 
  ylab("Relative abundance") + 
  scale_color_manual("Core ASVs",values="brown3") +
  guides(col = guide_legend(ncol = 3, nrow = 3))

ps.highlight.plot