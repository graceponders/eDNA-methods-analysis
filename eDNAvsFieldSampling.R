#eDNA community analysis 
  #Grace S, May 2026 

#Making figures - heatmap comparing eDNA to field surveys and PCoA of communities at each site. 

#nMDS
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(tidyverse)
library(pairwiseAdonis)

# 1: heatmaps -------------------------------------------------------------

# Read data
dat <- read.csv("heatmap_detections.csv") #this is a spreadsheet with the relative taxon detection (0-1) at each site by method eDNA vs method field surveys

#ggplot like long data
dat_long <- dat %>%
  pivot_longer(
    cols = !Taxa,       # pivot all columns except Taxa
    names_to = "Site",
    values_to = "Abundance") %>%
  mutate(Taxa = factor(Taxa, levels = unique(Taxa)))  # keep original Class order

# plot heatmap
ggplot(dat_long, aes(x = Site, y = Taxa, fill = Abundance)) +
  geom_tile() +
  scale_fill_viridis_c(option = "D", direction = 1) +
  theme_classic() 

#done - add annotations externally. 


# 2: PCoA -----------------------------------------------------------------
dat <- read.csv("OTU.csv") #this is my LotuS2 OTU abundance matrix containing all ~2000 ASV detections across all 9 samples. 

head(dat)
#1910 obs of 9 variables 

#vegan likes samples as rows
tdat <- t(dat)

# Calculate Bray-Curtis dissimilarity
bray <- vegdist(tdat, method = "bray")

# Run PCoA
pcoa <- cmdscale(bray, k = 2, eig = TRUE)

# Extract coordinates
pcoa_scores <- as.data.frame(pcoa$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")

# Add site grouping
pcoa_scores$Site <- c("Duck Creek", "Duck Creek", "Duck Creek",
                      "Dilli Swamp", "Dilli Swamp", "Dilli Swamp",
                      "Red Lagoon", "Red Lagoon", "Red Lagoon")
pcoa_scores$Sample <- rownames(pcoa_scores)

# Calculate % variance explained by each axis
eig_percent <- round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)

# Plot
ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2, colour = Site, label = Sample)) +
  geom_point() +
  labs(x = paste0("PCoA1 (", eig_percent[1], "%)"),
       y = paste0("PCoA2 (", eig_percent[2], "%)")) +
  theme_classic() 

#Dilli Swamp is quite distinct from the others; some overlap between Duck Creek and Red Lagoon.
#Fairly high percent variance: PCoA1 49.8%; PCoA2 19.8%.

#Run permanova  

# Define site groupings
groups <- c("DC", "DC", "DC", "DS", "DS", "DS", "RL", "RL", "RL")

# PERMANOVA
permanova <- adonis2(bray ~ groups, permutations = 999)
print(permanova)

# PERMDISP
disp <- betadisper(bray, groups)
permutest(disp)

#Also run pairwise analysis 

pairwise <- pairwise.adonis(bray, groups, p.adjust.m = "bonferroni")
print(pairwise)

#Differences between sites are not statistically significant due to the sample size. 
