#eDNA community analysis 
  #Grace 2025 

#eDNA community data is clean and ready to plot. Making figures - heatmap comparing eDNA to field surveys and PCA of communities at each site by method. 

#nMDS
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(tidyverse)

# Read data
dat <- read.csv("methods_comparison1.txt") #this is a spreadsheet with the species detection at each site by method eDNA vs method field surveys


# 1: heatmaps -------------------------------------------------------------

#ggplot like long data
dat_long <- dat %>%
  pivot_longer(
    cols = !Taxa,       # pivot all columns except Class
    names_to = "Site",
    values_to = "Value") %>%
  mutate(Taxa = factor(Taxa, levels = unique(Taxa)),  # keep original Class order
    Site  = factor(Site,  levels = colnames(dat)[-ncol(dat)]))  # keep original site order
  
# plot heatmap
ggplot(dat_long, aes(x = Site, y = Taxa, fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c(option = "D", direction = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#done - add annotations externally. 


# 2: nMDS -----------------------------------------------------------------

# Separate numeric data and Class
otu_matrix <- dat %>% 
  select(-Class) %>% 
  as.matrix()

rownames(otu_matrix) <- dat$Class

# Aggregate abundance by Class across replicates
otu_agg <- dat %>%
  group_by(Class) %>%
  summarise(across(X1DC:X2RL, sum)) %>%
  column_to_rownames("Class")

# Transpose: rows = sites, columns = classes
otu_t <- t(otu_agg)

# NMDS using Bray-Curtis dissimilarity
nmds <- metaMDS(otu_t, distance = "bray", k = 2, trymax = 100)
nmds

#extract scores
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$Site <- rownames(site_scores)

site_scores <- site_scores %>%
  mutate(
    Peatland = case_when(
      grepl("DC", Site) ~ "DC",
      grepl("DS", Site) ~ "DS",
      grepl("RL", Site) ~ "RL"),
    Method = case_when(
      grepl("X1", Site) ~ "X1",
      grepl("X2", Site) ~ "X2"))
    

# Plot
ggplot(site_scores,
       aes(x = NMDS1, y = NMDS2,
           shape = Method,
           color = Peatland)) +
  scale_shape_manual(values = c(
      X1 = 16,  # circle
      X2 = 17))+   # triangle
  scale_color_manual(
    values = c(
      RL = "#FDE726",  
      DC = "#440154",   
      DS = "#21908C"))+
  coord_equal() +
  theme_classic()

##PERMANOVA
Method <- factor(c("eDNA", "eDNA", "eDNA", "Field", "Field", "Field"))
Site <- factor(c("Duck Creek", "Dilli Swamp", "Red Lagoon",
                 "Duck Creek", "Dilli Swamp", "Red Lagoon"))

# Check
cbind(Method, Site)
permanova <- adonis2(otu_t ~ Method + Site, method = "bray", permutations = 999, by = "terms")
print(permanova)

#done!