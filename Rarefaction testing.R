#Rarefaction analysis on OTU dataset

  #Grace S, May 2026


  #In this page, I am testing rarefaction curves for the OTU richness vs Read Depth and sampling effort for each sample. 

library(vegan)
library(iNEXT)

# Bring in OTU data
OTU <- read.table("OTU.txt", sep = "\t", header = TRUE, row.names = 1) #this is the LotuS2 OTU output
head(OTU)

# Vegan likes samples as rows
tOTU <- t(OTU)

# Check read depth per sample
rowSums(tOTU)

#Total reads is very consistent across samples:
#DC1     DC2     DC3     DS1     DS2     DS3     RL1     RL2     RL3 
#904323  774877  918396  923767  904496  942479  802401  814962 1127996 

# Plot rarefaction curves
rarecurve(tOTU, 
          step = 50,        # low(ish) granularity (should be fairly smooth but will take a while to run)
          label = TRUE,      # labels each curve with sample name
          xlab = "Read Depth",
          ylab = "OTU Richness")

#To me, these curves look very normal. They plateau early and are relatively consistent. 

#I'm going to rarefy to the smallest sample, which is DC2 at 774877. 

rOTU <- rrarefy(tOTU, sample = 774877)
rowSums(rOTU)

#To confirm this worked: 
#DC1    DC2    DC3    DS1    DS2    DS3    RL1    RL2    RL3 
#774877 774877 774877 774877 774877 774877 774877 774877 774877 

# Compare OTU richness before and after rarefaction
comparison <- data.frame(
  ReadDepth_before = rowSums(OTU_t),
  Richness_before = specnumber(tOTU),
  ReadDepth_after = rowSums(rOTU),
  Richness_after = specnumber(rOTU),
  OTUs_lost = specnumber(tOTU) - specnumber(rOTU)
)

print(comparison)

#Sample ReadDepth_before Richness_before ReadDepth_after Richness_after OTUs_lost
#DC1    DC1           904323            1641          774877           1625        16
#DC2    DC2           774877            1457          774877           1457         0
#DC3    DC3           918396            1618          774877           1611         7
#DS1    DS1           923767            1505          774877           1499         6
#DS2    DS2           904496            1465          774877           1457         8
#DS3    DS3           942479            1147          774877           1136        11
#RL1    RL1           802401            1331          774877           1328         3
#RL2    RL2           814962            1398          774877           1392         6
#RL3    RL3          1127996            1443          774877           1398        45

#To me, this indicates that diversity estimates are not an artifact of sequencing depth. 

#Additionally, we will test whether additional samples would have been likely to detect a greater number of species/OTUs. 
#To do this, I will use the iNext package. 

#This will use the untransposed OTU table from earlier. 

#First: would additional replicates per site increase detection?  

#iNEXT inputs a table of OTU frequencies:
freq_matrix <- rowSums(OTU > 0)
input <- c(ncol(OTU), freq_matrix)

head(input)

#This can then be graphed using the iNEXT() function: 

out <- iNEXT(x = input, datatype = "incidence_freq")

ggiNEXT(out, type = 1, se = TRUE)

out$DataInfo

#Based on these outputs, we can see the curve plateus, giving a very high SC value (0.994)
  #((this means that an overall, an additional OTU would be 99.4% one that we have already detected))

#It is more fair to assess this by site, to determine whether 3 biological replicates is adequate. 

#Not sure what the most efficient way to do this would be; I'm making a separate input for each site and will pool at the end. 

OTU <- OTU > 0
  #DUCK CREEK
dc <- OTU[, c("DC1", "DC2", "DC3")]
dc_freq <- rowSums(dc)
dc_freq <- dc_freq[dc_freq > 0] # remove OTUs never detected
dc_iNEXT <- c(3, dc_freq)

  #DILLI SWAMP
ds <- OTU[, c("DS1", "DS2", "DS3")]
ds_freq <- rowSums(ds)
ds_freq <- ds_freq[ds_freq > 0]
ds_iNEXT <- c(3, ds_freq)

  #RED LAGOON
rl <- OTU[, c("RL1", "RL2", "RL3")]
rl_freq <- rowSums(rl)
rl_freq <- rl_freq[rl_freq > 0]
rl_iNEXT <- c(3, rl_freq)


site_list <- list(DC = dc_iNEXT, DS = ds_iNEXT, RL = rl_iNEXT)


out <- iNEXT(site_list, datatype = "incidence_freq", endpoint = 10) #I've chosen 10 to give the extrapolation more room. 

ggiNEXT(out, type = 1, se = TRUE)
out$DataInfo

#Curves do not reach a plateu until ~5 biological replicates; this is expected based on literature. SC values are still high. 
#This would indicate that our sample size was adequeate for detecting overall community patterns, but we likely missed rare species. 