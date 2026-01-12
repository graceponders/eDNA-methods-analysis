This is the visualisation component of an eDNA metabarcoding project looking at invertebrate community assembly in subtropical peatlands in southeast Queensland. 

Bioinformatics methodology; DNA was extracted and amplified from water samples using Zhan et al. (2014)'s v4 18S primers and sequenced using the Illumina MiSeq platform (Novogene).
Reads were processed in Galaxy Australia using with the LotuS2 pipeline to merge, demultiplex, and remove chimeras, generating operational taxonomic units (OTUs).
Resulting OTUs were queried against the SILVA SSU/LSU database as well as an inhouse custom 18S rRNA invertebrate database. 
In a separate R script, the two OTU classifications were joined by accession number and filtered by minimum percent identity as follows: ≥98% = genus-level (custom database), ≥95% = family-level (custom database), ≥90% = order-level (custom database), <90% = SILVA annotation.
OTUs with taxonomic assignment confidence <85% were discarded.

This script visualises community structure using Bray–Curtis nMDS ('vegan') and compares detections across methods using log-corrected abundance heat maps (eDNA vs. field surveys).
All figures generated using 'ggplot2' with the 'viridis' colour palate. 
For further annotation of figures, phylogenetic trees were generated using PhyloT and visualized in iTOL.
