This is the visualisation component of an eDNA metabarcoding project looking at invertebrate community assembly in subtropical peatlands in southeast Queensland. 

Bioinformatics methodology; DNA was extracted and amplified from water samples using Zhan et al. (2014)'s v4 18S primers and sequenced using the Illumina MiSeq platform (Novogene).
Reads were processed in Galaxy Australia using with the LotuS2 pipeline to merge, demultiplex, and remove chimeras, generating operational taxonomic units (OTUs).
Resulting OTUs were queried against the SILVA SSU/LSU database as well as an inhouse custom 18S rRNA invertebrate database. 
The two OTU classifications were joined by accession number and filtered by minimum percent identity as follows: ≥98% = genus-level (custom database), ≥95% = family-level (custom database), ≥90% = order-level (custom database), <90% = SILVA annotation.
OTUs with taxonomic assignment confidence <85% were discarded.

There are three sections to this page: Firstly, plotting log-corrected relative detection of taxa detected with metabarcoding surveys vs field surveys. Secondly, plotting community differences by site using PCoA based on Bray-Curtis nMDS ('vegan').
Finally, the Rarefaction Analyses sheet tests the effect of sequencing depth and sampling effort on OTU detection ('vegan' and 'iNEXT').
All figures generated using 'ggplot2' with the 'viridis' colour palate. 
For further annotation of figures, phylogenetic trees were generated using PhyloT and visualized in iTOL.
