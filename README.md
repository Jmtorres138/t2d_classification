# t2d_classification
**Goal:** Categorize T2D-associated loci according to tissue specificity

**Method A:** Islet, muscle, liver, adipose tissue chromatin states from Varshney _et_ al. 2016 are parititioned into tissue specific, islet shared, and insulin-responsive peripheral tissue shared annotations. fgwas was performed with CDS and distance to TSS (0 5000 dist) information to obtain fgwas credible sets. Scores for tissue classification based on PPA proportions attributable to tissues. Additional steps required to divvy up PPA for loci largely determined by CDS, distance, and shared annotations. 

**Method B:** A full set of chromatin states for each of the four key tissue of interests are used in individual fgwas runs to obtain to best joint model for each tissue. Note: this set includes genomic annotations (some of which redundant) and distance. Scores for tissue classification based on changes in credible set features (i.e. size, max SNP PPA) compared to null fgwas credible sets. 

**Method C:** GenoSkyline approach 
