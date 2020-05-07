

"%&%" <- function(a,b) paste0(a,b)
library("tidyverse")
library("data.table")
library("GenomicRanges")

fuse.path <- "/home/jason/science/servers/"
serv.dir <- fuse.path %&% "FUSE5/"
serv.dir0 <- fuse.path %&% "FUSE/"

proj.dir <- serv.dir %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "tactical_analyses/"
analysis.dir <- work.dir %&% "analysis_files/"

cred.df <- fread(work.dir %&% "genetic_credible_sets/gencred.txt")
class.df <- fread(analysis.dir %&% "classified-loci_weighted_with-shared.txt")



varsh.df <- serv.dir0 %&% "reference/chromatin_segmentation/varshney_2016/bed_files/Islets.chromatinStates.bed" %>% 
  fread(.)
thurn.df <- serv.dir0 %&% "reference/islet/Pancreat_islet_15_dense.reformatted_colours.bed" %>% 
  fread(.)


#Get proportions



(varsh.df$V4 %>% table(.)) / dim(varsh.df)[1] 
(thurn.df$V4 %>% table(.)) / dim(varsh.df)[1] 

