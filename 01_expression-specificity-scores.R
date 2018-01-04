
# Setup 



"%&%" <- function(a,b) paste0(a,b)

library("tidyverse")
library("data.table")

serv.dir <- "/well/mccarthy/users/jason/"
work.dir <- serv.dir %&% "projects/t2d_classification/"
afile.dir <- work.dir %&% "analysis_files/"
gtex.dir <- serv.dir %&% "datasets/GTEx/v7/"

args <- commandArgs(trailingOnly=TRUE)

my.start <- args[1]
my.end <- args[2]

ensid.v2.df <- read.table(afile.dir%&%"gencode.v19.gene_lengths.txt",header=TRUE,
                          stringsAsFactors=FALSE)


islet.tpm.df <- fread(afile.dir%&%"islet_expression_tpm.txt")
mus.tpm.df <- fread(afile.dir%&%"muscle"%&%"_expression_tpm.txt")  
adi.tpm.df <- fread(afile.dir%&%"adipose"%&%"_expression_tpm.txt")  
liv.tpm.df <- fread(afile.dir%&%"liver"%&%"_expression_tpm.txt")  


build_ess_df <- function(my.start,my.end){
  ess.df <- c() 
  ens.vec <- islet.tpm.df$GeneID
  pb <- txtProgressBar(min=0,max=(my.end-my.start),style=3)
  for (i in my.start:my.end){
    setTxtProgressBar(pb,i)
    ens <- ens.vec[i]
    ensname <- filter(islet.tpm.df,GeneID==ens)$GeneName
    islet.vec <- filter(islet.tpm.df,GeneID==ens) %>% dplyr::select(.,-one_of("GeneID","GeneName")) %>% as.numeric(.)
    mus.vec <- filter(mus.tpm.df,GeneID==ens) %>% dplyr::select(.,-one_of("GeneID","GeneName")) %>% as.numeric(.)
    adi.vec <- filter(adi.tpm.df,GeneID==ens) %>% dplyr::select(.,-one_of("GeneID","GeneName")) %>% as.numeric(.)
    liv.vec <- filter(liv.tpm.df,GeneID==ens) %>% dplyr::select(.,-one_of("GeneID","GeneName")) %>% as.numeric(.)
    tot <- median(islet.vec) + median(mus.vec) + median(adi.vec) + median(liv.vec) 
    islet.score <- median(islet.vec) / tot 
    muscle.score <- median(mus.vec) / tot 
    adipose.score <- median(adi.vec) / tot 
    liver.score <- median(liv.vec) / tot 
    build.df <- data.frame(GeneID=ens,GeneName=ensname,islet.score,muscle.score,adipose.score,liver.score,stringsAsFactors=FALSE)
    ess.df <- rbind(ess.df,build.df)
  }
  return(ess.df)
}

ess.df <- build_ess_df(my.start,my.end)

write.table(ess.df,file=afile.dir %&%"expression_specificity_scores_indices_" %&% my.start %&% "-" %&% my.end %&% ".txt",
            sep="\t",quote=F,row.names=F)

