# Setup



"%&%" <- function(a,b) paste0(a,b)

library("tidyverse")
library("data.table")
library("Homo.sapiens")

serv.dir <- "/well/mccarthy/users/jason/"
proj.dir <- serv.dir %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "method_C/"
out.dir <- proj.dir %&% "method_C/analysis_files/"
res.dir <- work.dir %&% "genetic_credible_sets/" #"null_results/"
input.dir <- proj.dir %&% "analysis_files/"

sym.ids <- unique(keys(Homo.sapiens, keytype = "SYMBOL"))
sym.df <- select(Homo.sapiens,key=sym.ids,keytype="SYMBOL",
                 columns=c("ENSEMBL"))

ess.dir <- proj.dir %&% "analysis_files/"
ess.df <- fread(ess.dir %&% "expression_specificity_scores.txt")
#ess.df <- fread(ess.dir %&% "expression_specificity_scores-rntransform.txt")

weight.all.df <- fread(proj.dir %&%  "method_A/analysis_files/weight-enrich-all.txt")





#blk.df <- fread(res.dir %&% "null_results_blocks.txt")
cred.df <- fread(res.dir %&% "gencred.txt")
bed1.df <- fread(input.dir %&% "all_shared.bed")
bed2.df <- fread(input.dir %&% "shared_processed.txt")
shared.df <- rbind(bed1.df,bed2.df); rm(bed1.df); rm(bed2.df)
specific.df <- fread(input.dir %&% "specific_processed.txt")
states.df <- rbind(shared.df,specific.df)


# Append necessary information 

#mult.df <- fread(proj.dir %&%  "method_A/multi_results/results_func-cred-sets.txt")

#fgwas.path <- "/well/got2d/jason/projects/t2d-integration/fgwas/" %&% 
#  "diagram_hrc/cross_tissue/multi_tissue_joint_analysis/fgwas_input/" %&% 
#  "ukbb_diamante-euro.fgwas.gz"

#fgwas.df <- fread("cat " %&% fgwas.path %&% " | zmore")
#names(fgwas.df)[names(fgwas.df)=="distance_tss"] <- "distance_tss_0_5000"
#fgwas.df <- dplyr::select(fgwas.df,one_of(names(mult.df))) %>% 
#  dplyr::select(.,-one_of("CHR","POS","Z"))

#cred.df <- inner_join(cred.df,fgwas.df,by="SNPID")

# Divvy function

locid_to_symbol <- function(loc.id){
  return(filter(cred.df,Locus.ID==loc.id)$symbol %>% as.character(.) %>% unique(.))
}


map_tissue <- function(tiss){
  num <- ifelse(tiss=="islet",1,
                ifelse(tiss=="muscle",2,
                       ifelse(tiss=="adipose",3,4)))
  return(num)
}

handle_coding <- function(loc.id,ppa){
  sym <- locid_to_symbol(loc.id)
  sub <- filter(ess.df,GeneName==sym)
  score.vec <- c(0.25,0.25,0.25,0.25)
  divvy.vec <- ppa * score.vec
  return(divvy.vec)
}

handle_annotation <- function(annot.df, annot.name, ppa){
  sub <- filter(annot.df,grepl(pattern=annot.name,V4))
  if (dim(sub)[1]>1){
    print(sub)
  }
  if (dim(sub)[1]==0){
    divvy.vec <- c(0,0,0,0)
  } else if(grepl(pattern="_specific_",x=sub$V4)){
    tiss <- strsplit(x=sub$V4,split="_specific_")[[1]][1]
    num <- map_tissue(tiss)
    divvy.vec <- c(0,0,0,0)
    divvy.vec[num] <- ppa
  } else if(grepl(pattern="all_shared",x=sub$V4)){
    divvy.vec <- rep(ppa/4,4)
  } else if(grepl(pattern="shared.",x=sub$V4,fixed=TRUE)){
    vec <- strsplit(x=sub$V4,split=".",fixed=TRUE)[[1]]
    tiss.vec <- vec[2:(length(vec)-1)]
    ppa.val <- ppa/length(tiss.vec)
    divvy.vec <- c(0,0,0,0)
    index.vec <- map(tiss.vec,function(t){map_tissue(t)}) %>% as.integer(.)
    divvy.vec[index.vec] <- ppa.val
  } else{
    divvy.vec <- NA
  }
  return(divvy.vec)
}



divvy_ppa_snp <- function(loc.id,snp.id,mode="full",weights=FALSE){
  # set weights=NULL if unweighted
  sub.df <- filter(cred.df,Locus.ID==loc.id,SNPID==snp.id)
  if (dim(sub.df)[1]>1){
    ppa <- sub.df$PPA %>% sum(.)
    sub.df <- sub.df[1,]
    print("\nduplicate snp: " %&% snp.id %&% " : Locus.ID " %&% loc.id)
    chrom <- sub.df$CHR; pos <- sub.df$POS;
    coding <- sub.df$coding  
  } else{
    chrom <- sub.df$CHR; pos <- sub.df$POS;
    ppa <- sub.df$PPA; coding <- sub.df$coding    
  }
  if (coding==1){
    coding.vec <- handle_coding(loc.id,ppa) # islet, muscle, adipose, liver
  } else{
    coding.vec <- c(0,0,0,0)
  }
  annot.df <- filter(states.df,V1==chrom,V2<=pos,V3>=pos)
  strongenh.vec <- handle_annotation(annot.df,"strong_enhancer",ppa)
  weakenh.vec <- handle_annotation(annot.df,"weak_enhancer",ppa)
  genetrans.vec <- handle_annotation(annot.df,"gene_transcription",ppa)
  prom.vec <- handle_annotation(annot.df,"promoter",ppa)
  genenh.vec <- handle_annotation(annot.df,"genic_enhancer",ppa)
  if (weights==TRUE){
    c.df <- filter(weight.all.df,annotation=="coding")
    c.scores <- c(filter(c.df,tissue=="islet")$weight, filter(c.df,tissue=="muscle")$weight,
                  filter(c.df,tissue=="adipose")$weight,filter(c.df,tissue=="liver")$weight)
    coding.vec <- coding.vec  * c.scores

    se.df <- filter(weight.all.df,annotation=="strong.enhancers")
    se.scores <- c(filter(se.df,tissue=="islet")$weight, filter(se.df,tissue=="muscle")$weight,
                   filter(se.df,tissue=="adipose")$weight,filter(se.df,tissue=="liver")$weight)
    strongenh.vec <- strongenh.vec * se.scores

    we.df <- filter(weight.all.df,annotation=="weak.enhancers")
    we.scores <- c(filter(we.df,tissue=="islet")$weight, filter(we.df,tissue=="muscle")$weight,
                   filter(we.df,tissue=="adipose")$weight,filter(we.df,tissue=="liver")$weight)
    weakenh.vec <- weakenh.vec * we.scores

    gt.df <- filter(weight.all.df,annotation=="gene.transcription")
    gt.scores <- c(filter(gt.df,tissue=="islet")$weight, filter(gt.df,tissue=="muscle")$weight,
                   filter(gt.df,tissue=="adipose")$weight,filter(gt.df,tissue=="liver")$weight)
    genetrans.vec <- genetrans.vec * gt.scores

    pr.df <- filter(weight.all.df,annotation=="promoters")
    pr.scores <- c(filter(pr.df,tissue=="islet")$weight, filter(pr.df,tissue=="muscle")$weight,
                   filter(pr.df,tissue=="adipose")$weight, filter(pr.df,tissue=="liver")$weight)
    prom.vec <- prom.vec * pr.scores

    ge.df <- filter(weight.all.df,annotation=="genic.enhancer")
    ge.scores <- c(filter(ge.df,tissue=="islet")$weight,filter(ge.df,tissue=="muscle")$weight,
                   filter(ge.df,tissue=="adipose")$weight,filter(ge.df,tissue=="liver")$weight)
    genenh.vec <- genenh.vec * ge.scores
  }
  if (mode=="full"){
    score.vec <- coding.vec + strongenh.vec + weakenh.vec + genenh.vec + prom.vec + genetrans.vec
  } else if(mode=="coding+strongenhancers"){
    score.vec <- coding.vec + strongenh.vec
  } else{
    print("Need to enter acceptable mode")
  }
  if (sum(score.vec)>0){
    score.vec <- score.vec / (sum(score.vec))
  } else{
    score.vec <- c(0,0,0,0)
  }
  ppa.vec <- ppa * score.vec
  names(ppa.vec) <- c("islet","muscle","adipose","liver")
  out.df <- as.data.frame(t(ppa.vec))
  return(out.df)
}

divvy_ppa_loc <- function(loc.id,mode="full",weights=FALSE,scaled=FALSE){
  sub.df <- filter(cred.df,Locus.ID==loc.id)
  out.df <- c()
  snp.vec <- sub.df$SNPID %>% unique(.)
  pb <- txtProgressBar(min=0,max=length(snp.vec),style=3)
  for (i in 1:length(snp.vec)){
    setTxtProgressBar(pb,i)
    snp.id <- snp.vec[i]
    build.df <- divvy_ppa_snp(loc.id,snp.id,mode,weights)
    out.df <- rbind(out.df,build.df)
  }
  mat <- as.matrix(out.df)
  vec <- colSums(mat)
  cumppa <- sub.df$PPA %>% sum(.)
  other <- ifelse((cumppa-sum(vec))<=1&(cumppa-sum(vec))>=0,cumppa-sum(vec),0)
  vec <- c(vec,other)
  names(vec)[length(vec)] <- "other"
  if (scaled==TRUE){
    # Scale ppa's to add up to 1
    vec <- vec/sum(vec)
  }
  out.df <- as.data.frame(t(vec))
  out.df <- data.frame(Locus.ID=loc.id,out.df)
  return(out.df)
}

build_ppa_partition_df <- function(mode="full",weights=FALSE,scaled=FALSE){
  loc.ids <- cred.df$Locus.ID %>% unique(.)
  out.df <- c()
  pb <- txtProgressBar(min=0,max=length(loc.ids),style=3)
  for (i in 1:length(loc.ids)){
    setTxtProgressBar(pb,i)
    loc.id <- loc.ids[i]
    print(loc.id)
    build.df <- divvy_ppa_loc(loc.id,mode,weights,scaled)
    out.df <- rbind(out.df,build.df)
  }
  return(out.df)
}



#Generate and save tables

part.fullUU.df <- build_ppa_partition_df(mode="full",weights=FALSE,scaled=FALSE)


write.table(x=part.fullUU.df,file=out.dir%&%"tissue_ppa_divvy-full-unweighted-unscaled.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


