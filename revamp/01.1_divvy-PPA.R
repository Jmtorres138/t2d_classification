# Setup



"%&%" <- function(a,b) paste0(a,b)

library("tidyverse")
library("data.table")
library("Homo.sapiens")

serv.dir <- "/well/mccarthy/users/jason/"
proj.dir <- serv.dir %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "revamp/"
out.dir <- proj.dir %&% "revamp/analysis_files/"
res.dir <- work.dir %&% "genetic_credible_sets/" 
#input.dir <- proj.dir %&% "analysis_files/"

sym.ids <- unique(keys(Homo.sapiens, keytype = "SYMBOL"))
sym.df <- select(Homo.sapiens,key=sym.ids,keytype="SYMBOL",
                 columns=c("ENSEMBL"))

ess.dir <- proj.dir %&% "analysis_files/"
ess.df <- fread(ess.dir %&% "expression_specificity_scores.txt")

weight.all.df <- fread(out.dir %&%  "weight-enrich-all.txt")


cred.df <- fread(res.dir %&% "gencred.txt")
bed1.df <- fread(out.dir %&% "all_shared.bed")
bed2.df <- fread(out.dir %&% "shared_processed.txt")
shared.df <- rbind(bed1.df,bed2.df); rm(bed1.df); rm(bed2.df)
specific.df <- fread(out.dir %&% "specific_processed.txt")
states.df <- rbind(shared.df,specific.df)


# Divvy function

locid_to_symbol <- function(loc.id){
  return(filter(cred.df,CondID==loc.id)$symbol %>% as.character(.) %>% unique(.))
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
  if (dim(sub)[1] == 1){
    score.vec <- c(sub$islet.score,sub$muscle.score,sub$adipose.score,sub$liver.score)
    if (any(is.na(score.vec)==TRUE)){
      score.vec <- c(0.25,0.25,0.25,0.25)
    }
  } else{
    score.vec <- c(0.25,0.25,0.25,0.25)
  }
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



divvy_ppa_snp <- function(loc.id,snp.id,mode="full",weights=TRUE){
  # set weights=NULL if unweighted
  sub.df <- filter(cred.df,CondID==loc.id,SNPID==snp.id)
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
  genenh.vec <- handle_annotation(annot.df,"genic_enhancer",ppa)
  strongprom.vec <- handle_annotation(annot.df,"strong_promoter",ppa)
  weakprom.vec <- handle_annotation(annot.df,"weak_promoter",ppa)
  flankprom.vec <- handle_annotation(annot.df,"flank_promoter",ppa)
  stronggenetrans.vec <- handle_annotation(annot.df,"strong_gene_transcription",ppa)
  weakgenetrans.vec <- handle_annotation(annot.df,"weak_gene_transcription",ppa)
  
  if (weights==TRUE){
    c.df <- filter(weight.all.df,annotation=="coding")
    c.scores <- c(2^(filter(c.df,tissue=="Islet")$weight), 2^(filter(c.df,tissue=="Muscle")$weight),
                  2^(filter(c.df,tissue=="Adipose")$weight),2^(filter(c.df,tissue=="Liver")$weight))
    coding.vec <- coding.vec  * c.scores

    se.df <- filter(weight.all.df,annotation=="strong.enhancers")
    se.scores <- c(2^(filter(se.df,tissue=="Islet")$weight), 2^(filter(se.df,tissue=="Muscle")$weight),
                   2^(filter(se.df,tissue=="Adipose")$weight),2^(filter(se.df,tissue=="Liver")$weight))
    strongenh.vec <- strongenh.vec * se.scores

    we.df <- filter(weight.all.df,annotation=="weak.enhancers")
    we.scores <- c(2^(filter(we.df,tissue=="Islet")$weight), 2^(filter(we.df,tissue=="Muscle")$weight),
                   2^(filter(we.df,tissue=="Adipose")$weight),2^(filter(we.df,tissue=="Liver")$weight))
    weakenh.vec <- weakenh.vec * we.scores

    sgt.df <- filter(weight.all.df,annotation=="strong.gene.transcription")
    sgt.scores <- c(2^(filter(sgt.df,tissue=="Islet")$weight), 2^(filter(sgt.df,tissue=="Muscle")$weight),
                   2^(filter(sgt.df,tissue=="Adipose")$weight),2^(filter(sgt.df,tissue=="Liver")$weight))
    stronggenetrans.vec <- stronggenetrans.vec * sgt.scores
    
    wgt.df <- filter(weight.all.df,annotation=="weak.gene.transcription")
    wgt.scores <- c(2^(filter(wgt.df,tissue=="Islet")$weight), 2^(filter(wgt.df,tissue=="Muscle")$weight),
                   2^(filter(wgt.df,tissue=="Adipose")$weight),2^(filter(wgt.df,tissue=="Liver")$weight))
    weakgenetrans.vec <- weakgenetrans.vec * wgt.scores

    
    spr.df <- filter(weight.all.df,annotation=="strong.promoter")
    spr.scores <- c(2^(filter(spr.df,tissue=="Islet")$weight), 2^(filter(spr.df,tissue=="Muscle")$weight),
                   2^(filter(spr.df,tissue=="Adipose")$weight), 2^(filter(spr.df,tissue=="Liver")$weight))
    strongprom.vec <- strongprom.vec * spr.scores

    wpr.df <- filter(weight.all.df,annotation=="weak.promoter")
    wpr.scores <- c(2^(filter(wpr.df,tissue=="Islet")$weight), 2^(filter(wpr.df,tissue=="Muscle")$weight),
                    2^(filter(wpr.df,tissue=="Adipose")$weight), 2^(filter(wpr.df,tissue=="Liver")$weight))
    weakprom.vec <- weakprom.vec * wpr.scores
  
  
    fpr.df <- filter(weight.all.df,annotation=="flank.promoter")
    fpr.scores <- c(2^(filter(fpr.df,tissue=="Islet")$weight), 2^(filter(fpr.df,tissue=="Muscle")$weight),
                    2^(filter(fpr.df,tissue=="Adipose")$weight), 2^(filter(fpr.df,tissue=="Liver")$weight))
    flankprom.vec <- flankprom.vec * fpr.scores  
    
    
    ge.df <- filter(weight.all.df,annotation=="genic.enhancer")
    ge.scores <- c(2^(filter(ge.df,tissue=="Islet")$weight),2^(filter(ge.df,tissue=="Muscle")$weight),
                   2^(filter(ge.df,tissue=="Adipose")$weight),2^(filter(ge.df,tissue=="Liver")$weight))
    genenh.vec <- genenh.vec * ge.scores
  }
  if (mode=="full"){
    score.vec <- coding.vec + strongenh.vec + weakenh.vec + genenh.vec + 
      strongprom.vec + weakprom.vec + flankprom.vec + stronggenetrans.vec + weakgenetrans.vec
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

divvy_ppa_loc <- function(loc.id,mode="full",weights=TRUE,scaled=FALSE){
  sub.df <- filter(cred.df,CondID==loc.id)
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
  out.df$Locus.ID <- out.df$Locus.ID %>% as.character(.)
  return(out.df)
}

build_ppa_partition_df <- function(mode="full",weights=TRUE,scaled=FALSE){
  loc.ids <- cred.df$CondID %>% unique(.)
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

part.fullWU.df <- build_ppa_partition_df(mode="full",weights=TRUE,scaled=FALSE)



write.table(x=part.fullWU.df,file=out.dir%&%"tissue_ppa_divvy-full-weighted-unscaled.txt",
            sep="\t",quote=FALSE,row.names=FALSE)



