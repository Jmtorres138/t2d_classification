# Setup



"%&%" <- function(a,b) paste0(a,b)

library("tidyverse")
library("data.table")
library("Homo.sapiens")

serv.dir <- "/well/mccarthy/users/jason/"
proj.dir <- serv.dir %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "method_C/"
out.dir <- proj.dir %&% "method_C/analysis_files/"
res.dir <- work.dir %&% "genetic_credible_sets/" 
input.dir <- proj.dir %&% "analysis_files/"

sym.ids <- unique(keys(Homo.sapiens, keytype = "SYMBOL"))
sym.df <- select(Homo.sapiens,key=sym.ids,keytype="SYMBOL",
                 columns=c("ENSEMBL"))


cred.df <- fread(work.dir %&% "genetic_credible_sets/gencred.txt")
class.df <- fread(work.dir %&% "analysis_files/classified-loci_weighted_renamed.txt")

map_condids <- function(df){
  Cond.ID <- map(1:dim(df)[1],function(i){
    filter(cred.df,Locus.ID==df$Locus.ID[i])$CondID %>% unique(.)
  }) %>% as.character(.)
  df$Locus.ID <- Cond.ID
  return(df)
}


atac.bed <- fread(serv.dir %&% "datasets/atac_islet_cells_kim/pancreas_cell_all.bed")


make_gr <- function(df){
  GRanges(seqnames = df$V1,IRanges(start=df$V2,end=df$V3))
}

cred.gr <- GRanges(seqnames=cred.df$CHR,IRanges(start=cred.df$POS,end = cred.df$POS))
alpha.gr <- filter(atac.bed,V4=="alpha_DOR") %>% make_gr(.)
beta.gr <- filter(atac.bed,V4=="beta_DOR") %>% make_gr(.)
duct.gr <- filter(atac.bed,V4=="duct_DOR") %>% make_gr(.)
endocrine.gr <- filter(atac.bed,V4=="endocrine_DOR") %>% make_gr(.)
exocrine.gr <- filter(atac.bed,V4=="exocrine_DOR") %>% make_gr(.)
islet.gr <- filter(atac.bed,V4=="islet_cell_atac") %>% make_gr(.)

cred.df$islet <- (cred.gr %over% islet.gr) %>% as.integer(.)
cred.df$alpha <- (cred.gr %over% alpha.gr) %>% as.integer(.)
cred.df$beta <- (cred.gr %over% beta.gr) %>% as.integer(.)
cred.df$duct <- (cred.gr %over% duct.gr) %>% as.integer(.)
cred.df$endocrine <- (cred.gr %over% endocrine.gr) %>% as.integer(.)
cred.df$exocrine <- (cred.gr %over% exocrine.gr) %>% as.integer(.)

annot.vec <- c("coding","islet", "alpha" ,"beta" ,"duct" ,"endocrine" ,"exocrine")

profile_locus <- function(id){
  sub <- filter(cred.df,CondID==id)
  vec <- c()
  for (a in annot.vec){
    member.vec <- dplyr::select(sub,one_of(a))[,1]
    ppa.vec <- dplyr::select(sub,one_of("PPA"))[,1]
    score <- (member.vec * ppa.vec) %>% sum(.)
    vec <- append(vec,score)
  }
  names(vec) <- annot.vec
  return(as.data.frame(t(vec)))
}


profile_all <- function(){
  id.vec <- cred.df$CondID %>% unique(.)
  out.df <- c() 
  pb <- txtProgressBar(min=0,max=length(id.vec),style = 3)
  for (i in 1:length(id.vec)){
    setTxtProgressBar(pb,i)
    id <- id.vec[i]
    build.df <- profile_locus(id)
    out.df <- rbind(out.df,build.df)
  }
  out.df$Locus.ID <- id.vec
  out.df$symbol <- map(out.df$Locus.ID,function(id){
    filter(cred.df,CondID==id)$symbol %>% unique(.)
  }) %>% as.character(.)
  out.df$maxppa <- map(out.df$Locus.ID,function(id){
    filter(cred.df,CondID==id)$PPA %>% max(.)
  }) %>% as.numeric(.)  
  out.df$assigned_20 <- map(out.df$Locus.ID,function(id){
    filter(class.df,Locus.ID==id)$assigned_20
  }) %>% as.character(.)
  return(out.df)
}


res.df <- profile_all()
isl.df <- filter(res.df,assigned_20=="islet",maxppa>0.50)

write.table(x=res.df,file=out.dir %&% "islet-profile-all.txt",sep="\t",quote=F,row.names=F)
write.table(x=isl.df,file=out.dir %&% "islet-profile-IsletMaxPPA50.txt",sep="\t",quote=F,row.names=F)

