set.seed(2)

"%&%" <- function(a,b) paste0(a,b)

library("tidyverse")
library("data.table")
library("Homo.sapiens")


serv.dir1 <- "/well/got2d/jason/" # 
#serv.dir1 <- "/home/jason/science/servers/FUSE/"
serv.dir2 <- "/well/mccarthy/users/jason/" #
#serv.dir2 <- "/home/jason/science/servers/FUSE5/"

proj.dir <- serv.dir2 %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "revamp/"

gwas.dir <- serv.dir1 %&% "reference/gwas/diamante-ukbb_hrc/"
gwas.df <- fread(gwas.dir %&% "ukbb_diamante-euro.bed") # will take a while to load, 0.681 GB file
names(gwas.df) <- c("CHR","POS0","POS","SNPID","Freq","Z","PVAL","NCASE","NCONTROL")

pb <- txtProgressBar(min=0,max=dim(gwas.df)[1],style=3)
gwas.df$MAF <- map(1:dim(gwas.df)[1],function(i){
  setTxtProgressBar(pb,i)
  freq <- gwas.df$Freq[i]
  return(min(freq,1-freq))
}) %>% as.numeric(.)

#g.df <- dplyr::select(gwas.df,SNPID,MAF); g.df <- g.df[!duplicated(g.df),]

rare.snps <- filter(gwas.df, MAF < 0.005)$SNPID
low.freq.snps <- filter(gwas.df, MAF >= 0.005 & MAF < 0.05)$SNPID
high.freq.snps <- filter(gwas.df, MAF >= 0.05)$SNPID


gtex.dir <- serv.dir2 %&% "datasets/GTEx/v7/eqtl/GTEx_Analysis_v7_eQTL/"

liv.df <- fread("cat " %&% gtex.dir %&% "Liver.v7.signif_variant_gene_pairs.txt.gz" %&% " | zmore")
mus.df <- fread("cat " %&% gtex.dir %&% "Muscle_Skeletal.v7.signif_variant_gene_pairs.txt.gz" %&% " | zmore")
adi.df <- fread("cat " %&% gtex.dir %&% "Adipose_Subcutaneous.v7.signif_variant_gene_pairs.txt.gz" %&% " | zmore")

pb <- txtProgressBar(min=0,max=dim(liv.df)[1],style=3)
liv.eqtls <- map(1:length(liv.df$variant_id),function(i){
  setTxtProgressBar(pb,i)
  vec <- liv.df$variant_id[i] %>% strsplit(.,split="_")
  id <- "chr" %&% vec[[1]][1] %&% ":" %&% vec[[1]][2]
}) %>% as.character(.)

pb <- txtProgressBar(min=0,max=dim(mus.df)[1],style=3)
mus.eqtls <- map(1:length(mus.df$variant_id),function(i){
  setTxtProgressBar(pb,i)
  vec <- mus.df$variant_id[i] %>% strsplit(.,split="_")
  id <- "chr" %&% vec[[1]][1] %&% ":" %&% vec[[1]][2]
}) %>% as.character(.)

pb <- txtProgressBar(min=0,max=dim(adi.df)[1],style=3)
adi.eqtls <- map(1:length(adi.df$variant_id),function(i){
  setTxtProgressBar(pb,i)
  vec <- adi.df$variant_id[i] %>% strsplit(.,split="_")
  id <- "chr" %&% vec[[1]][1] %&% ":" %&% vec[[1]][2]
}) %>% as.character(.)


islet.dir <- serv.dir1 %&% "reference/islet/eqtls/oxford/nominal_pass/output/"
islet.df <- fread("cat " %&% islet.dir %&% "snp_keyfile_fdr05.txt.gz" %&% " | zmore")
pb <- txtProgressBar(min=0,max=dim(islet.df)[1],style=3)
islet.eqtls <- map(1:dim(islet.df)[1],function(i){
  setTxtProgressBar(pb,i)
  vec <- c(islet.df$CHR[i],islet.df$POS[i])
  id <- vec[1] %&% ":" %&% vec[2]
}) %>% as.character(.)


isl.spec <- islet.eqtls[!(islet.eqtls %in% c(mus.eqtls,liv.eqtls,adi.eqtls))]
mus.spec <- mus.eqtls[!(mus.eqtls %in% c(islet.eqtls,liv.eqtls,adi.eqtls))]
adi.spec <- adi.eqtls[!(adi.eqtls %in% c(islet.eqtls,liv.eqtls,mus.eqtls))]
liv.spec <- liv.eqtls[!(liv.eqtls %in% c(islet.eqtls,adi.eqtls,mus.eqtls))]

# min number, liver, n=37391

min.val <- min(c(length(isl.spec),length(mus.spec),
                 length(adi.spec),length(liv.spec)))

islet.eqtls <- sample(x=isl.spec,replace=FALSE,size=min.val)
mus.eqtls <- sample(x=mus.spec,replace=FALSE,size=min.val)
liv.eqtls <- sample(x=liv.spec,replace=FALSE,size=min.val)
adi.eqtls <- sample(x=adi.spec,replace=FALSE,size=min.val)


get_overlap <- function(query.vec,eqtl.vec){
  sum(query.vec %in% eqtl.vec)
}

get_maf_matched_null_vec <- function(query.vec){
  sub.df <- filter(gwas.df,SNPID%in%query.vec)
  non.dup.df <- sub.df[!duplicated(sub.df$SNPID),]  # dup.df <- sub.df[duplicated(sub.df$SNPID),]
  maf.vec <- non.dup.df$MAF
  rare.count <- sum(maf.vec < 0.005)
  low.freq.count <- sum(maf.vec >= 0.005 & maf.vec < 0.05) # MAF >= 0.005 & MAF < 0.05
  high.freq.count <- sum(maf.vec >= 0.05) # MAF >= 0.05
  rare.samp <- sample(rare.snps,size=rare.count,replace=FALSE)
  low.freq.samp <- sample(low.freq.snps,size=low.freq.count,replace=FALSE)
  high.freq.samp <- sample(high.freq.snps,size=high.freq.count,replace=FALSE)
  null.vec <- c(rare.samp,low.freq.samp,high.freq.samp)
  return(null.vec)
}

enrich_test <- function(query.vec,eqtl.vec,iter){
  obs <- get_overlap(query.vec,eqtl.vec)
  null.counts <- c()
  pb <- txtProgressBar(min=0,max=iter,style=3)
  for (i in 1:iter){
    setTxtProgressBar(pb,i)
    null.vec <- get_maf_matched_null_vec(query.vec)
    nullcount <- get_overlap(null.vec,eqtl.vec)
    null.counts <- append(null.counts,nullcount)
  }
  pval <- (sum(null.counts >= obs) + 1)/ (iter + 1)
  fac <- obs / mean(null.counts)
  out.df <- data.frame(observed=obs,enrichment=fac,pvalue=pval,stringsAsFactors = FALSE)
}

enrich_test_across_tissue <- function(query.vec,iter){
  print("islets..."); df1 <- enrich_test(query.vec,islet.eqtls,iter)
  print("liver..."); df2 <- enrich_test(query.vec,liv.eqtls,iter)
  print("adipose..."); df3 <- enrich_test(query.vec,adi.eqtls,iter)
  print("muscle..."); df4 <- enrich_test(query.vec,mus.eqtls,iter)
  out.df <- rbind(tissue=c("islet","liver","adipose","mucle"),df1,df2,df3,df4)
  out.df$tissue <- as.character(out.df$tissue)
  return(out.df)
}

evaluate_threshold <- function(input.df,classified,iter,threshold){
  islet.loci <- input.df$Locus.ID[classified=="islet"]
  liver.loci <- input.df$Locus.ID[classified=="liver"]
  adipose.loci <- input.df$Locus.ID[classified=="adipose"]
  muscle.loci <- input.df$Locus.ID[classified=="muscle"]
  islet.snps <- filter(cred.df,Locus.ID %in% islet.loci)$SNPID
  liver.snps <- filter(cred.df,Locus.ID %in% liver.loci)$SNPID
  adipose.snps <- filter(cred.df,Locus.ID %in% adipose.loci)$SNPID
  muscle.snps <- filter(cred.df,Locus.ID %in% muscle.loci)$SNPID
  
  df1a <- enrich_test(islet.snps,islet.eqtls,iter)
  df1b <- enrich_test(islet.snps,liv.eqtls,iter)
  df1c <- enrich_test(islet.snps,adi.eqtls,iter)
  df1d <- enrich_test(islet.snps,mus.eqtls,iter)
  
  df2a <- enrich_test(liver.snps,islet.eqtls,iter)
  df2b <- enrich_test(liver.snps,liv.eqtls,iter)
  df2c <- enrich_test(liver.snps,adi.eqtls,iter)
  df2d <- enrich_test(liver.snps,mus.eqtls,iter)
  
  df3a <- enrich_test(adipose.snps,islet.eqtls,iter)
  df3b <- enrich_test(adipose.snps,liv.eqtls,iter)
  df3c <- enrich_test(adipose.snps,adi.eqtls,iter)
  df3d <- enrich_test(adipose.snps,mus.eqtls,iter)
  
  df4a <- enrich_test(muscle.snps,islet.eqtls,iter)
  df4b <- enrich_test(muscle.snps,liv.eqtls,iter)
  df4c <- enrich_test(muscle.snps,adi.eqtls,iter)
  df4d <- enrich_test(muscle.snps,mus.eqtls,iter)
  
  
  df0 <- data.frame(tissue=c(rep("islet",4),rep("liver",4),rep("adipose",4),rep("muscle",4)),
                    eQTL=rep(c("islet","liver","adipose","muscle"),4))
  out.df <- rbind(df1a,df1b,df1c,df1d,df2a,df2b,df2c,df2d,df3a,df3b,df3c,df3d,df4a,df4b,df4c,df4d)
  out.df <- cbind(df0,out.df)
  out.df$tissue <- as.character(out.df$tissue)
  out.df$threshold <- rep(threshold,dim(out.df)[1])
  return(out.df)
}

evaluate_thresholds <- function(input.df,iter){
  print("Threshold:");print("0.00")
  df00 <- evaluate_threshold(input.df,classified=input.df$assigned_00,iter,"0.00")
  print("Threshold:");print("0.20")
  df20 <- evaluate_threshold(input.df,classified=input.df$assigned_20,iter,"0.20")
  print("Threshold:");print("0.50")
  df50 <- evaluate_threshold(input.df,classified=input.df$assigned_50,iter,"0.50")
  print("Threshold:");print("0.80")
  df80 <- evaluate_threshold(input.df,classified=input.df$assigned_80,iter,"0.80")
  out.df <- rbind(df00,df20,df50,df80)
  return(out.df)
}



cred.dir <- serv.dir2 %&% "projects/t2d_classification/revamp/genetic_credible_sets/"
cred.df <- fread(cred.dir %&% "gencred.txt") %>% filter(.,PPA>=0.01)
cred.df$Locus.ID <- cred.df$CondID
# PPA 1% filter reduces number of query SNPs from 126K to 5,010 SNPs 

#input1.df <- fread(proj.dir %&% "revamp/analysis_files/classified-loci_unweighted.txt")
#thresh1.df <- evaluate_thresholds(input1.df,iter=1000)
#write.table(x=thresh1.df,file=proj.dir%&%"revamp/analysis_files/downsampled_eqtl-specific-validation_unweighted.txt",
#            sep="\t",quote=FALSE,row.names=FALSE)

input2.df <- fread(proj.dir %&% "revamp/analysis_files/classified-loci_weighted.txt")
thresh2.df <- evaluate_thresholds(input2.df,iter=1000)
write.table(x=thresh2.df,file=proj.dir%&%"revamp/analysis_files/downsampled_eqtl-specific-validation_weighted-seed2.txt",
            sep="\t",quote=FALSE,row.names=FALSE)

