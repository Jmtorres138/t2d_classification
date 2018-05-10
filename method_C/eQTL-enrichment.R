set.seed(1)

"%&%" <- function(a,b) paste0(a,b)

library("tidyverse")
library("data.table")
library("Homo.sapiens")


serv.dir1 <- "/well/got2d/jason/"
serv.dir2 <- "/well/mccarthy/users/jason/"

proj.dir <- serv.dir2 %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "method_A/"

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

evaluate_threshold <- function(x,fcred.df,input.df,iter){
  # x is evaluated threshold
  classified <- map(input.df$Locus.ID,function(id){
    sub <- filter(input.df,Locus.ID==id) %>% dplyr::select(-one_of("Locus.ID","other")) %>%
      sort(.,decreasing=TRUE) %>% as.data.frame(.)
    tiss <- names(sub)[1]
    val <- sub[,1]
    tiss <- ifelse(val>=x,tiss,"unclassified")
    return(tiss)
  }) %>% as.character(.)
  islet.loci <- input.df$Locus.ID[classified=="islet"]
  liver.loci <- input.df$Locus.ID[classified=="liver"]
  adipose.loci <- input.df$Locus.ID[classified=="adipose"]
  muscle.loci <- input.df$Locus.ID[classified=="muscle"]
  islet.snps <- filter(fcred.df,Locus.ID %in% islet.loci)$SNPID
  liver.snps <- filter(fcred.df,Locus.ID %in% liver.loci)$SNPID
  adipose.snps <- filter(fcred.df,Locus.ID %in% adipose.loci)$SNPID
  muscle.snps <- filter(fcred.df,Locus.ID %in% muscle.loci)$SNPID
  df1 <- enrich_test(islet.snps,islet.eqtls,iter)
  df2 <- enrich_test(liver.snps,liv.eqtls,iter)
  df3 <- enrich_test(adipose.snps,adi.eqtls,iter)
  df4 <- enrich_test(muscle.snps,mus.eqtls,iter)
  df0 <- data.frame(tissue=c("islet","liver","adipose","muscle"))
  out.df <- rbind(df1,df2,df3,df4)
  out.df <- cbind(df0,out.df)
  out.df$tissue <- as.character(out.df$tissue)
  out.df$threshold <- rep(x,dim(out.df)[1])
  return(out.df)
}

evaluate_thresholds <- function(x.vec=seq(0,1,0.05),fcred.df,input.df,iter){
  out.df <- c()
  for (x in x.vec){
    print("Threshold:");print(x)
    build.df <- evaluate_threshold(x,fcred.df,input.df,iter)
    out.df <- rbind(out.df,build.df)
  }
  return(out.df)
}


fcred.dir <- serv.dir2 %&% "projects/t2d_classification/method_A/multi_results/"
fcred.df <- fread(fcred.dir %&% "results_func-cred-sets.txt")

input1.df <- fread(proj.dir %&% "method_A/analysis_files/tissue_ppa_divvy-full-weighted-scaled.txt")
thresh1.df <- evaluate_thresholds(x.vec=seq(0,1,0.05),fcred.df,input1.df,iter=1000)
write.table(x=thresh1.df,file=proj.dir%&%"method_A/analysis_files/select-thresh-eqtl_full-LocusScaled.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


input2.df <- fread(proj.dir %&% "method_A/analysis_files/tissue_ppa_divvy-full-weighted-unscaled.txt")
thresh2.df <- evaluate_thresholds(x.vec=seq(0,1,0.05),fcred.df,input2.df,iter=1000)
write.table(x=thresh2.df,file=proj.dir%&%"method_A/analysis_files/select-thresh-eqtl_full-unscaled.txt",
            sep="\t",quote=FALSE,row.names=FALSE)

#input2.df <- fread(proj.dir %&% "method_A/analysis_files/tissue_ppa_divvy-coding-strongEnhancers.txt")
#thresh2.df <- evaluate_thresholds(x.vec=seq(0,1,0.05),fcred.df,input2.df,iter=100)
#write.table(x=thresh2.df,file=proj.dir%&%"method_A/analysis_files/select-thresh-eqtl_cse.txt",
#            sep="\t",quote=FALSE,row.names=FALSE)
