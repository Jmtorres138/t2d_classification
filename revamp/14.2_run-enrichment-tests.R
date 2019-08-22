
# Setup

#Script for performing eQTL enrichment test using two different approaches for generating null snp sets: 1) SNPsnap using CEU 1K Genomes reference and 2) MAF-based sampling from HRC-imputed GWAS data frame


"%&%" <- function(a,b) paste0(a,b)

library("rlang")
library("dplyr")
library("purrr")
library("wrapr")
library("data.table")

serv.dir1 <- "/well/got2d/jason/"
serv.dir2 <- "/well/mccarthy/users/jason/"
proj.dir <- serv.dir2 %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "revamp/"
out.dir <- work.dir %&% "enrichment_files/eqtls/"

toa.df <- fread(work.dir %&% "analysis_files/classified-loci_weighted_with-shared.txt")
cred.df <- fread(work.dir %&% "genetic_credible_sets/gencred.txt")

snpsnap.df <- fread(out.dir %&% "SNPsnap_oneK_5_20_20_50_rsq5/matched_snps.txt") # Note that only 359/380 (94%) signals were matched with SNPSNAP


isl.spec <- fread(out.dir %&% "islet-specific-esnps.txt",header=F)$V1
mus.spec <- fread(out.dir %&% "muscle-specific-esnps.txt",header=F)$V1
adi.spec <- fread(out.dir %&% "adipose-specific-esnps.txt",header=F)$V1
liv.spec <- fread(out.dir %&% "liver-specific-esnps.txt",header=F)$V1
insp.exon.spec <- fread(out.dir %&%"inspExon-specific-esnps.txt",header=F)$V1
insp.gene.spec <- fread(out.dir %&% "inspGene-specific-esnps.txt",header=F)$V1
mvdb.spec <- fread(out.dir %&% "mvdbExon-specific-esnps.txt",header=F)$V1

#rare.snps <- fread("cat " %&% out.dir %&% "bin_rare.txt.gz" %&% " | zmore",header=F)$V1
#low.snps <- fread("cat " %&% out.dir %&% "bin_low.txt.gz"  %&% " | zmore", header=F)$V1
#high.snps <- fread("cat " %&% out.dir %&% "bin_high.txt.gz" %&% " | zmore",header=F)$V1

maf.df <- fread(out.dir %&% "leadSNP_maf.txt")
names(maf.df) <- c("SNP","MAF")
maf.df$SNPID <- map(1:dim(maf.df)[1],function(i){
  strsplit(maf.df$SNP[i],split="chr")[[1]][2]
}) %>% as.character(.)


# Enrichment functions


stringToQuoser <- function(varName) {
  wrapr::let(c(VARNAME = varName), quo(VARNAME))
}

extract_query_vec <- function(tissue,threshold){
  # threshold must be: "00","20","50", or "80"
  cname <- ("assigned_" %&% threshold) %>% stringToQuoser(.)
  id.vec <- dplyr::filter(toa.df, !!cname==tissue)$Locus.ID
  query.vec <- map(1:length(id.vec),function(i){
    id <- id.vec[i]
    snp <- (filter(cred.df,CondID==id) %>% arrange(.,desc(PPA)))[1,]$IndexSNP
    snp <- strsplit(snp,split="chr")[[1]][2]
  }) %>% as.character(.)
  return(query.vec)
}

get_overlap <- function(query.vec,eqtl.vec){
  sum(query.vec %in% eqtl.vec)
}

get_maf_matched_null_vec <- function(query.vec){
  sub.df <- filter(maf.df,SNPID%in%query.vec)
  non.dup.df <- sub.df[!duplicated(sub.df$SNPID),]
  maf.vec <- non.dup.df$MAF
  rare.count <- sum(maf.vec < 0.005)
  low.freq.count <- sum(maf.vec >= 0.005 & maf.vec < 0.05) # MAF >= 0.005 & MAF < 0.05
  high.freq.count <- sum(maf.vec >= 0.05) # MAF >= 0.05
  rare.samp <- sample(rare.snps,size=rare.count,replace=FALSE)
  low.freq.samp <- sample(low.snps,size=low.freq.count,replace=FALSE)
  high.freq.samp <- sample(high.snps,size=high.freq.count,replace=FALSE)
  null.vec <- c(rare.samp,low.freq.samp,high.freq.samp)
  return(null.vec)
}


enrich_test_snpsnap <- function(query.vec, eqtl.vec){
  sub.df <- filter(snpsnap.df,Input_SNP %in% query.vec)
  num.sigs <- length(query.vec); num.miss <- num.sigs - dim(sub.df)[1]
  obs <- get_overlap(query.vec,eqtl.vec)
  null.counts <- c()
  pb <- txtProgressBar(min=2,max=(dim(sub.df)[2]),style=3)
  for (i in 2:(dim(sub.df)[2])){
    setTxtProgressBar(pb,i)
    null.vec <- sub.df[,i]
    nullcount <- get_overlap(null.vec,eqtl.vec)
    null.counts <- append(null.counts,nullcount)
  }
  pval <- (sum(null.counts >= obs) + 1)/ ((dim(sub.df)[2]-1) + 1)
  fac <- obs / mean(null.counts)
  out.df <- data.frame(observed=obs,enrichment=fac,pvalue=pval,num_sigs=num.sigs,
                       num_missing=num.miss,stringsAsFactors = FALSE)
}

enrich_test_sampGWAS <- function(query.vec,eqtl.vec,iter=1000){
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



get_build_df <- function(tissue,threshold){
  query.vec <- extract_query_vec(tissue,threshold)
  eqtl.list <- list(isl.spec,mus.spec,adi.spec,liv.spec,insp.exon.spec,
                    insp.gene.spec,mvdb.spec)
  eqtl.names <- c("islet","muscle","adipose","liver",
                  "inspire.exon","inspire.gene","mvdb.exon")
  out.df <- c()
  for (i in 1:length(eqtl.list)){
    eqtl.vec = eqtl.list[[i]]
    ename = eqtl.names[i]
    print(eqtl.vec)
    print(ename)
    #df1 <- enrich_test_sampGWAS(query.vec,eqtl.vec,iter=1000)
    #df2 <- enrich_test_snpsnap(query.vec,eqtl.vec)
    #names(df2) <- c("observed_snpsnap","enrichment_snpsnap","pvalue_snpsnap","num_sigs","num_missing")
    #df <- data.frame("tissue_toa"=tissue,tissue_eqtl=ename,"threshold"=threshold,stringsAsFactors=F)
    #build.df <- cbind(df,df1,df2)
    #out.df <- rbind(out.df,build.df)
  }
  #return(out.df)
}


# RUN

test <- get_build_df("liver","00")
test
str(test)
