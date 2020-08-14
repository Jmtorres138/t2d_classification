
# Setup

#Script for performing eQTL enrichment test using two different approaches for generating null snp sets: 1) SNPsnap using CEU 1K Genomes reference and 2) MAF-based sampling from HRC-imputed GWAS data frame


"%&%" <- function(a,b) paste0(a,b)

library("rlang")
library("wrapr")
library("dplyr")
library("purrr")
library("data.table")

serv.dir1 <- "/well/got2d/jason/" ## "/home/jason/science/servers/FUSE/" #
serv.dir2 <- "/well/mccarthy/users/jason/" # "/home/jason/science/servers/FUSE5/" #
proj.dir <- serv.dir2 %&% "projects/t2d_classification/"
work.dir0 <- proj.dir %&% "tactical_analyses/"
work.dir <- proj.dir %&% "tactical_analyses/revisions/reviewer_1/"
out.dir <- work.dir %&% "enrichment_files/eqtls/"

toa.df <- fread(work.dir0 %&% "analysis_files/classified-loci_weighted_with-shared.txt")
cred.df <- fread(work.dir0 %&% "genetic_credible_sets/gencred.txt")

snpsnap.df <- fread(out.dir %&% "SNPsnap_oneK_5_20_20_50_rsq5/matched_snps.txt") # Note that only 359/380 (94%) signals were matched with SNPSNAP


isl.spec <- fread(out.dir %&% "islet-specific-esnps.txt",sep="\t",header=FALSE)$V1
mus.spec <- fread(out.dir %&% "muscle-specific-esnps.txt",sep="\t",header=FALSE)$V1
liv.spec <- fread(out.dir %&% "liver-specific-esnps.txt",sep="\t",header=FALSE)$V1
adi.visc.spec <- fread(out.dir %&% "adipose-visc-specific-esnps.txt",sep="\t",header=FALSE)$V1
adi.sub.spec <- fread(out.dir %&% "adipose-sub-specific-esnps.txt",sep="\t",header=FALSE)$V1
adi.union.spec <- fread(out.dir %&% "adipose-union-specific-esnps.txt",sep="\t",header=FALSE)$V1


stringToQuoser <- function(varName) {
  wrapr::let(c(VARNAME = varName), quo(VARNAME))
}

extract_query_vec <- function(tissue,threshold){
  # threshold must be: "00","20","50", or "80"
  cname <- ("assigned_" %&% threshold) %>% stringToQuoser(.)
  id.vec <- dplyr::filter(toa.df, !!cname==tissue)$Locus.ID
  if (length(id.vec)>0){
    query.vec <- map(1:length(id.vec),function(i){
      id <- id.vec[i]
      snp <- (filter(cred.df,CondID==id) %>% arrange(.,desc(PPA)))[1,]$IndexSNP
      snp <- strsplit(snp,split="chr")[[1]][2]
    }) %>% as.character(.)
  } else{
    query.vec <- c()
  }
  return(query.vec)
}

get_overlap <- function(query.vec,eqtl.vec){
  sum(query.vec %in% eqtl.vec)
}



enrich_test_snpsnap <- function(query.vec, eqtl.vec){
  if (length(query.vec)==0){
    out.df <- data.frame(observed=NA,enrichment=NA,pvalue=NA,num_sigs=NA,
                         num_missing=NA,stringsAsFactors = FALSE)
  } else{
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
  return(out.df)
}


get_build_df <- function(tissue,threshold){
  query.vec <- extract_query_vec(tissue,threshold)
  eqtl.list <- list(isl.spec,mus.spec,liv.spec,adi.visc.spec,adi.sub.spec,adi.union.spec)
  eqtl.names <- c("islet","muscle","liver","adipose.visc","adipose.sub","adipose.both")
  out.df <- c()
  for (i in 1:length(eqtl.list)){
    eqtl.vec = eqtl.list[[i]]
    ename = eqtl.names[i]
    print("eQTL set: " %&% ename)
    df2 <- enrich_test_snpsnap(query.vec,eqtl.vec)
    names(df2) <- c("observed_snpsnap","enrichment_snpsnap","pvalue_snpsnap","num_sigs","num_missing")
    df <- data.frame("tissue_toa"=tissue,tissue_eqtl=ename,"threshold"=threshold,stringsAsFactors=F)
    build.df <- cbind(df,df2)
    out.df <- rbind(out.df,build.df)
  }
  return(out.df)
}

build_enrich_df <- function(){
  out.df <- c()
  tiss.vec <- c("islet","liver","muscle","adipose","shared","unclassified")
  thresh.vec <- c("00","20","50","80")
  out.df <- c()
  for (tiss in tiss.vec){
    for (thresh in thresh.vec){
      print("TOA: " %&% tiss %&% "; Threshold: " %&% thresh)
      build.df <- get_build_df(tiss,thresh)
      print(build.df)
      out.df <- rbind(out.df,build.df)
    }
  }
  return(out.df)
}


# RUN
enrich.df <- build_enrich_df()
write.table(x=enrich.df,file=out.dir %&% "eqtl_enrichment.txt",sep="\t",quote=F,row.names=F)
