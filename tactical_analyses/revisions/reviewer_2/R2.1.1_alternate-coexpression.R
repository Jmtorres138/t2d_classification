# Setup
set.seed(1)
"%&%" <- function(a,b) paste0(a,b)
#library("tidyverse")
library("dplyr")
library("purrr")
library("GenomicRanges")
library("viridis")
library("data.table")
library("rtracklayer")
library("Hmisc")
library("plyr")

serv.dir <- "/gpfs2/well/mccarthy/users/jason/"  #/home/jason/science/servers/FUSE5/"
proj.dir <- serv.dir %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "tactical_analyses/revisions/reviewer_2/"

gtf <- rtracklayer::import(serv.dir %&% "datasets/gencode.v30lift37.annotation.gtf.gz")
gtf.df<-as.data.frame(gtf) %>% filter(.,type=="gene",gene_type=="protein_coding")

cred.df <- fread(proj.dir %&% "tactical_analyses/genetic_credible_sets/gencred.txt")
cred.df$Locus.ID <- cred.df$CondID

weighted.df <- fread(proj.dir %&% "tactical_analyses/analysis_files/" %&%
                       "classified-loci_weighted_with-shared.txt")
weighted.df$assigned_00 <- gsub("other","unclassified",weighted.df$assigned_00)
weighted.df$assigned_20 <- gsub("other","unclassified",weighted.df$assigned_20)
weighted.df$assigned_50 <- gsub("other","unclassified",weighted.df$assigned_50)
weighted.df$assigned_80 <- gsub("other","unclassified",weighted.df$assigned_80)

annot.prof.df <- fread(proj.dir %&% "tactical_analyses/analysis_files/annotation-scores-signals.txt")
names(annot.prof.df)[1] <- "Locus.ID"
annot.prof.df$coding <- annot.prof.df$coding.Islet + annot.prof.df$coding.Adipose +
  annot.prof.df$coding.Muscle + annot.prof.df$coding.Liver
zero.vec <- filter(annot.prof.df,coding==0)$Locus.ID %>% unique(.)
p10.vec <- filter(annot.prof.df,coding<=0.1)$Locus.ID %>% unique(.) #358

tpm.dir <- proj.dir %&% "tactical_analyses/revisions/reviewer_1/analysis_files/"
isl.tpm.df <- fread(tpm.dir %&% "islet_expression_tpm_v2.txt")
liv.tpm.df <- fread(tpm.dir %&% "liver_expression_tpm.txt")
mus.tpm.df <- fread(tpm.dir %&% "muscle_expression_tpm.txt")
adisub.tpm.df <- fread(tpm.dir %&% "adipose_expression_tpm.txt")
adivis.tpm.df <- fread(tpm.dir %&% "adipose-visceral_expression_tpm.txt")
tpm.list <- list(isl.tpm.df,liv.tpm.df,mus.tpm.df,adisub.tpm.df,adivis.tpm.df)
tpm.names <- c("islet","liver","muscle","adipose-sub","adipose-visc")
# Functions



get_meansq_corr_fast <- function(gene.set,tpm.df){
  sub.df <- filter(tpm.df,GeneName%in%gene.set)
  sub.df <- sub.df[!duplicated(sub.df$GeneName),]
  sub.mat <- select(sub.df,-one_of(c("GeneID","GeneName"))) %>%
                    as.matrix(.) %>% t(.)
  res <- rcorr(sub.mat,type="spearman")$r
  corr.vec <- res[upper.tri(res,diag=FALSE)] %>% as.numeric(.) %>% na.omit(.)
  msc <- ((corr.vec)^2) %>% mean(.)
  return(msc)
}

get_num_sig_pos_corr <- function(gene.set,tpm.df){
  sub.df <- filter(tpm.df,GeneName%in%gene.set)
  sub.df <- sub.df[!duplicated(sub.df$GeneName),]
  sub.mat <- select(sub.df,-one_of(c("GeneID","GeneName"))) %>%
                    as.matrix(.) %>% t(.)
  rc.obj <- rcorr(sub.mat,type="spearman")
  res.r <- rc.obj$r
  res.p <- rc.obj$P
  upp.r <- res.r[upper.tri(res.r,diag=FALSE)]
  upp.p <- res.p[upper.tri(res.p,diag=FALSE)]
  count.val <- upp.r[upp.p < 0.05 & upp.r > 0] %>% na.omit(.) %>% length(.)
  return(count.val)
}


enrichment <- function(geneset, perms){
  out.df <- c()
  for (i in 1:length(tpm.list)){
    tpm.df <- tpm.list[[i]]
    tiss.name <- tpm.names[i]
    obs.msc <- get_meansq_corr_fast(geneset,tpm.df)
    null.vec <- c()
    pb <- txtProgressBar(1,perms,style=3)
    for (i in 1:perms){
      setTxtProgressBar(pb,i)
      nullset <- nullset <-  sample(nearest380,length(geneset)) #sample(tpm.df$GeneName,length(geneset))
      null.msc <- suppressWarnings(get_meansq_corr_fast(nullset,tpm.df))
      null.vec <- append(null.vec,null.msc)
    }
    count <- sum(null.vec >= obs.msc)
    perm.pvalue <- (sum(null.vec >= obs.msc) + 1)  / (perms + 1)
    log_pval <- -1*log(perm.pvalue,base=10)
    enrich.factor <- obs.msc / mean(null.vec)
    build.df <- data.frame("tissue"=tiss.name,"count"=count,
                         "perm.pvalue"=perm.pvalue,"log_pval"=log_pval,
                       "enrich.factor"=enrich.factor)
    out.df <- rbind(out.df,build.df)
  }
  return(out.df)
}

enrichment_v2 <- function(geneset, perms){
  out.df <- c()
  for (i in 1:length(tpm.list)){
    tpm.df <- tpm.list[[i]]
    tiss.name <- tpm.names[i]
    obs.nsp <- get_num_sig_pos_corr(geneset,tpm.df)
    null.vec <- c()
    pb <- txtProgressBar(1,perms,style=3)
    for (i in 1:perms){
      setTxtProgressBar(pb,i)
      nullset <- sample(nearest380,length(geneset)) #sample(tpm.df$GeneName,length(geneset)) #
      null.nsp <- suppressWarnings(get_num_sig_pos_corr(nullset,tpm.df))
      null.vec <- append(null.vec,null.nsp)
    }
    count <- sum(null.vec >= obs.nsp)
    perm.pvalue <- (sum(null.vec >= obs.nsp) + 1)  / (perms + 1)
    log_pval <- -1*log(perm.pvalue,base=10)
    enrich.factor <- obs.nsp / mean(null.vec)
    build.df <- data.frame("tissue"=tiss.name,"count"=count,
                         "perm.pvalue"=perm.pvalue,"log_pval"=log_pval,
                       "enrich.factor"=enrich.factor)
    out.df <- rbind(out.df,build.df)
  }
  return(out.df)
}

build_plot_df <- function(islet.genes,muscle.genes,
                          liver.genes,adipose.genes,
                          shared.genes,unclassified.genes,iter){
  if(length(islet.genes)>2){
    enrich.isl.df <- enrichment(islet.genes, iter)
  } else{
    enrich.isl.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                count=NA,log_pval=NA,enrich.factor=NA)
  }
  if(length(muscle.genes)>2){
    enrich.mus.df <- enrichment(muscle.genes, iter)
  } else{
    enrich.mus.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                count=NA,log_pval=NA,enrich.factor=NA)
  }
  if(length(adipose.genes)>2){
    enrich.adi.df <- enrichment(adipose.genes, iter)
  } else{
    enrich.adi.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                count=NA,log_pval=NA,enrich.factor=NA)
  }
  if(length(liver.genes)>2){
    enrich.liv.df <- enrichment(liver.genes, iter)
  } else{
    enrich.liv.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                count=NA,log_pval=NA,enrich.factor=NA)
  }
  if(length(shared.genes)>2){
    enrich.share.df <- enrichment(shared.genes, iter)
  } else{
    enrich.share.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                count=NA,log_pval=NA,enrich.factor=NA)
  }
#  if(length(unclassified.genes)>2){
  if(sum(unclassified.genes %in% gtf.df$gene_name)>2){
    enrich.unclass.df <- enrichment(unclassified.genes, iter)
  } else{
    enrich.unclass.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                  count=NA,log_pval=NA,enrich.factor=NA)
  }
  enrich.isl.df$geneset <- "islet"
  enrich.mus.df$geneset <- "muscle"
  enrich.liv.df$geneset <- "liver"
  enrich.adi.df$geneset <- "adipose"
  enrich.share.df$geneset <- "shared"
  enrich.unclass.df$geneset <- "unclassified"
  plot.df <- rbind(enrich.isl.df,enrich.mus.df,
                   enrich.adi.df,enrich.liv.df,enrich.share.df,enrich.unclass.df)
  return(plot.df)
}


build_plot_df_v2 <- function(islet.genes,muscle.genes,
                          liver.genes,adipose.genes,
                          shared.genes,unclassified.genes,iter){
  if(length(islet.genes)>2){
    enrich.isl.df <- enrichment_v2(islet.genes, iter)
  } else{
    enrich.isl.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                count=NA,log_pval=NA,enrich.factor=NA)
  }
  if(length(muscle.genes)>2){
    enrich.mus.df <- enrichment_v2(muscle.genes, iter)
  } else{
    enrich.mus.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                count=NA,log_pval=NA,enrich.factor=NA)
  }
  if(length(adipose.genes)>2){
    enrich.adi.df <- enrichment_v2(adipose.genes, iter)
  } else{
    enrich.adi.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                count=NA,log_pval=NA,enrich.factor=NA)
  }
  if(length(liver.genes)>2){
    enrich.liv.df <- enrichment_v2(liver.genes, iter)
  } else{
    enrich.liv.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                count=NA,log_pval=NA,enrich.factor=NA)
  }
  if(length(shared.genes)>2){
    enrich.share.df <- enrichment_v2(shared.genes, iter)
  } else{
    enrich.share.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                count=NA,log_pval=NA,enrich.factor=NA)
  }
#  if(length(unclassified.genes)>2){
  if(sum(unclassified.genes %in% gtf.df$gene_name)>2){
    enrich.unclass.df <- enrichment_v2(unclassified.genes, iter)
  } else{
    enrich.unclass.df <- data.frame(tissue=NA,perm.pvalue=NA,
                                  count=NA,log_pval=NA,enrich.factor=NA)
  }
  enrich.isl.df$geneset <- "islet"
  enrich.mus.df$geneset <- "muscle"
  enrich.liv.df$geneset <- "liver"
  enrich.adi.df$geneset <- "adipose"
  enrich.share.df$geneset <- "shared"
  enrich.unclass.df$geneset <- "unclassified"
  plot.df <- rbind(enrich.isl.df,enrich.mus.df,
                   enrich.adi.df,enrich.liv.df,enrich.share.df,enrich.unclass.df)
  return(plot.df)
}



find_specified_nearest_gene <- function(indexsnp,nearest.rank){
  vec <- strsplit(indexsnp,split=":")[[1]]
  pos <- vec[2] %>% as.integer(.)
  c <- vec[1]  #strsplit(vec[1],split="chr")[[1]][2]
  g.df <- filter(gtf.df,seqnames==c)
  g.df$absdist <- abs(g.df$start - pos )
  g.df <- arrange(g.df,absdist)
  g.df$gene_name[nearest.rank]
}

build_complete_df_extendedGenes <- function(group.df,iter,nearest.rank){
  thresh.vec <- c(0,0.2,0.5,0.8)
  out.df <- c()
  for (t in thresh.vec){
    print(t)
    var <- ifelse(t==0,"assigned_00",
                  ifelse(t==0.2,"assigned_20",
                         ifelse(t==0.5,"assigned_50",
                                ifelse(t=0.80,"assigned_80"))))
    sub <- dplyr::select(group.df,one_of(var,"Locus.ID"))
    sub$indexsnp <- map(sub$Locus.ID,function(id){
      filter(cred.df,CondID==id)$IndexSNP %>% unique(.)
    }) %>% as.character(.)
    names(sub)[1] <- "assigned"

    gene <- c()
    for (i in 1:dim(sub)[1]){
      gene <- append(gene,
                     find_specified_nearest_gene(sub$indexsnp[i],nearest.rank))
    }
    sub$gene <- gene
    islet.genes <- filter(sub,assigned=="islet")$gene %>% unique(.)
    muscle.genes <-  filter(sub,assigned=="muscle")$gene %>% unique(.)
    liver.genes <-  filter(sub,assigned=="liver")$gene %>% unique(.)
    adipose.genes <-  filter(sub,assigned=="adipose")$gene %>% unique(.)
    shared.genes <-  filter(sub,assigned=="shared")$gene %>% unique(.)
    unclassified.genes <-  filter(sub,assigned=="unclassified")$gene %>% unique(.)

    plot.df <- build_plot_df(islet.genes,muscle.genes,
                             liver.genes,adipose.genes,shared.genes,unclassified.genes,iter)
    plot.df$threshold <- t
    out.df <- rbind(out.df,plot.df)
  }
  return(out.df)
}


build_complete_df_extendedGenes_v2 <- function(group.df,iter,nearest.rank){
  thresh.vec <- c(0,0.2,0.5,0.8)
  out.df <- c()
  for (t in thresh.vec){
    print(t)
    var <- ifelse(t==0,"assigned_00",
                  ifelse(t==0.2,"assigned_20",
                         ifelse(t==0.5,"assigned_50",
                                ifelse(t=0.80,"assigned_80"))))
    sub <- dplyr::select(group.df,one_of(var,"Locus.ID"))
    sub$indexsnp <- map(sub$Locus.ID,function(id){
      filter(cred.df,CondID==id)$IndexSNP %>% unique(.)
    }) %>% as.character(.)
    names(sub)[1] <- "assigned"

    gene <- c()
    for (i in 1:dim(sub)[1]){
      gene <- append(gene,
                     find_specified_nearest_gene(sub$indexsnp[i],nearest.rank))
    }
    sub$gene <- gene
    islet.genes <- filter(sub,assigned=="islet")$gene %>% unique(.)
    muscle.genes <-  filter(sub,assigned=="muscle")$gene %>% unique(.)
    liver.genes <-  filter(sub,assigned=="liver")$gene %>% unique(.)
    adipose.genes <-  filter(sub,assigned=="adipose")$gene %>% unique(.)
    shared.genes <-  filter(sub,assigned=="shared")$gene %>% unique(.)
    unclassified.genes <-  filter(sub,assigned=="unclassified")$gene %>% unique(.)

    plot.df <- build_plot_df_v2(islet.genes,muscle.genes,
                             liver.genes,adipose.genes,shared.genes,unclassified.genes,iter)
    plot.df$threshold <- t
    out.df <- rbind(out.df,plot.df)
  }
  return(out.df)
}

full.df <- filter(weighted.df,Locus.ID %in% p10.vec)
full.df$indexsnp <- map(full.df$Locus.ID,function(id){
  filter(cred.df,CondID==id)$IndexSNP %>% unique(.)
}) %>% as.character(.)
nearest380 <- c()
for (i in 1:dim(full.df)[1]){
  nearest380 <- append(nearest380,
                 find_specified_nearest_gene(full.df$indexsnp[i],nearest.rank=1))
}



#df1st <- build_complete_df_extendedGenes(filter(weighted.df,Locus.ID %in% p10.vec),iter=10000,nearest.rank = 1)
#write.table(x=df1st,file=work.dir%&%"analysis_files/coexpress-msc_1stNearest_SPEARMAN.txt",
#            sep="\t",row.names=F,quote=F)
#
#df2nd <- build_complete_df_extendedGenes(filter(weighted.df,Locus.ID %in% p10.vec),iter=10000,nearest.rank = 2)
#write.table(x=df2nd,file=work.dir%&%"analysis_files/coexpress-msc_2ndNearest.txt",
#            sep="\t",row.names=F,quote=F)

#df3rd <- build_complete_df_extendedGenes(filter(weighted.df,Locus.ID %in% p10.vec),iter=10000,nearest.rank = 3)
#write.table(x=df3rd,file=work.dir%&%"analysis_files/coexpress-msc_3rdNearest.txt",
#            sep="\t",row.names=F,quote=F)

#print("Number of positively correlated, signficant gene pairs")

#df1stB <- build_complete_df_extendedGenes_v2(filter(weighted.df,Locus.ID %in% p10.vec),iter=10000,nearest.rank = 1)
#write.table(x=df1stB,file=work.dir%&%"analysis_files/coexpress-nsp_1stNearest_SPEARMAN.txt",
#            sep="\t",row.names=F,quote=F)

#df2ndB <- build_complete_df_extendedGenes_v2(filter(weighted.df,Locus.ID %in% p10.vec),iter=10000,nearest.rank = 2)
#write.table(x=df2ndB,file=work.dir%&%"analysis_files/coexpress-nsp_2ndNearest.txt",
#            sep="\t",row.names=F,quote=F)

#df3rdB <- build_complete_df_extendedGenes_v2(filter(weighted.df,Locus.ID %in% p10.vec),iter=10000,nearest.rank = 3)
#write.table(x=df3rdB,file=work.dir%&%"analysis_files/coexpress-nsp_3rdNearest.txt",
#            sep="\t",row.names=F,quote=F)

#print("Null set 380, switched nullset line in enrichment functions")

df1stC <- build_complete_df_extendedGenes(filter(weighted.df,Locus.ID %in% p10.vec),iter=10000,nearest.rank = 1)
write.table(x=df1stC,file=work.dir%&%"analysis_files/coexpress-msc_1stNearest_null380_SPEARMAN.txt",
            sep="\t",row.names=F,quote=F)

df1stD <- build_complete_df_extendedGenes_v2(filter(weighted.df,Locus.ID %in% p10.vec),iter=10000,nearest.rank = 1)
write.table(x=df1stD,file=work.dir%&%"analysis_files/coexpress-nsp_1stNearest_null380_SPEARMAN.txt",
            sep="\t",row.names=F,quote=F)
