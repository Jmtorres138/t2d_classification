# Setup 



"%&%" <- function(a,b) paste0(a,b)
library("tidyverse")
library("GenomicRanges")
library("viridis")
library("data.table")

library("plyr")

serv.dir <- "/gpfs2/well/mccarthy/users/jason/"#/home/jason/science/servers/FUSE5/"
proj.dir <- serv.dir %&% "projects/t2d_classification/"

work.dir <- proj.dir %&% "method_A/"
analysis.dir <- work.dir %&% "analysis_files/"
plot.dir <- work.dir %&% "plots/"

keep.df <- fread(proj.dir %&% "input_data/" %&% "380.locus.ID.txt")
keep.vec <- keep.df$Locus.ID_inCREDS

block.df <- fread(work.dir %&% "multi_results/results_blocks.txt")




#Read data files 



block.df <- fread(work.dir %&%"multi_results/results_blocks.txt",sep="\t")
fcred.df <- fread(work.dir %&%"multi_results/results_func-cred-sets.txt",sep="\t")






TPM.merged.ranks <- readRDS(proj.dir%&%"input_data/TPM.merged.ranks.rds")
TPMs.tissues.rank.matrix <- readRDS(proj.dir%&%"input_data/TPMs.tissues.rank.matrix.rds")




# Functions 






enrichment <- function(geneset, perms){
  #geneset_ens <- TPM.merged.ranks[TPM.merged.ranks$Symbol %in% geneset, "Gene"]
  #names(geneset_ens) <- geneset
  
  sub.df <- filter(TPM.merged.ranks,Symbol%in%geneset)
  geneset_ens <- sub.df$Gene
  names(geneset_ens) <- sub.df$Symbol
  
  # Calculate average expression to make null and random distributions
  TPM.merged.ranks_geneset <- TPM.merged.ranks[TPM.merged.ranks$Symbol %in% geneset, "ranks2"]
  names(TPM.merged.ranks_geneset) <- TPM.merged.ranks[TPM.merged.ranks$Symbol %in% geneset, "Gene"]
  #h(TPM.merged.ranks_geneset)
  
  # create null bg genes for creating null gene_sets
  random.bg.ens <- sapply(TPM.merged.ranks_geneset, function(x){
    up <- x-150
    down <- x+150
    sample_from <- TPM.merged.ranks[up:down, "Gene"]
    sample_from2 <- setdiff(sample_from, geneset_ens)
  }, simplify=FALSE)
  
  # Create random gene sets using ranks + and - 100 
  random.geneList.ens <- list()
  
  for (i in 1:perms){
    gene_set <- lapply(random.bg.ens, function(x){
      rand <- sample(x, 1)
    })
    random.geneList.ens[[i]] <- unlist(gene_set)
  }
  
  geneset.list.test.ens <- list()
  
  pb <- txtProgressBar(min = 1, max = perms, style = 3)
  for (i in 1:perms){
    setTxtProgressBar(pb, i)
    gene_set <- random.geneList.ens[[i]]
    gene_set <- gene_set[!is.na(gene_set)] # remove NAs
    lista <- lapply(TPMs.tissues.rank.matrix, function(tissue){
      tis <- tissue
      df2 <- tis[match(gene_set, rownames(tis)), ]
      df3 <- apply(df2[,-c(1:2)], 1, function(x){sum(x, na.rm=TRUE)})
    })
    geneset.list.test.ens[[i]] <- lista
  }
  close(pb)
  
  geneset.list.test.ens <- lapply(geneset.list.test.ens, function(x)ldply(x))
  geneset.list.test.ens.sum <- lapply(geneset.list.test.ens, function(x) apply(x[,-1], 1, sum) )
  
  geneset.true <- lapply(TPMs.tissues.rank.matrix, function(tissue){
    tis <- tissue
    data <- tis[match(geneset_ens, rownames(tis)), ]
    data <- apply(data[,-c(1:2)], 1, sum)
  })
  
  geneset.true.sum <- lapply(geneset.true, function(x){sum(x, na.rm=T)})
  
  res.lista = list()
  
  for (i in 1:54){
    perm.value = c()
    for (j in 1:perms){
      perm = geneset.list.test.ens.sum[[j]][i]
      perm.value = c(perm.value, perm)
    }
    res.lista[[i]] = perm.value
  }
  
  names(res.lista) = names(geneset.true)
  zim <- data.frame(tissue = names(geneset.true), perm.pvalue = -9, count = -9, log_pval = -9,
                    enrich.factor=-9)
  for ( i in 1:length(res.lista)){
    zim[i,2] <- (sum(res.lista[[i]] <= geneset.true.sum[[i]])+1)/(length(res.lista[[i]])+1)
    zim[i,3] <- (sum(geneset.true.sum[[i]] >=res.lista[[i]])+1)
    zim[i,4] <- -log10(zim[i,2])
    zim[i,5] <- (res.lista[[i]] %>% mean(.)) / geneset.true.sum[[i]] # May need to ammend
  }
  zim <- zim[order(zim$log_pval, decreasing = TRUE),]
  #ht(zim)
  return(zim)
}

build_plot_df <- function(islet.genes,muscle.genes,
                          liver.genes,adipose.genes,iter){
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
  enrich.isl.df$geneset <- "islet"
  enrich.mus.df$geneset <- "muscle"
  enrich.liv.df$geneset <- "liver"
  enrich.adi.df$geneset <- "adipose"
  plot.df <- rbind(enrich.isl.df,enrich.mus.df,
                   enrich.adi.df,enrich.liv.df)  
  return(plot.df)
}


plot_enrich <- function(plot.df){
  plt <- ggplot(data=plot.df,aes(x=tissue,y=log_pval)) + 
    geom_point(shape=21,size=2,aes(fill=geneset)) + 
    geom_hline(yintercept = -log(0.05/54,base=10)) + 
    geom_hline(yintercept = -log(0.05,base=10),linetype=2) + 
    scale_fill_manual("Gene Set",values=c("gold1","olivedrab2","brown","red"))+
    ylab("-log10(p)") + xlab("Tissue") + 
    coord_flip()
  return(plt)
}


#Get co-expression enrichments enrichments 


work.dir2 <- proj.dir %&% "revamp/"
cred.df <- fread(work.dir2 %&% "genetic_credible_sets/gencred.txt")
cred.df$Locus.ID <- cred.df$CondID


weighted.df <- fread(work.dir2 %&% "analysis_files/" %&% 
                       "classified-loci_weighted.txt") 


weighted.df$assigned_00 <- gsub("other","unclassified",weighted.df$assigned_00)
weighted.df$assigned_20 <- gsub("other","unclassified",weighted.df$assigned_20)
weighted.df$assigned_50 <- gsub("other","unclassified",weighted.df$assigned_50)
weighted.df$assigned_80 <- gsub("other","unclassified",weighted.df$assigned_80)






build_complete_df <- function(group.df,iter){
  thresh.vec <- c(0,0.2,0.5,0.8)
  out.df <- c()
  for (t in thresh.vec){
    print(t)
    var <- ifelse(t==0,"assigned_00",
                  ifelse(t==0.2,"assigned_20",
                         ifelse(t==0.5,"assigned_50",
                                ifelse(t=0.80,"assigned_80"))))
    sub <- dplyr::select(group.df,one_of(var,"symbol"))
    names(sub)[1] <- "assigned"
    islet.genes <- filter(sub,assigned=="islet")$symbol %>% unique(.)
    muscle.genes <- filter(sub,assigned=="muscle")$symbol %>% unique(.)
    liver.genes <- filter(sub,assigned=="liver")$symbol %>% unique(.)
    adipose.genes <- filter(sub,assigned=="adipose")$symbol %>% unique(.) 
    plot.df <- build_plot_df(islet.genes,muscle.genes,
                             liver.genes,adipose.genes,iter)
    plot.df$threshold <- t
    out.df <- rbind(out.df,plot.df)
  }
  return(out.df)
}


annot.prof.df <- fread(proj.dir %&% "revamp/analysis_files/annotation-divvy-weighted-unscaled.txt")

annot.prof.df$coding <- annot.prof.df$coding_islet + annot.prof.df$coding_adipose + 
  annot.prof.df$coding_muscle + annot.prof.df$coding_liver



zero.vec <- filter(annot.prof.df,coding==0)$Locus.ID %>% unique(.)
p10.vec <- filter(annot.prof.df,coding<=0.1)$Locus.ID %>% unique(.)



df4 <- build_complete_df(filter(weighted.df,Locus.ID %in% p10.vec),iter=10000)
write.table(x=df4,file=work.dir2%&%"analysis_files/coexpress-enrich_weighted_p10Coding.txt",
            sep="\t",row.names=F,quote=F)










