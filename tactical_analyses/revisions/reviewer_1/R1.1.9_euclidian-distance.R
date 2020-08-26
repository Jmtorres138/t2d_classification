
"%&%" <- function(a,b) paste0(a,b)
library("tidyverse")
library("data.table")


serv.dir <- "/well/mccarthy/users/jason/" #"/Users/jasont/science/servers/FUSE5/"
proj.dir <- serv.dir %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "tactical_analyses/revisions/reviewer_1/"

tpm.dir <- proj.dir %&% "tactical_analyses/revisions/reviewer_1/analysis_files/"
isl.tpm.df <- fread(tpm.dir %&% "islet_expression_tpm_v2.txt")
liv.tpm.df <- fread(tpm.dir %&% "liver_expression_tpm.txt")
mus.tpm.df <- fread(tpm.dir %&% "muscle_expression_tpm.txt")
adisub.tpm.df <- fread(tpm.dir %&% "adipose_expression_tpm.txt")
adivis.tpm.df <- fread(tpm.dir %&% "adipose-visceral_expression_tpm.txt")
tpm.list <- list(isl.tpm.df,liv.tpm.df,mus.tpm.df,adisub.tpm.df,adivis.tpm.df)
tpm.names <- c("islet","liver","muscle","adipose-sub","adipose-visc")


process_tpm_df <- function(tpm.df,samp.name){
  gene.vec <- tpm.df$GeneName
  df <- tpm.df %>% select(-one_of("GeneID","GeneName")) %>%
    as.matrix(.) %>% t(.) %>% as.data.frame(.)
  names(df) <- gene.vec
  sample <- samp.name %&% "." %&%1:dim(df)[1]
  out.df <- data.frame("Sample"=sample,df)
  names(out.df)[2:dim(out.df)[2]] <- gene.vec
  return(out.df)
}


isl.tpm.df <- process_tpm_df(isl.tpm.df,"isl")
liv.tpm.df <- process_tpm_df(liv.tpm.df,"liv")
mus.tpm.df <- process_tpm_df(mus.tpm.df,"mus")
sat.tpm.df <- process_tpm_df(adisub.tpm.df,"sat")
vat.tpm.df <- process_tpm_df(adivis.tpm.df,"vat")



gene.vec <- names(liv.tpm.df)[2:dim(liv.tpm.df)[2]]
col.vec <- c("Sample",gene.vec)
dup.vec <- col.vec[duplicated(col.vec)]
col.vec <- col.vec[!(col.vec %in% dup.vec)]

isl.tpm.df <- isl.tpm.df[,names(isl.tpm.df) %in% col.vec]
liv.tpm.df <- liv.tpm.df[,names(liv.tpm.df) %in% col.vec]
mus.tpm.df <- mus.tpm.df[,names(mus.tpm.df) %in% col.vec]
sat.tpm.df <- sat.tpm.df[,names(sat.tpm.df) %in% col.vec]
vat.tpm.df <- vat.tpm.df[,names(vat.tpm.df) %in% col.vec]

isl.tpm.df <- dplyr::select(isl.tpm.df,one_of(col.vec))
liv.tpm.df <- dplyr::select(liv.tpm.df,one_of(col.vec))
mus.tpm.df <- dplyr::select(mus.tpm.df,one_of(col.vec))
sat.tpm.df <- dplyr::select(sat.tpm.df,one_of(col.vec))
vat.tpm.df <- dplyr::select(vat.tpm.df,one_of(col.vec))

full.df <- rbind(isl.tpm.df,liv.tpm.df)
full.df <- rbind(full.df,mus.tpm.df)
full.df <- rbind(full.df,sat.tpm.df)
full.df <- rbind(full.df,vat.tpm.df)


full.df <- full.df[, colSums(full.df != 0) > 0]


hc.vec <- 1:dim(full.df)[1]
hc.df <- full.df[hc.vec,] %>% select(-one_of("Sample"))
hc.names <- full.df$Sample[hc.vec]
row.names(hc.df) <- hc.names


hc.df <- scale(hc.df)
dist.mat <- dist(hc.df,method="euclidean")

euc.mat <- dist.mat %>% as.matrix(.)



euc.df <- c()
pb <- txtProgressBar(min=0,max=dim(euc.mat)[1],style = 3) #
track.vec <- c()
for (i in 1:dim(euc.mat)[1]){
  setTxtProgressBar(pb,i)
  samp1 <- row.names(euc.mat)[i]
  group1 <- (samp1 %>% strsplit(.,split=".",fixed=TRUE))[[1]][1]
  for (j in 1:dim(euc.mat)[2]){
    samp2 <- colnames(euc.mat)[j]
    group2 <- (samp2 %>% strsplit(.,split=".",fixed=TRUE))[[1]][1]
    if (samp2 != samp1 & !(samp2 %in% track.vec)){
      euc.dist <- euc.mat[i,j]
      build.df <- data.frame("sample1"=samp1,"sample2"=samp2,
                             "group1"=group1,"group2"=group2,"distance"=euc.dist,
                             stringsAsFactors = FALSE)
      euc.df <- rbind(euc.df,build.df)
    }
  }
  track.vec <- append(track.vec,samp1)
}

write.table(x=euc.df,file=tpm.dir %&% "euclidean-distances.txt",sep="\t",
            row.names=F,quote=F)


diff.df <- c()
group.vec1 <- euc.df$group1 %>% unique(.)
group.vec2 <- euc.df$group2 %>% unique(.)
for (g1 in group.vec1){
  for (g2 in group.vec2){
    sub <- filter(euc.df,group1==g1,group2==g2)
    euc.avg <- sub$distance %>% mean(.)
    build.df <- data.frame("group1"=g1,"group2"=g2,"avg.euc.dist"=euc.avg,
                           stringsAsFactors = FALSE)
    diff.df <- rbind(diff.df,build.df)
  }
}
diff.df <- na.omit(diff.df)

write.table(x=diff.df,file=tpm.dir %&% "average-euclidean-distances.txt",sep="\t",
            row.names=F,quote=F)
