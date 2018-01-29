"%&%" <- function(a,b) paste0(a,b)

library("tidyverse")
library("data.table")
library("Homo.sapiens")


serv.dir1 <- "/Users/jtorres/FUSE/"
serv.dir2 <- "/Users/jtorres/FUSE5/"

proj.dir <- serv.dir2 %&% "projects/t2d_classification/"

ppa.df <- fread(serv.dir2 %&% "projects/t2d_classification/method_A/analysis_files/" %&% "full_best-joint_ppa-prop_divvy-coding.txt")

collapse_ppa_df <- function(ppa.df){
  part.ppa.df <- c()
  pb <- txtProgressBar(min=0,max=dim(ppa.df)[1],style=3)
  for (i in 1:dim(ppa.df)[1]){
    setTxtProgressBar(pb,i)
    r <- ppa.df[i,] %>% as.data.frame(.)
    na.vec <- (1:dim(r)[2])[is.na(r) %>% as.logical(.)]
    r[,na.vec] <- 0
    islet <- dplyr::select(r,one_of("islet.coding","islet_specific_strong_enhancer",
                                    "islet_specific_gene_transcription","islet_specific_bivalent_tss",
                                    "islet_shared_strong_enhancer","islet_shared_promoter")) %>% as.numeric(.) %>% sum(.)
    adipose <- dplyr::select(r,one_of("adipose.coding","adipose_specific_strong_enhancer")) %>% as.numeric(.) %>% sum(.)
    liver <- dplyr::select(r,one_of("liver.coding","liver_specific_strong_enhancer",
                                    "liver_specific_genic_enhancer")) %>% as.numeric(.) %>% sum(.)
    muscle <- dplyr::select(r,one_of("muscle.coding","muscle_specific_weak_enhancer",
                                     "muscle_specific_genic_enhancer")) %>% as.numeric(.) %>% sum(.)
    build.df <- data.frame(Locus.ID=r$Locus.ID,islet,liver,adipose,muscle,stringsAsFactors=FALSE)
    part.ppa.df <- rbind(part.ppa.df,build.df)
  } 
  return(part.ppa.df)
}

df <- collapse_ppa_df(ppa.df)
write.table(df,serv.dir2 %&% "projects/t2d_classification/method_A/analysis_files/" %&% "test_input.txt", quote=FALSE,sep="\t",row.names=FALSE)
