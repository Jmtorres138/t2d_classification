"%&%" <- function(a,b) paste0(a,b)
library("tidyverse")
library("GenomicRanges")
library("viridis")
library("data.table")
library("RColorBrewer")

#serv.dir <- "/Users/jtorres/FUSE5/"
serv.dir <- "/home/jason/science/servers/FUSE5/"
proj.dir <- serv.dir %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "method_C/"

keep.vec <- (fread(work.dir %&% "analysis_files/maxppa.txt") %>% filter(.,max.ppa >=0.5))$Locus.ID
df <- fread(work.dir %&% "analysis_files/classified-loci_weighted_renamed.txt") %>% 
  filter(.,Locus.ID %in% keep.vec)

# Islet 
isl <-filter(df,islet>=0.9)
liv <-filter(df,liver>=0.9)
adi <-filter(df,adipose>=0.9)
mus <-filter(df,muscle>=0.9)

mat <- select(df,one_of("islet","muscle","adipose","liver")) %>% as.matrix(.)

check.vec <- c()
for (i in 1:dim(mat)[1]){
  vec <- mat[i,]
  eval <- c()
  for (e in vec){
    if (e >=0.25){
      eval <- append(eval,1)
    } else{
      eval <- append(eval,0)
    }
  }
  check.vec <- append(check.vec,sum(eval))
}
(check.vec >=2) %>% sum(.)


df2 <- fread(work.dir %&% "analysis_files/classified-loci_weighted_renamed.txt") 
arrange(df2,desc(muscle)) %>% head(.)


## PPI Project Check 


mult.df <- df[(check.vec >=2),] # TOA screen 

get_med_exp <- function(tiss.df){
  pb <- txtProgressBar(min=0,max=dim(tiss.df)[1],style=3)
  m <- select(tiss.df, -one_of("GeneID","GeneName")) %>% as.matrix(.)
  med.vec <- map(1:dim(m)[1],function(i){
    setTxtProgressBar(pb,i)
    median(m[i,])
  }) %>% as.numeric(.)
  return(med.vec)
}

#expression 
isl.exp <- fread(proj.dir %&% "analysis_files/islet_expression_tpm.txt")
mus.exp <- fread(proj.dir %&% "analysis_files/muscle_expression_tpm.txt")
adi.exp <- fread(proj.dir %&% "analysis_files/adipose_expression_tpm.txt")
liv.exp <- fread(proj.dir %&% "analysis_files/liver_expression_tpm.txt")

islet.medexp <- get_med_exp(isl.exp)
adipose.medexp <- get_med_exp(adi.exp)
muscle.medexp <- get_med_exp(mus.exp)
liver.medexp <- get_med_exp(liv.exp)

isl.df <- cbind(select(isl.exp,one_of("GeneID","GeneName")),medexp=get_med_exp(isl.exp))
adi.df <- cbind(select(adi.exp,one_of("GeneID","GeneName")),medexp=get_med_exp(adi.exp))
mus.df <- cbind(select(mus.exp,one_of("GeneID","GeneName")),medexp=get_med_exp(mus.exp))
liv.df <- cbind(select(liv.exp,one_of("GeneID","GeneName")),medexp=get_med_exp(liv.exp))


mult.df$symbol[c(7,23,43)] <- c("COBLL1","CAMK1D","APOE")#c("COBLL1","CAMK1D","TOMM40")
mult.df2 <- c()
for (i in 1:dim(mult.df)[1]){
  print(i)
  sub <- mult.df[i,]
  g <- sub$symbol
  islet.expr <- filter(isl.df,GeneName==g)$medexp;   adipose.expr <- filter(adi.df,GeneName==g)$medexp
  liver.expr <- filter(liv.df,GeneName==g)$medexp;  muscle.expr <- filter(mus.df,GeneName==g)$medexp
  if (length(islet.expr)==0){
    islet.expr <- NA
  }
  if (length(adipose.expr)==0){
    adipose.expr <- NA
  }
  if (length(liver.expr)==0){
    liver.expr <- NA
  }
  if (length(muscle.expr)==0){
    muscle.expr <- NA
  }
  e.df <- cbind(islet.expr,muscle.expr,adipose.expr,liver.expr)
  build.df <- cbind(sub,e.df)
  mult.df2 <- rbind(mult.df2,build.df)
}

mult.df2 <- na.omit(mult.df2)

# filter expression median TPM > 1 in more than one tissue 

keep.index <- c() 
for (i in 1:dim(mult.df2)[1]){
  sub <- select(mult.df2[i,],contains("expr")) %>% as.numeric(.)
  check.vec <- c()
  for (e in sub){
    if (e > 1){
      check.vec <- append(check.vec,1)
    } else{
      check.vec <- append(check.vec,0) 
    }
  }
  if (sum(check.vec)>=2){
    keep.index <- append(keep.index,i)
  }
}

mult.df3 <- mult.df2[keep.index,] %>% select(.,-contains("assigned"))


# Coding check 

targ.df <- fread(work.dir%&%"analysis_files/results-maxppa50_weighted.txt") %>% 
  select(.,one_of("Locus.ID","coding","assigned_20","eGene","CLPP","mixed.locus",
                  "secondary.tissue","secondary.eGene","secondary.CLPP"))
df4 <- inner_join(mult.df3,targ.df,by="Locus.ID")
keep.index2 <- c()
for (i in 1:dim(df4)[1]){
  sub <- df4[i,]
  if (sub$coding > 0.10 | sub$eGene != ""){
    keep.index2 <- append(keep.index2,i)
  }
}

df5 <- df4[keep.index2,]
write.table(x=df5,file=work.dir %&% "analysis_files/nih-PPI.txt",sep="\t",quote=F,row.names=F)
