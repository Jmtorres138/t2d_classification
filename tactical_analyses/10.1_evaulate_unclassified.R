
"%&%" <- function(a,b) paste0(a,b)
library("tidyverse")
library("data.table")
library("GenomicRanges")

fuse.path <- "/home/jason/science/servers/"
serv.dir <- fuse.path %&% "FUSE5/"
serv.dir0 <- fuse.path %&% "FUSE/"

proj.dir <- serv.dir %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "tactical_analyses/"
analysis.dir <- work.dir %&% "analysis_files/"

cred.df <- fread(work.dir %&% "genetic_credible_sets/gencred.txt")
class.df <- fread(analysis.dir %&% "classified-loci_weighted_with-shared.txt")

states.df <- fread(serv.dir0 %&% "reference/chromatin_segmentation/varshney_2016/bed_files/complete_varshney_chromHMM_states.bed")



df_to_gr <- function(df){
  gr <- GRanges(seqnames=df$V1,IRanges(start=df$V2,end=df$V3))
  return(gr)
}


t2d.sub <- filter(states.df,grepl(pattern = "Liver",x=V4) | 
                    grepl(pattern = "Islet",x=V4) | 
                    grepl(pattern = "Adipose",x=V4) | 
                    grepl(pattern = "SkeletalMuscle",x=V4))


t2d.qui <- filter(t2d.sub,grepl(pattern = "_18_",x=V4))
t2d.qui.isl <- filter(t2d.qui,grepl(pattern = "Islet",x=V4)) %>% df_to_gr(.)
t2d.qui.liv <- filter(t2d.qui,grepl(pattern = "Liver",x=V4)) %>% df_to_gr(.)
t2d.qui.adi <- filter(t2d.qui,grepl(pattern = "Adipose",x=V4)) %>% df_to_gr(.)
t2d.qui.mus <- filter(t2d.qui,grepl(pattern = "SkeletalMuscle",x=V4)) %>%
  df_to_gr(.)
t2d.qui <- intersect(t2d.qui.isl,t2d.qui.liv) %>% intersect(.,t2d.qui.adi) %>% 
  intersect(.,t2d.qui.mus) %>% as.data.frame()
names(t2d.qui)[1:3] <- c("V1","V2","V3")


t2d.rep <- filter(t2d.sub,grepl(pattern = "epressed_",x=V4))
t2d.rep.isl <- filter(t2d.rep,grepl(pattern = "Islet",x=V4)) %>% df_to_gr(.)
t2d.rep.liv <- filter(t2d.rep,grepl(pattern = "Liver",x=V4)) %>% df_to_gr(.)
t2d.rep.adi <- filter(t2d.rep,grepl(pattern = "Adipose",x=V4)) %>% df_to_gr(.)
t2d.rep.mus <- filter(t2d.rep,grepl(pattern = "SkeletalMuscle",x=V4)) %>%
  df_to_gr(.)
t2d.rep <- intersect(t2d.rep.isl,t2d.rep.liv) %>% intersect(.,t2d.rep.adi) %>% 
  intersect(.,t2d.rep.mus) %>% as.data.frame()
names(t2d.rep)[1:3] <- c("V1","V2","V3")


brain.enh <- filter(states.df,grepl(pattern = "Active_enhancer",x=V4)) %>% 
  filter(.,grepl(pattern = "SubstantiaNigra",x=V4) | 
           grepl(pattern = "CingulateGyrus",x=V4) | 
           grepl(pattern = "HippocampusMiddle",x=V4) | 
           grepl(pattern = "MidFrontalLobe",x=V4) | 
           grepl(pattern = "AnteriorCaudate",x=V4) | 
           grepl(pattern = "InferiorTemporalLobe",x=V4))

brain.tss <- filter(states.df,grepl(pattern = "Active_TSS",x=V4)) %>% 
  filter(.,grepl(pattern = "SubstantiaNigra",x=V4) | 
           grepl(pattern = "CingulateGyrus",x=V4) | 
           grepl(pattern = "HippocampusMiddle",x=V4) | 
           grepl(pattern = "MidFrontalLobe",x=V4) | 
           grepl(pattern = "AnteriorCaudate",x=V4) | 
           grepl(pattern = "InferiorTemporalLobe",x=V4))                      



get_cum_score <- function(assigned="unclassified",query.df){
  query.gr <- GRanges(seqnames = query.df$V1, IRanges(query.df$V2,query.df$V3))
  sigs <- filter(class.df,assigned_20==assigned)$Locus.ID
  cumsum.total <- filter(cred.df,CondID%in%sigs)$PPA %>% sum(.)
  cumsum.state <- 0 
  out.df <- c()
  pb <- txtProgressBar(min=1,max=length(sigs),style=3)
  for (i in 1:length(sigs)){
    setTxtProgressBar(pb,i)
    sig <- sigs[i]
    sig.sum <- 0 
    cred.sub <- filter(cred.df,CondID==sig)
    for (e in 1:dim(cred.sub)[1]){
      row.df <- cred.sub[e,]
      row.gr <- GRanges(seqnames = row.df$CHR,IRanges(row.df$POS,row.df$POS))
      if ((row.gr %over% query.gr)==TRUE){
        sig.sum <- sig.sum + row.df$PPA; 
        cumsum.state <- cumsum.state + row.df$PPA
      }
    }
    build.df <- data.frame("Locus.ID"=sig,"CumPPA"=sig.sum,stringsAsFactors = F)
    out.df <- rbind(out.df,build.df)
  }
  return(list(cumsum.total,cumsum.state,out.df))
}



unc.t2d.qui <- get_cum_score("unclassified",t2d.qui)
unc.t2d.rep <- get_cum_score("unclassified",t2d.rep)
unc.brain.enh <- get_cum_score("unclassified",brain.enh)
unc.brain.tss <- get_cum_score("unclassified",brain.tss)

isl.t2d.qui <- get_cum_score("islet",t2d.qui)
isl.t2d.rep <- get_cum_score("islet",t2d.rep)
isl.brain.enh <- get_cum_score("islet",brain.enh)
isl.brain.tss <- get_cum_score("islet",brain.tss)

adi.t2d.qui <- get_cum_score("adipose",t2d.qui)
adi.t2d.rep <- get_cum_score("adipose",t2d.rep)
adi.brain.enh <- get_cum_score("adipose",brain.enh)
adi.brain.tss <- get_cum_score("adipose",brain.tss)

mus.t2d.qui <- get_cum_score("muscle",t2d.qui)
mus.t2d.rep <- get_cum_score("muscle",t2d.rep)
mus.brain.enh <- get_cum_score("muscle",brain.enh)
mus.brain.tss <- get_cum_score("muscle",brain.tss)

liv.t2d.qui <- get_cum_score("liver",t2d.qui)
liv.t2d.rep <- get_cum_score("liver",t2d.rep)
liv.brain.enh <- get_cum_score("liver",brain.enh)
liv.brain.tss <- get_cum_score("liver",brain.tss)

sha.t2d.qui <- get_cum_score("shared",t2d.qui)
sha.t2d.rep <- get_cum_score("shared",t2d.rep)
sha.brain.enh <- get_cum_score("shared",brain.enh)
sha.brain.tss <- get_cum_score("shared",brain.tss)



# Repressed regions in the four T2D tissues



unc.t2d.rep[[3]]$assigned_20 <- "unclassified"
isl.t2d.rep[[3]]$assigned_20 <- "islet"
mus.t2d.rep[[3]]$assigned_20 <- "muscle"
liv.t2d.rep[[3]]$assigned_20 <- "liver"
adi.t2d.rep[[3]]$assigned_20 <- "adipose"
sha.t2d.rep[[3]]$assigned_20 <- "shared"

rep.df <- rbind(unc.t2d.rep[[3]],isl.t2d.rep[[3]],mus.t2d.rep[[3]],
                liv.t2d.rep[[3]],adi.t2d.rep[[3]],sha.t2d.rep[[3]])
rep.df$classified <- map(rep.df$assigned_20,function(s){
  ifelse(s=="unclassified","unclassified","classified")
}) %>% as.character(.)

rep.plot <- ggplot(data=rep.df,aes(x=classified,y=CumPPA)) + 
  geom_boxplot() + 
  geom_jitter(shape=15,color="steelblue",position=position_jitter(0.21))

filter(rep.df,classified=="classified")$CumPPA %>% mean(.)
filter(rep.df,classified=="unclassified")$CumPPA %>% mean(.)
wilcox.test(filter(rep.df,classified=="classified")$CumPPA,
            filter(rep.df,classified=="unclassified")$CumPPA)



# Quiescent regions in the four T2D tissues 



unc.t2d.qui[[3]]$assigned_20 <- "unclassified"
isl.t2d.qui[[3]]$assigned_20 <- "islet"
mus.t2d.qui[[3]]$assigned_20 <- "muscle"
liv.t2d.qui[[3]]$assigned_20 <- "liver"
adi.t2d.qui[[3]]$assigned_20 <- "adipose"
sha.t2d.qui[[3]]$assigned_20 <- "shared"

qui.df <- rbind(unc.t2d.qui[[3]],isl.t2d.qui[[3]],mus.t2d.qui[[3]],
                liv.t2d.qui[[3]],adi.t2d.qui[[3]],sha.t2d.qui[[3]])
qui.df$classified <- map(qui.df$assigned_20,function(s){
  ifelse(s=="unclassified","unclassified","classified")
}) %>% as.character(.)

qui.plot <- ggplot(data=qui.df,aes(x=classified,y=CumPPA)) + 
  geom_boxplot() + 
  geom_jitter(shape=15,color="steelblue",position=position_jitter(0.21))
filter(qui.df,classified=="classified")$CumPPA %>% mean(.)
filter(qui.df,classified=="unclassified")$CumPPA %>% mean(.)
wilcox.test(filter(qui.df,classified=="classified")$CumPPA,
            filter(qui.df,classified=="unclassified")$CumPPA)


filter(qui.df,classified=="unclassified")$CumPPA %>% mean(.)
filter(rep.df,classified=="unclassified")$CumPPA %>% mean(.)
wilcox.test(filter(qui.df,classified=="unclassified")$CumPPA,
            filter(rep.df,classified=="unclassified")$CumPPA)



# Brain enhancers 


unc.brain.enh[[3]]$assigned_20 <- "unclassified"
isl.brain.enh[[3]]$assigned_20 <- "islet"
mus.brain.enh[[3]]$assigned_20 <- "muscle"
liv.brain.enh[[3]]$assigned_20 <- "liver"
adi.brain.enh[[3]]$assigned_20 <- "adipose"
sha.brain.enh[[3]]$assigned_20 <- "shared"

brain.enh.df <- rbind(unc.brain.enh[[3]],isl.brain.enh[[3]],mus.brain.enh[[3]],
                      liv.brain.enh[[3]],adi.brain.enh[[3]],sha.brain.enh[[3]])
brain.enh.df$classified <- map(brain.enh.df$assigned_20,function(s){
  ifelse(s=="unclassified","unclassified","classified")
}) %>% as.character(.)

brain.enh.plot <- ggplot(data=brain.enh.df,aes(x=assigned_20,y=CumPPA)) + 
  geom_boxplot() + 
  geom_jitter(shape=15,color="steelblue",position=position_jitter(0.21))
brain.enh.plot.v2 <- ggplot(data=brain.enh.df,aes(x=classified,y=CumPPA)) + 
  geom_boxplot() + 
  geom_jitter(shape=15,color="steelblue",position=position_jitter(0.21))
filter(brain.enh.df,classified=="classified")$CumPPA %>% mean(.)
filter(brain.enh.df,classified=="unclassified")$CumPPA %>% mean(.)
wilcox.test(filter(brain.enh.df,classified=="classified")$CumPPA,
            filter(brain.enh.df,classified=="unclassified")$CumPPA)


# Brain promoters 


unc.brain.tss[[3]]$assigned_20 <- "unclassified"
isl.brain.tss[[3]]$assigned_20 <- "islet"
mus.brain.tss[[3]]$assigned_20 <- "muscle"
liv.brain.tss[[3]]$assigned_20 <- "liver"
adi.brain.tss[[3]]$assigned_20 <- "adipose"
sha.brain.tss[[3]]$assigned_20 <- "shared"

brain.tss.df <- rbind(unc.brain.tss[[3]],isl.brain.tss[[3]],mus.brain.tss[[3]],
                      liv.brain.tss[[3]],adi.brain.tss[[3]],sha.brain.tss[[3]])
brain.tss.df$classified <- map(brain.tss.df$assigned_20,function(s){
  ifelse(s=="unclassified","unclassified","classified")
}) %>% as.character(.)

brain.tss.plot <- ggplot(data=brain.tss.df,aes(x=assigned_20,y=CumPPA)) + 
  geom_boxplot() + 
  geom_jitter(shape=15,color="steelblue",position=position_jitter(0.21))
brain.tss.plot.v2 <- ggplot(data=brain.tss.df,aes(x=classified,y=CumPPA)) + 
  geom_boxplot() + 
  geom_jitter(shape=15,color="steelblue",position=position_jitter(0.21))

filter(brain.tss.df,classified=="classified")$CumPPA %>% mean(.)
filter(brain.tss.df,classified=="unclassified")$CumPPA %>% mean(.)
wilcox.test(filter(brain.tss.df,classified=="classified")$CumPPA,
            filter(brain.tss.df,classified=="unclassified")$CumPPA)


# Evaluate Fine-mapping resolution 


class.df$classified <- map(class.df$assigned_20,function(s){
  ifelse(s=="unclassified","unclassified","classified")
}) %>% as.character(.)
class.df$maxppa <- map(class.df$Locus.ID,function(id){
  (filter(cred.df,CondID==id) %>% arrange(.,desc(PPA)))$PPA[1]
}) %>% as.numeric(.)
class.df$ncredsnps <- map(class.df$Locus.ID,function(id){
  (filter(cred.df,CondID==id) %>% dim(.))[1]
}) %>% as.integer(.)


filter(class.df,classified=="classified")$maxppa %>% mean(.)
filter(class.df,classified=="classified")$maxppa %>% median(.)
filter(class.df,classified=="unclassified")$maxppa %>% mean(.)
filter(class.df,classified=="unclassified")$maxppa %>% median(.)
wilcox.test(filter(class.df,classified=="classified")$maxppa,
            filter(class.df,classified=="unclassified")$maxppa)

filter(class.df,classified=="classified")$ncredsnps %>% mean(.)
filter(class.df,classified=="classified")$ncredsnps %>% median(.)
filter(class.df,classified=="unclassified")$ncredsnps %>% mean(.)
filter(class.df,classified=="unclassified")$ncredsnps %>% median(.)
wilcox.test(filter(class.df,classified=="classified")$ncredsnps,
            filter(class.df,classified=="unclassified")$ncredsnps)



rep.class <- filter(rep.df,classified=="classified")$CumPPA %>% mean(.)
rep.unclass <- filter(rep.df,classified=="unclassified")$CumPPA %>% mean(.)
rep.pval <- wilcox.test(filter(rep.df,classified=="classified")$CumPPA,
                        filter(rep.df,classified=="unclassified")$CumPPA)$p.value

qui.class <- filter(qui.df,classified=="classified")$CumPPA %>% mean(.)
qui.unclass <- filter(qui.df,classified=="unclassified")$CumPPA %>% mean(.)
qui.pval <- wilcox.test(filter(qui.df,classified=="classified")$CumPPA,
                        filter(qui.df,classified=="unclassified")$CumPPA)$p.value

brain.enh.class <- filter(brain.enh.df,classified=="classified")$CumPPA %>% mean(.)
brain.enh.unclass <- filter(brain.enh.df,classified=="unclassified")$CumPPA %>% mean(.)
brain.enh.pval <- wilcox.test(filter(brain.enh.df,classified=="classified")$CumPPA,
                              filter(brain.enh.df,classified=="unclassified")$CumPPA)$p.value

brain.tss.class <- filter(brain.tss.df,classified=="classified")$CumPPA %>% mean(.)
brain.tss.unclass <- filter(brain.tss.df,classified=="unclassified")$CumPPA %>% mean(.)
brain.tss.pval <- wilcox.test(filter(brain.tss.df,classified=="classified")$CumPPA,
                              filter(brain.tss.df,classified=="unclassified")$CumPPA)$p.value


out.df <- data.frame(annotation_class=c("Repressed.T2Dtiss","Quiescent.T2Dtiss",
                                        "Enhancer.Brain","ActiveTSS.Brain"),
                     mean.cum.ppa.classified=c(rep.class,qui.class,brain.enh.class,brain.tss.class),
                     mean.cum.ppa.unclassified=c(rep.unclass,qui.unclass,brain.enh.unclass,brain.tss.unclass),
                     pval=c(rep.pval,qui.pval,brain.enh.pval,brain.tss.pval),stringsAsFactors = F)

write.table(x=out.df,file=analysis.dir%&%"brain-enrich.txt",sep="\t",quote=F,row.names=F)



