
"%&%" <- function(a,b) paste0(a,b)
library("tidyverse")
library("data.table")
library("GenomicRanges")

fuse.path <- "/Users/jasont/science/servers/"
serv.dir <- fuse.path %&% "FUSE5/"
serv.dir0 <- fuse.path %&% "FUSE/"

proj.dir <- serv.dir %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "tactical_analyses/revisions/reviewer_1/"
analysis.dir <- work.dir %&% "analysis_files/"

cred.df <- fread(proj.dir %&%"tactical_analyses/genetic_credible_sets/gencred.txt")
class.df <- fread(proj.dir %&% "tactical_analyses/analysis_files/"
                  %&% "classified-loci_weighted_with-shared.txt")
het.df <- fread(work.dir %&% "bmi-het-signals.txt")
het.sigs <- het.df$Locus.ID 
non.het.sigs <- class.df$Locus.ID[!(class.df$Locus.ID %in% het.sigs)]

states.df <- fread(serv.dir0 %&%
 "reference/chromatin_segmentation/varshney_2016/bed_files/complete_varshney_chromHMM_states.bed")

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



get_cum_score <- function(sigs,query.df){
  query.gr <- GRanges(seqnames = query.df$V1, IRanges(query.df$V2,query.df$V3))
  #sigs <- filter(class.df,assigned_20==assigned)$Locus.ID
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



het.t2d.qui <- get_cum_score(het.sigs,t2d.qui)
saveRDS(het.t2d.qui,file=work.dir %&% "het.t2d.qui.RDS")
non.t2d.qui <- get_cum_score(non.het.sigs,t2d.qui)
saveRDS(non.t2d.qui,file=work.dir %&% "non.t2d.qui.RDS")

het.t2d.rep <- get_cum_score(het.sigs,t2d.rep)
saveRDS(het.t2d.rep,file=work.dir %&% "het.t2d.rep.RDS")
non.t2d.rep <- get_cum_score(non.het.sigs,t2d.rep)
saveRDS(non.t2d.rep,file=work.dir %&% "non.t2d.rep.RDS")

het.brain.enh <- get_cum_score(het.sigs,brain.enh)
saveRDS(het.brain.enh,file=work.dir %&% "het.brain.enh.RDS")
non.brain.enh <- get_cum_score(non.het.sigs,brain.enh)
saveRDS(non.brain.enh,file=work.dir %&% "non.brain.enh.RDS")

het.brain.tss <- get_cum_score(het.sigs,brain.tss)
saveRDS(het.brain.tss,file=work.dir %&% "het.brain.tss.RDS")
non.brain.tss <- get_cum_score(non.het.sigs,brain.tss)
saveRDS(non.brain.tss,file=work.dir %&% "non.brain.tss.RDS")


# Repressed regions in the four T2D tissues

rep.df <- rbind(het.t2d.rep[[3]],non.t2d.rep[[3]])
rep.df$classified <- map(rep.df$Locus.ID,function(s){
  ifelse(s%in%het.sigs,"hetero","non.hetero")
}) %>% as.character(.)

rep.plot <- ggplot(data=rep.df,aes(x=classified,y=CumPPA)) + 
  geom_boxplot() + 
  geom_jitter(shape=15,color="steelblue",position=position_jitter(0.21))

filter(rep.df,classified=="hetero")$CumPPA %>% mean(.)
filter(rep.df,classified=="non.hetero")$CumPPA %>% mean(.)
wilcox.test(filter(rep.df,classified=="hetero")$CumPPA,
            filter(rep.df,classified=="non.hetero")$CumPPA)

# Quiescent regions in the four T2D tissues 



qui.df <- rbind(het.t2d.qui[[3]],non.t2d.qui[[3]])
qui.df$classified <- map(qui.df$Locus.ID,function(s){
  ifelse(s%in%het.sigs,"hetero","non.hetero")
}) %>% as.character(.)

qui.plot <- ggplot(data=qui.df,aes(x=classified,y=CumPPA)) + 
  geom_boxplot() + 
  geom_jitter(shape=15,color="steelblue",position=position_jitter(0.21))

filter(qui.df,classified=="hetero")$CumPPA %>% mean(.)
filter(qui.df,classified=="non.hetero")$CumPPA %>% mean(.)
wilcox.test(filter(qui.df,classified=="hetero")$CumPPA,
            filter(qui.df,classified=="non.hetero")$CumPPA)




# Brain enhancers 

brain.enh.df <- rbind(het.brain.enh[[3]],non.brain.enh[[3]])
brain.enh.df$classified <- map(brain.enh.df$Locus.ID,function(s){
  ifelse(s%in%het.sigs,"hetero","non.hetero")
}) %>% as.character(.)

brain.enh.plot <- ggplot(data=brain.enh.df,aes(x=classified,y=CumPPA)) + 
  geom_boxplot() + 
  geom_jitter(shape=15,color="steelblue",position=position_jitter(0.21))

filter(brain.enh.df,classified=="hetero")$CumPPA %>% mean(.)
filter(brain.enh.df,classified=="non.hetero")$CumPPA %>% mean(.)
wilcox.test(filter(brain.enh.df,classified=="hetero")$CumPPA,
            filter(brain.enh.df,classified=="non.hetero")$CumPPA)



# Brain promoters 

brain.tss.df <- rbind(het.brain.tss[[3]],non.brain.tss[[3]])
brain.tss.df$classified <- map(brain.tss.df$Locus.ID,function(s){
  ifelse(s%in%het.sigs,"hetero","non.hetero")
}) %>% as.character(.)

brain.tss.plot <- ggplot(data=brain.tss.df,aes(x=classified,y=CumPPA)) + 
  geom_boxplot() + 
  geom_jitter(shape=15,color="steelblue",position=position_jitter(0.21))

filter(brain.tss.df,classified=="hetero")$CumPPA %>% mean(.)
filter(brain.tss.df,classified=="non.hetero")$CumPPA %>% mean(.)
wilcox.test(filter(brain.tss.df,classified=="hetero")$CumPPA,
            filter(brain.tss.df,classified=="non.hetero")$CumPPA)


# Evaluate Fine-mapping resolution 

class.df$classified <- map(class.df$Locus.ID,function(s){
  ifelse(s%in%het.sigs,"hetero","non.hetero")
}) %>% as.character(.)
class.df$maxppa <- map(class.df$Locus.ID,function(id){
  (filter(cred.df,CondID==id) %>% arrange(.,desc(PPA)))$PPA[1]
}) %>% as.numeric(.)
class.df$ncredsnps <- map(class.df$Locus.ID,function(id){
  (filter(cred.df,CondID==id) %>% dim(.))[1]
}) %>% as.integer(.)


filter(class.df,classified=="hetero")$maxppa %>% mean(.)
filter(class.df,classified=="hetero")$maxppa %>% median(.)
filter(class.df,classified=="non.hetero")$maxppa %>% mean(.)
filter(class.df,classified=="non.hetero")$maxppa %>% median(.)
wilcox.test(filter(class.df,classified=="hetero")$maxppa,
            filter(class.df,classified=="non.hetero")$maxppa)

filter(class.df,classified=="hetero")$ncredsnps %>% mean(.)
filter(class.df,classified=="hetero")$ncredsnps %>% median(.)
filter(class.df,classified=="non.hetero")$ncredsnps %>% mean(.)
filter(class.df,classified=="non.hetero")$ncredsnps %>% median(.)
wilcox.test(filter(class.df,classified=="hetero")$ncredsnps,
            filter(class.df,classified=="non.hetero")$ncredsnps)


rep.class <- filter(rep.df,classified=="hetero")$CumPPA %>% mean(.)
rep.unclass <- filter(rep.df,classified=="non.hetero")$CumPPA %>% mean(.)
rep.pval <- wilcox.test(filter(rep.df,classified=="hetero")$CumPPA,
                        filter(rep.df,classified=="non.hetero")$CumPPA)$p.value

qui.class <- filter(qui.df,classified=="hetero")$CumPPA %>% mean(.)
qui.unclass <- filter(qui.df,classified=="non.hetero")$CumPPA %>% mean(.)
qui.pval <- wilcox.test(filter(qui.df,classified=="hetero")$CumPPA,
                        filter(qui.df,classified=="non.hetero")$CumPPA)$p.value

brain.enh.class <- filter(brain.enh.df,classified=="hetero")$CumPPA %>% mean(.)
brain.enh.unclass <- filter(brain.enh.df,classified=="non.hetero")$CumPPA %>% mean(.)
brain.enh.pval <- wilcox.test(filter(brain.enh.df,classified=="hetero")$CumPPA,
                              filter(brain.enh.df,classified=="non.hetero")$CumPPA)$p.value

brain.tss.class <- filter(brain.tss.df,classified=="hetero")$CumPPA %>% mean(.)
brain.tss.unclass <- filter(brain.tss.df,classified=="non.hetero")$CumPPA %>% mean(.)
brain.tss.pval <- wilcox.test(filter(brain.tss.df,classified=="hetero")$CumPPA,
                              filter(brain.tss.df,classified=="non.hetero")$CumPPA)$p.value


out.df <- data.frame(annotation_class=c("Repressed.T2Dtiss","Quiescent.T2Dtiss",
                                        "Enhancer.Brain","ActiveTSS.Brain"),
                     mean.cum.ppa.hetero=c(rep.class,qui.class,brain.enh.class,brain.tss.class),
                     mean.cum.ppa.non.hetero=c(rep.unclass,qui.unclass,brain.enh.unclass,brain.tss.unclass),
                     pval=c(rep.pval,qui.pval,brain.enh.pval,brain.tss.pval),stringsAsFactors = F)

write.table(x=out.df,file=work.dir%&%"heteroBMI-signal-enrich.txt",sep="\t",quote=F,row.names=F)



