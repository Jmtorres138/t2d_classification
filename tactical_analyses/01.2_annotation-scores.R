
"%&%" <- function(a,b) paste0(a,b)

library("data.table");library("dplyr")

serv.dir <- "/well/mccarthy/users/jason/" # "/home/jason/science/servers/FUSE5/" #
proj.dir <- serv.dir %&% "projects/t2d_classification/"
work.dir <- proj.dir %&% "tactical_analyses/"
file.dir <- work.dir %&% "analysis_files/"
input.dir <- work.dir %&% "input_files/"
tissue_annotation_file <- input.dir%&%"tissue-annotations.txt"
genomic_annotation_file <- input.dir%&%"genomic-annotations.txt"
ess.annot <-  "coding"
ess.file <- input.dir%&%"gene-expression-specificity-scores.txt"

snp.annotated.df <- fread(file.dir%&%"annotated_snp-file.txt")


# Partition function 

calculate_annotation_scores <- function(snp.annotated.df, tissue_annotation_file, genomic_annotation_file, 
          ess.annot = NULL, ess.file = NULL) 
{
  "%&%" <- function(a, b) paste0(a, b)
  "%>%" <- magrittr::"%>%"
  tiss.annot.df <- data.table::fread(tissue_annotation_file)
  gen.annot.df <- data.table::fread(genomic_annotation_file)
  if (!is.null(ess.annot)) {
    ess.df <- data.table::fread(ess.file)
  }
  out.df <- c()
  pb <- txtProgressBar(min = 1, max = dim(snp.annotated.df)[1], 
                       style = 3)
  for (r in 1:dim(snp.annotated.df)[1]) {
    full.row.df <- snp.annotated.df[r, ]
    row.df <- full.row.df %>% dplyr::select(., -one_of("SIGNAL", 
                                                       "SNPID", "CHR", "POS", "VALUE"))
    name.vec <- names(row.df)
    tiss.vec <- purrr::map(name.vec, function(s) {
      vec <- strsplit(x = s, split = ".", fixed = T)[[1]]
      ifelse(length(vec) == 2, vec[1], NA)
    }) %>% unlist(.) %>% na.omit(.) %>% as.character(.) %>% 
      unique(.) %>% sort(.)
    tiss.annot.vec <- purrr::map(name.vec, function(s) {
      vec <- strsplit(x = s, split = ".", fixed = T)[[1]]
      ifelse(length(vec) == 2, vec[2], NA)
    }) %>% unlist(.) %>% na.omit(.) %>% as.character(.) %>% 
      unique(.) %>% sort(.)
    gen.annot.vec <- purrr::map(name.vec, function(s) {
      vec <- strsplit(x = s, split = ".", fixed = T)[[1]]
      ifelse(length(vec) == 1, vec[1], NA)
    }) %>% unlist(.) %>% na.omit(.) %>% as.character(.) %>% 
      unique(.) %>% sort(.)
    annot.vec <- c(tiss.annot.vec, gen.annot.vec)
    annot.matrix <- matrix(nrow = length(annot.vec), ncol = length(tiss.vec))
    row.names(annot.matrix) <- annot.vec
    colnames(annot.matrix) <- tiss.vec
    for (i in 1:length(annot.vec)) {
      annot <- annot.vec[i]
      print.warning <- TRUE
      for (e in 1:length(tiss.vec)) {
        tiss <- tiss.vec[e]
        aname <- ifelse(annot %in% tiss.annot.vec, tiss %&% 
                          "." %&% annot, annot)
        wgt <- ifelse(annot %in% tiss.annot.vec, dplyr::filter(tiss.annot.df, 
                                                               V1 == tiss, V2 == annot)$V3, dplyr::filter(gen.annot.df, 
                                                                                                          V1 == annot)$V2)
        if (!is.null(ess.annot)) {
          if (aname == ess.annot) {
            eval.snp <- dplyr::select(row.df, one_of(aname)) %>% 
              as.numeric(.)
            if (eval.snp == 1) {
              snp.chr <- full.row.df$CHR
              snp.pos <- full.row.df$POS
              sub.df <- dplyr::filter(ess.df, CHR == 
                                        snp.chr, START <= snp.pos, END >= snp.pos)
              if (dim(sub.df)[1] == 1) {
                annot.matrix[i, e] <- (dplyr::select(sub.df, 
                                                     one_of(tiss)) %>% as.numeric(.)) * 
                  wgt
              }
              else if (dim(sub.df)[1] == 0) {
                if (print.warning == TRUE) {
                  write("\nWARNING: SNP " %&% full.row.df$SNPID %&% 
                          " at SIGNAL " %&% full.row.df$SIGNAL %&% 
                          " maps to annotation " %&% aname %&% 
                          " but SNP does not map to feature in provided specificity file, please inspect", 
                        stdout())
                  print.warning <- FALSE
                }
                default.val <- 1/length(tiss.vec)
                annot.matrix[i, e] <- default.val * wgt
              }
              else if (dim(sub.df)[1] > 1) {
                if (print.warning == TRUE) {
                  write("\nWARNING: SNP " %&% full.row.df$SNPID %&% 
                          " at SIGNAL " %&% full.row.df$SIGNAL %&% 
                          " maps to annotation " %&% aname %&% 
                          " but SNP maps to multiple features in provided specificity file, please inspect", 
                        stdout())
                  print.warning <- FALSE
                }
                default.val <- 1/length(tiss.vec)
                annot.matrix[i, e] <- default.val * wgt
              }
              else {
                write("\nError linking annotated SNP to specificity file, please inspect", 
                      stderr())
              }
            }
            else {
              annot.matrix[i, e] <- (dplyr::select(row.df, 
                                                   one_of(aname)) %>% as.numeric(.)) * wgt
            }
          }
          else {
            annot.matrix[i, e] <- (dplyr::select(row.df, 
                                                 one_of(aname)) %>% as.numeric(.)) * wgt
          }
        }
        else {
          annot.matrix[i, e] <- (dplyr::select(row.df, 
                                               one_of(aname)) %>% as.numeric(.)) * wgt
        }
      }
    }
    if (sum(colSums(annot.matrix)) > 0) {
      wgt.matrix <- annot.matrix/sum(colSums(annot.matrix))
    } else {
      wgt.matrix <- annot.matrix
    }
    score.matrix <- full.row.df$VALUE * wgt.matrix
    a.vec <- row.names(score.matrix)
    t.vec <- colnames(score.matrix)
    col.names.vec <- c()
    col.vec <- c() 
    for (a in 1:length(a.vec)){
      for (t in 1:length(t.vec)){
        nme <- a.vec[a] %&% "." %&% t.vec[t]
        val <- score.matrix[a,t]
        col.names.vec <- append(col.names.vec,nme)
        col.vec <- append(col.vec,val)
      }
    }
    names(col.vec) <- col.names.vec
    score.df <- col.vec %>% t(.) %>% as.data.frame(.)

    sub.df <- dplyr::select(full.row.df, one_of("SIGNAL", 
                                                "SNPID", "CHR", "POS", "VALUE"))
    build.df <- cbind(sub.df, score.df)
    out.df <- rbind(out.df, build.df)
    setTxtProgressBar(pb, r)
  }
  return(out.df)
}

collapse_annotation_scores <- function(ascore.df){
  sig.vec <- ascore.df$SIGNAL %>% unique(.)
  out.df <- c() 
  pb <- txtProgressBar(min=0,max=length(sig.vec),style=3)
  for (i in 1:length(sig.vec)){
    setTxtProgressBar(pb,i)
    sig <- sig.vec[i]
    sub.df <- filter(ascore.df,SIGNAL==sig) %>% 
      dplyr::select(.,-one_of("SIGNAL","SNPID","CHR","POS"))
    header.vec <- names(sub.df)
    sub.mat <- as.matrix(sub.df)
    collapse.df <- colSums(sub.mat) %>% t(.) %>% as.data.frame(.) 
    build.df <- cbind(data.frame("SIGNAL"=sig,stringsAsFactors = F),collapse.df)
    out.df <- rbind(out.df,build.df)
  }
  return(out.df)
}

ascore.df <- calculate_annotation_scores(snp.annotated.df, tissue_annotation_file, genomic_annotation_file, 
                                       ess.annot, ess.file)
write.table(x=ascore.df,file=file.dir%&%"annotation-scores-allSNPs.txt",quote=F,sep="\t",row.names=F)

collapse.df <- collapse_annotation_scores(ascore.df)
write.table(x=collapse.df,file=file.dir%&%"annotation-scores-signals.txt",quote=F,sep="\t",row.names=F)



