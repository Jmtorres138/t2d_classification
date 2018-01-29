set.seed(1)

"%&%" <- function(a,b) paste0(a,b)

library("tidyverse")
library("data.table")
library("Homo.sapiens")


serv.dir1 <- "/well/got2d/jason/"
serv.dir2 <- "/well/mccarthy/users/jason/"

proj.dir <- serv.dir2 %&% "projects/t2d_classification/"

