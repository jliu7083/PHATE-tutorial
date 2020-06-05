#!/usr/bin/env Rscript
########---------  optparse------------------------------------
library("optparse")
option_list<- list(
  make_option(c("--Input"), type="character", default=NULL, help="input count matrix file", metavar="character"),
  make_option(c("--Output"), type="character", default=NULL, help="output the 2D PHATE plot", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

#########--------------PHATE--------------------------######################
library(phateR)
library(ggplot2)
library(readr)
library(viridis)
#install.packages("Rmagic")
#library(Rmagic)
## load count matrix
counts<-read.table(opt$Input, sep = '\t', header = TRUE)
rownames(counts) <-counts$Feature
counts$Feature <- NULL
counts_transpose <- as.data.frame(t(as.matrix(counts)))
counts_transpose[1:10,1:10]
## normalize
count <- library.size.normalize(counts_transpose)
count <- sqrt(count)
#count[1:10,1:10]
count_PHATE <- phate(counts_transpose,t=10)
#plot
ggplot(count_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color=count$LTBR),shape=20) +
  labs(color="LTBR") +
  scale_color_viridis(option="C")
ggsave(opt$Output, width=4, height=4)