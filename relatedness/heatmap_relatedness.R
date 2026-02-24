rm(list=ls())

rel_data <- read.table("concat_chr_GQ_relatedness_set4.relatedness_matrix",header=TRUE, row.names=1)
heatmap(as.matrix(rel_data),Rowv=NA,Colv=NA,scale="none",col=hcl.colors(20, palette = "reds 2", alpha = 1,rev=TRUE),main="Set4_all")



