rm(list=ls())

recomb_replicates0 <- read.table("recomb_replicates_v2.txt",header=TRUE) #parent sex offspring chrom left right
recomb_replicates <- subset(recomb_replicates0,parent!="KT_AL1")
recomb_replicates <- subset(recomb_replicates,offspring!="PK_02II_6")

indv_list <- unique(sort(recomb_replicates$offspring))
colony_list <- unique(sort(recomb_replicates$parent))

num_indv <- length(indv_list)

chr_list <- unique(sort(recomb_replicates$chrom))
num_chr <- length(chr_list)

par(mfrow=c(2,2))

col_list <- rep(c("blue4","blue4","sienna3","sienna3"),length.out=(num_indv+1)) #[-16]

diff_count <- 0
diff_count_all <- 0

for (i in 1:num_chr){
  chr_i <- chr_list[i]
  chr_subset <- subset(recomb_replicates,chrom==chr_i)
  x_max <- 1.1*max(chr_subset$right,na.rm = TRUE)
  y_max <- num_indv+1
  plot(0,0,type="p",pch=20,col="NA",xlim=c(-0.08*x_max,x_max),ylim=c(0,y_max),main=chr_i,xlab="",ylab="",yaxt="n")
  #axis(side=2,at=1:num_indv,labels=indv_list,cex.axis=0.7,col=col_list,las=2)
  text(-0.06*x_max,1:num_indv,labels=indv_list,cex=0.7,col=col_list,bg="white")
  for (j in seq(1,num_indv,by=2)){
    indv_j <- indv_list[j]
    indv_chr_subset <- subset(chr_subset,offspring==indv_j)
    points(indv_chr_subset$left,rep(j,nrow(indv_chr_subset)),type="p",pch=8,col=col_list[j])
    points(indv_chr_subset$right,rep(j,nrow(indv_chr_subset)),type="p",pch=18,col=col_list[j])
    jj <- j+1
    indv_jj <- indv_list[jj]
    indv_chr_subset_jj <- subset(chr_subset,offspring==indv_jj)
    points(indv_chr_subset_jj$left,rep(jj,nrow(indv_chr_subset_jj)),type="p",pch=8,col=col_list[jj])
    points(indv_chr_subset_jj$right,rep(jj,nrow(indv_chr_subset_jj)),type="p",pch=18,col=col_list[jj])
    if (nrow(indv_chr_subset)>0 && nrow(indv_chr_subset_jj)>0){
      for (ij in 1:min(c(nrow(indv_chr_subset),nrow(indv_chr_subset_jj)))){
        if ((indv_chr_subset[ij,"left"]!=indv_chr_subset_jj[ij,"left"]) || (indv_chr_subset[ij,"right"]!=indv_chr_subset_jj[ij,"right"])){
          diff_count_all <- diff_count_all+1
          if ((min(c(indv_chr_subset[ij,"right"],indv_chr_subset_jj[ij,"right"]))<max(c(indv_chr_subset[ij,"left"],indv_chr_subset_jj[ij,"left"])))){ # not overlapping
            print(indv_chr_subset[ij,])
            print(indv_chr_subset_jj[ij,])
            diff_count <- diff_count+1
            print("*")
          }
        }
      }
    }
    if (nrow(indv_chr_subset)>nrow(indv_chr_subset_jj)){
      print(indv_chr_subset[nrow(indv_chr_subset_jj):nrow(indv_chr_subset),])
      print("*****")
      diff_count <- diff_count + (nrow(indv_chr_subset) - nrow(indv_chr_subset_jj))
      diff_count_all <- diff_count_all + (nrow(indv_chr_subset) - nrow(indv_chr_subset_jj))
    }
    else if (nrow(indv_chr_subset)<nrow(indv_chr_subset_jj)){
      print(indv_chr_subset_jj[nrow(indv_chr_subset):nrow(indv_chr_subset_jj),])
      print("*****")
      diff_count <- diff_count + (nrow(indv_chr_subset_jj) - nrow(indv_chr_subset))
      diff_count_all <- diff_count_all + (nrow(indv_chr_subset_jj) - nrow(indv_chr_subset))
    }
  }
}
