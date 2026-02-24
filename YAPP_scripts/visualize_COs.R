rm(list=ls())
library(plot3D)

# start from a vcf with the drones and the queen from one colony (all encoded as diploids), the queen should be in the last column
# then extract all the genotypes (gatk VariantsToTable)
# then the queen genotypes should be split to two different columns 
# each drone genotype should be replaced by "1" if it is homozygous for the first queen allele and replaced by "2" if it is homozygous for the second queen allele
# all this could be done with the commands below (in the terminal/Bash)

# tabix queen_phased.vcf.gz
# /Users/turev988/gatk-4.6.1.0/gatk VariantsToTable -V input_vcf_queen.vcf -F CHROM -F POS -GF GT -O input_vcf_queen.vcf_GT
# !!!! make sure queen is the last column!!!!
# !!!! make sure missing encoded as ".|." and not just "." !!!!
# if missing="." only in last column:
# cat owh_from_hieu_SLU3_NC50.vcf_GT | awk '{if ($NF=="\.") {print $0"|."} else {print $0}}' > owh_from_hieu_SLU3_NC50.vcf_GT2

# cat input_vcf_queen.vcf_GT | head -n 1 | awk '{print $0"\t"$NF"_2"}' | sed s%".GT"%""%g > input_vcf_queen.vcf_GT_mod 
# cat input_vcf_queen.vcf_GT | awk '{print $NF}' | sed s/"|"/"\t"/g > input_vcf_queen.vcf_GT_queen
# paste input_vcf_queen.vcf_GT input_vcf_queen.vcf_GT_queen | tail -n +2 | sed s%"|"%"\/"%g | awk '{if ($NF!="\." && $(NF-1)!="\.") {printf $1"\t"$2; for (i=3;i<(NF-2);i++) {if ($i==$(NF-1)"/"$(NF-1)) {printf "\t1"} else if ($i==$NF"/"$NF) {printf "\t2"} else {printf "\tNA"} }; printf "\t"$(NF-1)"\t"$NF"\n" }}' >> input_vcf_queen.vcf_GT_mod

# !!! incorrect handling of sites where queen is homozygous - exclude those!!!!
# cat owh_from_hieu_SLU3_NC50.vcf_GT2_mod | awk '{ if ($NF!=$(NF-1)) {print $0} }' > owh_from_hieu_SLU3_NC50.vcf_GT2_mod2

colony <- "PK_0805_NC38"

vcf_gt <- read.table(paste("1224samples_BiSNPs_MAC_DP_Fmissing_HETall_",colony,"_ID_GQfilt_infonly_queen.vcf_GT_mod",sep=""),header=TRUE)
recomb <- read.table(paste("1224samples_BiSNPs_MAC_DP_Fmissing_HETall_",colony,"_ID_GQfilt_infonly_queen_yapp_recombinations.txt",sep=""),header=TRUE) #parent sex offspring chrom left right

num_drones <- ncol(vcf_gt)-4

vcf_gt$count_hplt1 <- rep(0,nrow(vcf_gt))
vcf_gt$count_hplt2 <- rep(0,nrow(vcf_gt))
vcf_gt$count_NA <- rep(0,nrow(vcf_gt))

for (i in 1:nrow(vcf_gt)){
  num1 <- 0
  num2 <- 0
  numNA <- 0
  for (ii in 1:num_drones){
    if (is.na(vcf_gt[i,(ii+2)])){
      numNA <- numNA+1
    }
    else {
      if (vcf_gt[i,(ii+2)]==1){
        num1 <- num1+1
      }
      else if (vcf_gt[i,(ii+2)]==2){
        num2 <- num2+1
      }
    }
  }
  vcf_gt$count_hplt1[i] <- num1
  vcf_gt$count_hplt2[i] <- num2
  vcf_gt$count_NA[i] <- numNA
}

for (i in 1:num_drones) {
  vcf_gt[,(i+2)] <- vcf_gt[,(i+2)] + 2*(i-1)
}

# Plot all drone genotypes
col_vec <- c('grey57','cadetblue3','aquamarine4', 'darkgreen', 'olivedrab3', 'khaki','darkgoldenrod1', 'tomato', 'red3', 'violetred4','darkmagenta','deeppink2', 'plum3','deepskyblue1','slateblue','navy')
col_vec <- c(col_vec,col_vec)
x_min <- min(vcf_gt$POS,na.rm=TRUE) 
x_max <- max(vcf_gt$POS,na.rm=TRUE) 
plot_type <- "b"

#layout(matrix(c(1,1,1,1,1,1,2,2),nrow=4,ncol=2,byrow=TRUE)) #if plotting extra stats, last four lines

# set "add" to another value if the samples in the colony need to be divided into two plots 
# (then one has to run the script twice: first with add=0 and num_drones=(tot_num_drones/2) and then with add=(tot_num_drones/2) and num_drones=tot_num_drones)
add <- 0 

plot(vcf_gt$POS, vcf_gt[,(3+add)], type=plot_type, lty=3, pch=20, cex=1, lwd=1.5, ylim=c(add*2,(1+num_drones*2)), xlim=c(x_min,x_max), col=col_vec[(1+add)], xlab="Chr pos", ylab="Haplotypes relative to the queen; Black=inferred recomb events",axes=FALSE)
recomb_indv <- subset(recomb,offspring==colnames(vcf_gt)[3+add])
if (nrow(recomb_indv)>0){
  for (j in 1:nrow(recomb_indv)){
    lines(c(recomb_indv$left[j],recomb_indv$right[j]),rep((1.5+(add*2)),2),lty=1,col="black",lwd=1.6)
  }
  points(recomb_indv$right,rep((1.5+(add*2)),nrow(recomb_indv)),pch="|",col="black",cex=1.3)
}
for (i in (2+add):num_drones){ 
  lines(vcf_gt$POS, vcf_gt[,(i+2)], type=plot_type, lty=3, pch=20, cex=1, lwd=1.5, col=col_vec[i])
  recomb_indv <- subset(recomb,offspring==colnames(vcf_gt)[(i+2)])
  if (nrow(recomb_indv)>0){
    for (j in 1:nrow(recomb_indv)){
      lines(c(recomb_indv$left[j],recomb_indv$right[j]),rep((2*i-0.5),2),lty=1,col="black",lwd=1.6)
    }
    points(recomb_indv$right,rep((2*i-0.5),nrow(recomb_indv)),pch="|",col="black",cex=1.3)
  }
}
axis(1)
#axis(2, cex.axis=1)
text(x=x_min-(0.01*x_max), y=(2*((1+add):num_drones)-0.5), labels = colnames(vcf_gt)[(3+add):(num_drones+2)], cex=1)
title(paste(colony,": Haplotype for each drone, yapp 0.4.2",sep=""),cex.main=1.1)

# plot(vcf_gt$POS,vcf_gt$count_hplt1,type="l",lty=1,lwd=2,col=alpha.col(col="grey",alpha=0.8),ylim=c(0,num_drones),xlab="Chr pos",ylab="Num drones per hplt; Red=num miss. if >2")
# lines(vcf_gt$POS,vcf_gt$count_hplt2,lty=1,lwd=2,col=alpha.col(col="burlywood",alpha=0.8))
# vcf_highNA <- subset(vcf_gt,count_NA>2)
# points(vcf_highNA$POS,vcf_highNA$count_NA,pch=20,cex=0.8,col=alpha.col(col="red3",alpha=0.8))
