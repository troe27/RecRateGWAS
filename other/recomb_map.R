rm(list=ls())

# Use the recombinations.txt output from YAPP in order to calculate the recombination rate per bp, summed across all drones

chrom <- "NC_037653.1" #repeat for each chrom

out_file <- paste(chrom,"_CO_per_SNP.txt",sep="")

snp_pos_list <- read.table(paste("chr_pos_",chrom,".txt",sep=""),header=FALSE) # chr=V1 pos=V2 (list of positions of all SNPs included in YAPP analysis)
snp_pos_list$start_pos <- c(0,snp_pos_list$V2[-nrow(snp_pos_list)])
snp_pos_list$CO_per_bp <- rep(0,nrow(snp_pos_list))
snp_pos_list$IDX <- 1:nrow(snp_pos_list)

in_file000 <- read.table(paste(chrom,"_ID_queen_biallelic_yapp_recombinations.txt",sep=""),header=TRUE) # parent sex offspring chrom left right

to_exclude <- read.table("drones_GWAS/final_dataset/excluded_samp.txt",sep="\t",header=TRUE)
in_file00 <- subset(in_file000,!(offspring %in% to_exclude$indv))
rm(in_file000)

#Outlier if any drone in colony has sum CO:s =<10 or >=100
outliers <- c("HH_VB1_queen","IA_OS27_queen","KT_SA5_queen","UBF_J12_queen","UBF_J19_queen","UBF_J1_queen","UBF_J20_queen","UBF_J2_queen","UBF_J4_queen","UB_2_queen")

in_file0 <- subset(in_file00, !(parent %in% outliers))
rm(in_file00)

list_all_indv <- sort(unique(in_file0$offspring))
num_indv <- length(list_all_indv)

in_file <- in_file0[order(in_file0$left),]
rm(in_file0)

# Calculate the recombination rate (CO/bp) from the sum of all overlapping COs from the YAPP output
for (i in 1:nrow(in_file)){
  recomb_start <- in_file$left[i]
  recomb_end <- in_file$right[i]
  recomb_dist_bp <- recomb_end - recomb_start
  to_update <- subset(snp_pos_list, start_pos >= recomb_start & V2 <= recomb_end)
  if (nrow(to_update>0)){
    dist_bp <- to_update[nrow(to_update),"V2"] - to_update[1,"start_pos"]
    if (dist_bp != recomb_dist_bp){
      stop("Incorrect dist")
    }
    snp_pos_list$CO_per_bp[to_update$IDX] <- snp_pos_list$CO_per_bp[to_update$IDX] + (1/dist_bp)
  }
}

#write.table(snp_pos_list[,c("V1","start_pos","V2","CO_per_bp")],file=out_file,append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

#multiply CO/bp with 10^8 and divide by the number of drones in order to get cM/Mb
plot((snp_pos_list$start_pos+snp_pos_list$V2)/2, (1/length(unique(in_file$offspring)))*(10^8)*snp_pos_list$CO_per_bp,type="b",pch=20,lty=1,lwd=0.2,cex=0.2,ylab="cM/Mb",xlab=paste(chrom,", pos (bp)"))

