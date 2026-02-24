rm(list=ls())

eig_vec0 <- read.table("plink2_hieu.eigenvec",header=TRUE)
eig_val <- read.table("plink2_hieu.eigenval",header=FALSE)

to_exclude <- read.table("drones_GWAS/final_dataset/excluded_samp.txt",sep="\t",header=TRUE)

eig_vec <- subset(eig_vec0, !(IID %in% to_exclude$indv))

colony_list <- unique(eig_vec$FID)

pc1_expl <- eig_val$V1[1]/sum(eig_val$V1,na.rm=TRUE)
pc2_expl <- eig_val$V1[2]/sum(eig_val$V1,na.rm=TRUE)

eig_vec$IDX <- 1:nrow(eig_vec)

color_list <- rep("gray50",nrow(eig_vec))

colors_repeat <- c("yellow4","cadetblue3","lightpink4")
colors_repeat <- rep(colors_repeat,length.out=length(colony_list))

# assign different colors per colony
i <- 1
for (colony in colony_list){
  to_update <- subset(eig_vec,FID==colony)
  color_list[to_update$IDX] <- colors_repeat[i]
  i <- i+1
}


# to_update <- subset(eig_vec,FID=="Cgroup")
# color_list[to_update$IDX] <- "red2"
# 
# to_update <- subset(eig_vec,FID=="Mgroup")
# color_list[to_update$IDX] <- "slateblue"

par(mfrow=c(1,1))

plot(eig_vec$PC1,eig_vec$PC2,col=color_list,type="p",pch=20,cex=0.5,xlab=paste("PC1 = ",round(pc1_expl*100,digits=2)," %",sep=""),ylab=paste("PC2 = ",round(pc2_expl*100,digits=2)," %",sep=""))
to_update <- subset(eig_vec,FID=="Cgroup")
points(to_update$PC1,to_update$PC2,type="p",pch=20,cex=0.7,col="orange2")
to_update <- subset(eig_vec,FID=="Mgroup")
points(to_update$PC1,to_update$PC2,type="p",pch=20,cex=0.7,col="slateblue")

# ******** Recomb rate vs PC1

recomb_data00 <- read.table("filt_input_YAPP/YAPP_output/all_chr_all_drones_CO_per_bp_corrected_genome_size2.txt",header=TRUE)
recomb_data0 <- subset(recomb_data00,!(indv %in% to_exclude$indv))

recomb_data <- subset(recomb_data0,num_CO<100) # remove outliers
recomb_data$cM_Mb <- (10^8)*recomb_data$CO_per_bp

recomb_pca_merge <- merge(eig_vec,recomb_data,by.x="IID",by.y="indv")

colony_list_filt_uniq <- unique(recomb_pca_merge$FID)

recomb_pca_colony <- data.frame("colony"=colony_list_filt_uniq,"mean_PC1"=rep(0,length(colony_list_filt_uniq)),"mean_cMMb"=rep(0,length(colony_list_filt_uniq)),"mean_num_CO"=rep(0,length(colony_list_filt_uniq)),"mean_corr_len"=rep(0,length(colony_list_filt_uniq)))

for (ii in 1:length(colony_list_filt_uniq)){
  colony_ii <- recomb_pca_colony$colony[ii]
  data_subset <- subset(recomb_pca_merge,FID==colony_ii)
  recomb_pca_colony$mean_PC1[ii] <- mean(data_subset$PC1,na.rm=TRUE)
  recomb_pca_colony$mean_cMMb[ii] <- mean(data_subset$cM_Mb,na.rm=TRUE)
  recomb_pca_colony$mean_num_CO[ii] <- mean(data_subset$num_CO,na.rm=TRUE)
  recomb_pca_colony$mean_corr_len[ii] <- mean(data_subset$corrected_genome_bp,na.rm=TRUE)
}

par(mfrow=c(1,1))

cex_val_ax <- 1.5
cex_val_lab <- 1.5

# plot per colony
plot(recomb_pca_colony$mean_PC1,recomb_pca_colony$mean_cMMb,type="p",pch=20,xlab="PC1",ylab="cM/Mb",cex=0.9,cex.axis=cex_val_ax,cex.lab=cex_val_lab)
pc_corr <- cor.test(recomb_pca_colony$mean_PC1,recomb_pca_colony$mean_cMMb,alternative="two.sided",method="spearman")

# plot per drone
# plot(recomb_pca_merge$PC1,recomb_pca_merge$cM_Mb,type="p",pch=20,xlab="PC1",ylab="cM/Mb",cex=0.6)
# pc_corr <- cor.test(recomb_pca_merge$PC1,recomb_pca_merge$cM_Mb,alternative="two.sided",method="spearman")




