rm(list=ls())
library(scales)

to_exclude <- read.table("HoneyBees/drones_GWAS/final_dataset/excluded_samp.txt",sep="\t",header=TRUE)

recomb_data00 <- read.table("all_chr_all_drones_CO_per_bp_corrected_genome_size2.txt",header=TRUE)

recomb_data0 <- subset(recomb_data00,!(indv %in% to_exclude$indv))

# histogram including outliers
par(mfrow=c(1,1))
hist(recomb_data0$num_CO,col="lightpink3",breaks=150,xlab="Number of CO:s per drone",main="")
abline(v=100,lty=2)

recomb_data <- subset(recomb_data0,num_CO<100) # remove outliers

# histograms excluding outliers
par(mfrow=c(1,3))
q_list1 <- unique(recomb_data$colony)
drones_per_colony <- rep(0,length(q_list1))
rr_per_colony <- rep(0,length(q_list1))
CO_per_drone_per_colony <- rep(0,length(q_list1))

for (qi in 1:length(q_list1)){
  recomb_qi <- subset(recomb_data,colony==q_list1[qi])
  drones_per_colony[qi] <- nrow(recomb_qi)
  recomb_qi$cMMb <- recomb_qi$CO_per_bp*(10^8)
  rr_per_colony[qi] <- mean(recomb_qi$cMMb,na.rm=TRUE)
  CO_per_drone_per_colony[qi] <- mean(recomb_qi$num_CO,na.rm=TRUE)
}
hist(drones_per_colony,col="lightpink3",breaks=(4:20)+0.5,include.lowest = FALSE,xlab="Number of drones per colony",main="")
hist(recomb_data$num_CO,col="lightpink3",breaks=20,xlab="Number of CO:s per drone",main="")
#hist(recomb_data$CO_per_bp*(10^8),col="lightpink3",xlab="Recomb. rate per drone (cM/Mb)",main="")
hist(rr_per_colony,col="lightpink3",breaks=20,xlab="Mean recomb. rate per colony (cM/Mb)",main="")
#hist(CO_per_drone_per_colony,col="lightpink3",breaks=15,xlab="Mean number of COs per drone per colony",main="")

# broad sense heritability based on ANOVA
aov_lm <- anova(lm(CO_per_bp ~ colony, data=recomb_data))

h2 <- aov_lm$`Sum Sq`[1]/(aov_lm$`Sum Sq`[1] + aov_lm$`Sum Sq`[2])


# Age analysis

age_list <- read.table("/Users/turev988/Sync/HoneyBees/drones_GWAS/age_list.txt",header=TRUE)

age_list$meanCO <- rep(0,nrow(age_list))
age_list$varCO <- rep(0,nrow(age_list))

for (i in 1:nrow(age_list)){
  colony_i <- age_list$Colony_ID[i]
  rec_colony <- subset(recomb_data,colony==colony_i)
  if (nrow(rec_colony)>0){
    age_list$meanCO[i] <- mean(rec_colony$CO_per_bp,na.rm=TRUE)
    age_list$varCO[i] <- var(rec_colony$CO_per_bp,na.rm=TRUE)
  }
  else {
    age_list$meanCO[i] <- NA
    age_list$varCO[i] <- NA
  }
}

age_meanCO_cor_obs <- cor.test(age_list$Queen_age_inferred,age_list$meanCO,alternative="two.sided",na.rm=TRUE,method="spearman")
age_varCO_cor_obs <- cor.test(age_list$Queen_age_inferred,age_list$varCO,alternative="two.sided",na.rm=TRUE,method="spearman")

age_list$Queen_age_inferred_factor <- as.factor(age_list$Queen_age_inferred)

age_meanCO_aov_obs <- anova(lm(meanCO ~ Queen_age_inferred_factor, data=age_list))
age_varCO_aov_obs <- anova(lm(varCO ~ Queen_age_inferred_factor, data=age_list))

# Random permutations of age labels
niter <- 100

age_meanCO_cor_rand <- rep(0,niter)
age_varCO_cor_rand <- rep(0,niter)
age_meanCO_aov_rand <- rep(0,niter)
age_varCO_aov_rand <- rep(0,niter)

age_list$Queen_age_rand <- age_list$Queen_age_inferred

for (i in 1:niter){
  age_list$Queen_age_rand <- sample(age_list$Queen_age_inferred,size=nrow(age_list),replace=FALSE)
  age_meanCO_cor <- cor.test(age_list$Queen_age_rand,age_list$meanCO,alternative="two.sided",na.rm=TRUE,method="spearman")
  age_meanCO_cor_rand[i] <- age_meanCO_cor$estimate
  age_varCO_cor <- cor.test(age_list$Queen_age_rand,age_list$varCO,alternative="two.sided",na.rm=TRUE,method="spearman")
  age_varCO_cor_rand[i] <- age_varCO_cor$estimate
  age_list$Queen_age_rand_factor <- as.factor(age_list$Queen_age_rand)
  age_meanCO_aov <- anova(lm(meanCO ~ Queen_age_rand_factor, data=age_list))
  age_meanCO_aov_rand[i] <- age_meanCO_aov$`F value`[1]
  age_varCO_aov <- anova(lm(varCO ~ Queen_age_rand_factor, data=age_list))
  age_varCO_aov_rand[i] <- age_varCO_aov$`F value`[1]
}

# calculate p-values by comparison of results from random permutations and the observed data

p_cor_meanCO <- length(subset(age_meanCO_cor_rand,age_meanCO_cor_rand>=abs(age_meanCO_cor_obs$estimate) | age_meanCO_cor_rand<=-(abs(age_meanCO_cor_obs$estimate))))/niter
p_cor_varCO <- length(subset(age_varCO_cor_rand,age_varCO_cor_rand>=abs(age_varCO_cor_obs$estimate) | age_varCO_cor_rand<=-(abs(age_varCO_cor_obs$estimate))))/niter

p_aov_meanCO <- length(subset(age_meanCO_aov_rand,age_meanCO_aov_rand>=age_meanCO_aov_obs$`F value`[1]))/niter
p_aov_varCO <- length(subset(age_varCO_aov_rand,age_varCO_aov_rand>=age_varCO_aov_obs$`F value`[1]))/niter


par(mfrow=c(2,2))

hist(age_meanCO_cor_rand)
abline(v=age_meanCO_cor_obs$estimate,lty=2,col="red")
hist(age_varCO_cor_rand)
abline(v=age_varCO_cor_obs$estimate,lty=2,col="red")

hist(age_meanCO_aov_rand)
abline(v=age_meanCO_aov_obs$`F value`[1],lty=2,col="red")
hist(age_varCO_aov_rand)
abline(v=age_varCO_aov_obs$`F value`[1],lty=2,col="red")


# ****************** Boxplot of recomb rate per colony and q age

age_drones_merge <- merge(age_list,recomb_data,by.x="Colony_ID",by.y="colony",all.y=TRUE)

age_drones_merge_sorted <- age_drones_merge[order(age_drones_merge$Queen_age_inferred),]
age_drones_merge_sorted$IDX <- 1:nrow(age_drones_merge_sorted)
col_list <- rep("gray60",length(unique(age_drones_merge_sorted$Colony_ID)))
a1 <- 0.7
col_list_uniq <- c(alpha("yellow4",a1),alpha("lightpink4",a1),alpha("cadetblue3",a1),alpha("orange3",a1),alpha("slateblue",a1))

i <- 1
for (colony_i in unique(age_drones_merge_sorted$Colony_ID)){
  colony_subset <- subset(age_drones_merge_sorted, Colony_ID==colony_i)
  y <- colony_subset$Queen_age_inferred[1]
  if (!is.na(y)){
    col_list[i] <- col_list_uniq[(y+1)]
  }
  i <- i+1
}

par(mfrow=c(1,1))

boxplot((CO_per_bp)*(10^8) ~ Colony_ID, data=age_drones_merge_sorted,col=col_list,pch=21,bg=col_list,cex=0.7,names=FALSE,xaxt="n",xlab="Colony",ylab="cM/Mb",boxlwd=0.6,medlwd=0.7,whisklwd=0.6,outlwd=0.5)
legend("topright",c("0 years","1 years","2 years","3 years","4 years","unknown"),col=c(col_list_uniq,"gray60"),title="Queen age",pch=20,pt.cex=1.2,cex=0.9,inset=0.002,x.intersp=0.65,y.intersp=0.75)




