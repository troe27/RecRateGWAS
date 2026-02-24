rm(list=ls())

library(MM4LMM)
library(dplyr)

# Calculate heritability of CO-freq and intra-chromosomal genetic shuffling (r_intra - with and without controlling for CO-freq)

# Calculate heritability based on GRM from LDAK
# Commands for generating GRM:
# plink2 --vcf NC_XXX.vcf.gz --make-bed --vcf-half-call missing --out NC_XXX_queen_phased
# cat Amel_interval.list | xargs -n 1 -P 16 -I {} bash -c "ldak6 --cut-weights {} --bfile {}_queen_phased" (to do ldak6 on multiple chromosomes at the same time,)
# cat Amel_interval.list | xargs -n 1 -P 16 -I {} bash -c "ldak6 --calc-weights-all {} --bfile {}_queen_phased"
# Merge the above together.
# ldak6 --calc-kins-direct LDAKgmat --bfile all_queen_phased --weights weights.short --power -.25 --kinship-raw YES
#####

# import r_intra data, already filtered for quality and recombination outliers
r_intra_init0 <- read.table("drones_GWAS/final_dataset/r_intra_CO_interference/r_intra_per_drone_filt.txt",header=TRUE) #parent  offspring       r_intra colony_id       total_CO

# import observed recomb data
# not filtered, but filtering by combining with r_intra
obs_recomb00 <- read.table("all_chr_all_drones_CO_per_bp_corrected_genome_size2.txt",header=TRUE) #colony	indv	num_CO	indv_same	corrected_genome_bp	CO_per_bp

r_intra <- left_join(r_intra_init,obs_recomb00,by=join_by(offspring==indv)) # contains both r_intra and recomb rate

# To control for CO-freq when analysing r_intra:
# linear regression of r_intra on CO-freq, then take the residuals of that as the r_intra phenotype controlled for CO-freq. 
lm_intra_recomb <- lm(r_intra ~ CO_per_bp, data=r_intra)
r_intra$linreg_residuals <- lm_intra_recomb$residuals

num_drones <- nrow(r_intra)
queen_list <- unique(r_intra$colony_id)
num_queens <- length(queen_list)

# import GRM

grm_table <- read.table("unimputed_data/LDAKgmat.grm.raw",header=FALSE)
grm_queens_all <- read.table("unimputed_data/all_queen_phased.fam",header=FALSE) #colony names corresponding to GRM matrix

rownames(grm_table) <- grm_queens_all$V1 #rownames and colnames for relatedness plot
colnames(grm_table) <- grm_queens_all$V1
grm_table$queens <- grm_queens_all$V1
grm_table$IDX <- 1:nrow(grm_table)

# only include the same colonies in the GRM as in the (already filtered) r_intra table
grm_include <- subset(grm_table,queens %in% queen_list)

grm_table_filt <- grm_table[grm_include$IDX,grm_include$IDX]
grm_matrix <- as.matrix(grm_table_filt)
#heatmap(grm_matrix,Rowv=NA,Colv=NA,scale="none",main="LDAK")

#sort queens in r_intra in same order as in grm_include$queens

r_intra$dummy <- rep(0,nrow(r_intra))
r_intra$IDX <- 1:nrow(r_intra)
for (i in 1:nrow(grm_include)){
  queen_i <- grm_include$queens[i]
  ri_subset <- subset(r_intra,colony_id==queen_i)
  r_intra$dummy[ri_subset$IDX] <- i
}

r_intra_order <- order(r_intra$dummy)

r_intra_sorted <- r_intra[r_intra_order,]
queen_list_sorted <- unique(r_intra_sorted$colony_id)

r_intra_sorted$IDX <- 1:nrow(r_intra_sorted)

# create all input matrices for MMEst

VL_m <- diag(num_queens)
VL_a <- grm_matrix 
VL_e <- diag(num_drones)

ZL_e <- diag(num_drones)

r_intra_sorted$colony_id = factor(r_intra_sorted$colony_id,levels=queen_list_sorted)

## use model.matrix to build the design matrix with correct contrast for a random effect
ZL_m = model.matrix( ~ 0 + r_intra_sorted$colony_id)
ZL_a <- ZL_m

VL <- list(VL_m,VL_a,VL_e)
ZL <- list(ZL_m,ZL_a,ZL_e)

# Just scale the phenotype to a Gaussian distribution to avoid numerical problems
r_intra_sorted$Y_COfreq = qqnorm(r_intra_sorted$CO_per_bp)$x 
r_intra_sorted$Y_r_intra = qqnorm(r_intra_sorted$r_intra)$x 
r_intra_sorted$Y_r_intra_resid = qqnorm(r_intra_sorted$linreg_residuals)$x 

# Animal Model fitted with ML4LMM
reml_COfreq <- MMEst(Y=r_intra_sorted$Y_COfreq,VarList=VL,ZList=ZL,Method="Reml")
reml_r_intra <- MMEst(Y=r_intra_sorted$Y_r_intra,VarList=VL,ZList=ZL,Method="Reml")
reml_r_intra_resid <- MMEst(Y=r_intra_sorted$Y_r_intra_resid,VarList=VL,ZList=ZL,Method="Reml") #controlling for CO-freq

### Heritability CO-freq

sigma_m <- as.numeric(reml_COfreq$NullModel$Sigma2[1])
sigma_a <- as.numeric(reml_COfreq$NullModel$Sigma2[2])
sigma_e <- as.numeric(reml_COfreq$NullModel$Sigma2[3])

H2_broad_COfreq <- (sigma_m+sigma_a)/(sigma_m+sigma_a+sigma_e)
h2_narrow_COfreq <- (sigma_a)/(sigma_m+sigma_a+sigma_e)
perm_env_effect_COfreq = sigma_m/(sigma_m+sigma_a+sigma_e)

### Heritability r_intra

sigma_m <- as.numeric(reml_r_intra$NullModel$Sigma2[1])
sigma_a <- as.numeric(reml_r_intra$NullModel$Sigma2[2])
sigma_e <- as.numeric(reml_r_intra$NullModel$Sigma2[3])

H2_broad_r_intra <- (sigma_m+sigma_a)/(sigma_m+sigma_a+sigma_e)
h2_narrow_r_intra <- (sigma_a)/(sigma_m+sigma_a+sigma_e)
perm_env_effect_r_intra = sigma_m/(sigma_m+sigma_a+sigma_e)

### Heritability r_intra, controlling for CO-freq

sigma_m <- as.numeric(reml_r_intra_resid$NullModel$Sigma2[1])
sigma_a <- as.numeric(reml_r_intra_resid$NullModel$Sigma2[2])
sigma_e <- as.numeric(reml_r_intra_resid$NullModel$Sigma2[3])

H2_broad_r_intra_resid <- (sigma_m+sigma_a)/(sigma_m+sigma_a+sigma_e)
h2_narrow_r_intra_resid <- (sigma_a)/(sigma_m+sigma_a+sigma_e)
perm_env_effect_r_intra_resid = sigma_m/(sigma_m+sigma_a+sigma_e)

