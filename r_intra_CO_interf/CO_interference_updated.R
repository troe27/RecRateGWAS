rm(list=ls())

library(xoi)

CO_all_cM_dist <- read.table("drones_GWAS/final_dataset/r_intra_CO_interference/concat_ID_queen_biallelic_yapp_recombinations_filtered_CO_dist_cM.txt",header=TRUE)
#parent  sex     offspring       chrom   left    right   midpoint        cM_to_previous_CO

# cM_to_previous_CO = 0 in some cases, leads to problems in the fitStahl function, therefore replace those values with 0.01 (minimum non-zero distance is otherwise 0.066)

for (i in 1:nrow(CO_all_cM_dist)){
  if (CO_all_cM_dist$cM_to_previous_CO[i] < 0.01){
    CO_all_cM_dist$cM_to_previous_CO[i] <- 0.01
  }
}

num_drones <- 1509

chrom_list <- unique(CO_all_cM_dist$chrom)
num_chrom <- length(chrom_list)
chrom_cM_lengths <- rep(0,num_chrom)

for (i in 1:num_chrom){
  chrom_i <- chrom_list[i]
  chrom_cM_lengths[i] <- nrow(subset(CO_all_cM_dist,chrom==chrom_i))*100/num_drones
}


# for each queen, make list of (num_drones_in_colony)*16 vectors of CO positions
# make corresponding vector of chrom lengths
# run fitStahl for each queen

queens_all <- unique(CO_all_cM_dist$parent)
num_queens <- length(queens_all)

fitstahl_results <- data.frame(queen=queens_all,nu=rep(0,num_queens),p_nonint=rep(0,num_queens),log_lik=rep(0,num_queens),nu0=rep(0,num_queens),log_lik0=rep(0,num_queens),log_lik_ratio=rep(0,num_queens))

for (i in 1:num_queens){
  queen_i <- queens_all[i]
  CO_cM_queen_i <- subset(CO_all_cM_dist, parent==queen_i)
  drones_colony <- unique(CO_cM_queen_i$offspring)
  num_drones_colony <- length(drones_colony)
  chrom_lengths_per_queen <- rep(0,num_chrom*num_drones_colony)
  COs_per_queen <- list()
  COs_per_queen[[(num_chrom*num_drones_colony+1)]] <- c(0,0,0) # dummy entry to make list of the right length
  print(length(COs_per_queen))
  k <- 1
  for (j in 1:num_drones_colony){
    drone_j <- drones_colony[j]
    CO_cM_drone_j <- subset(CO_cM_queen_i,offspring==drone_j)
    for (ij in 1:num_chrom){
      chrom_ij <- chrom_list[ij]
      CO_cM_drone_chrom_ij <- subset(CO_cM_drone_j,chrom==chrom_ij)
      if (nrow(CO_cM_drone_chrom_ij)>0){
        cum_dists <- as.vector(CO_cM_drone_chrom_ij$cM_to_previous_CO)
        if (length(cum_dists)>1){
          for (ii in 2:length(cum_dists)){
            cum_dists[ii] <- cum_dists[ii] + cum_dists[(ii-1)]
          }
        }
        COs_per_queen[[k]] <- cum_dists
      }
      else {
        print("No COs")
        # if no COs: COs_per_queen[[k]] left empty, only chrom_lengths_per_queen[k] updated
      }
      chrom_lengths_per_queen[k] <- chrom_cM_lengths[ij]
      k <- k+1
    }
  }
  if (length(COs_per_queen)<length(chrom_lengths_per_queen)){
    print("!!!!!!!!!!!!!!!!!Error: Different lengths**************")
  }
  CO_interf <- fitStahl(COs_per_queen[1:(k-1)],chrlen=chrom_lengths_per_queen,verbose=TRUE)
  fitstahl_results$nu[i] <- CO_interf["nu"]
  fitstahl_results$p_nonint[i] <- CO_interf["p"]
  fitstahl_results$log_lik[i] <- CO_interf["loglik"]
  fitstahl_results$nu0[i] <- CO_interf["nu0"]
  fitstahl_results$log_lik0[i] <- CO_interf["loglik0"]
  fitstahl_results$log_lik_ratio[i] <- CO_interf["ln\ LR\ testing\ p=0"]
}

#write.table(fitstahl_results, file="fitstahl_results.txt",col.names = TRUE,row.names = FALSE, quote=FALSE, sep="\t",append=FALSE)
