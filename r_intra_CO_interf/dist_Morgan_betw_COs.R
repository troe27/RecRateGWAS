rm(list=ls())

# For input to CO interference calculations
# For each CO in each indv, count the number of COs in the whole dataset that are located between this CO and the previous CO in the same indv (or the start of the chromosome)

# output from YAPP, concatenated for all chromosomes
CO_all_concat <- read.table("concat_ID_queen_biallelic_yapp_recombinations_filtered.txt",header=TRUE) # parent  sex     offspring       chrom   left    right
CO_all_concat$midpoint <- (CO_all_concat$left + CO_all_concat$right)/2

# CO_all_concat is sorted on 1: chromosome, 2: individual and 3: position
# Make a version of this table that is sorted only by CO-position instead (mid of left and right) - used later for counting intervening COs

CO_pos_order <- order(CO_all_concat$chrom,CO_all_concat$midpoint)
CO_all_sorted_pos <- CO_all_concat[CO_pos_order,]

CO_all_concat$cM_to_previous_CO <- rep(0,nrow(CO_all_concat))

num_drones <- length(unique(sort(CO_all_concat$offspring)))

for (i in 1:nrow(CO_all_concat)){
  chr_i <- CO_all_concat$chrom[i]
  indv_i <- CO_all_concat$offspring[i]
  bp_i <- CO_all_concat$midpoint[i]
  # check if there is a CO before on the same chrom in the same indv
  # then count the number of COs after that CO
  # otherwise count the number of COs from the start of the chrom
  if (i>1){
    if (chr_i == CO_all_concat$chrom[(i-1)] & indv_i == CO_all_concat$offspring[(i-1)]){
      bp_before <- CO_all_concat$midpoint[(i-1)]
      #print(bp_before)
    }
    else {
      bp_before <- 0
      #print(bp_before)
    }
  }
  else {
    bp_before <- 0
  }
  CO_all_concat$cM_to_previous_CO[i] <- nrow(subset(CO_all_sorted_pos, CO_all_sorted_pos$chrom==chr_i & CO_all_sorted_pos$midpoint>bp_before & CO_all_sorted_pos$midpoint<bp_i))*100/num_drones
}

#write.table(CO_all_concat,file="drones_GWAS/final_dataset/r_intra_CO_interference/concat_ID_queen_biallelic_yapp_recombinations_filtered_CO_dist_cM.txt",col.names=TRUE,row.names = FALSE,quote=FALSE,sep="\t",append=FALSE)

