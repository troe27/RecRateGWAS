rm(list=ls())

library(ggplot2)

gwas_data0 <- read.table("Demetris_lmm.assoc.txt",header=TRUE)

gwas_data <- gwas_data0[,c("chr","ps","p_lrt")]

rm(gwas_data0)

# Get total length of each contig
contig_lengths <- read.table("contigs_lengths.txt",header=TRUE,row.names=1)
num_contigs <- nrow(contig_lengths)
contig_lengths$cum <- rep(0,num_contigs)
contig_lengths$cum_mid <- rep(0,num_contigs)
contig_lengths$cum_mid[1] <- 0.5*contig_lengths$Length[1]
contigs_list <- rownames(contig_lengths)

# Calcualte the cumulative length of the contigs (based on their order in the file)
# for the start position of each window and the middle position
for (i in 2:num_contigs){
  contig_lengths$cum[i] <- sum(contig_lengths$Length[1:(i-1)])
  contig_lengths$cum_mid[i] <- contig_lengths$cum[i] + 0.5*contig_lengths$Length[i]
}


# Calculate the total position based on the chromosomal position and the lengths of the contigs before
gwas_data$tot_pos <- gwas_data$ps
gwas_data$IDX <- 1:nrow(gwas_data)
for (contig in contigs_list){
  gwas_to_update <- subset(gwas_data, chr==contig, select=IDX)
  gwas_data[gwas_to_update$IDX,"tot_pos"] <- gwas_data[gwas_to_update$IDX,"tot_pos"] + contig_lengths[contig,"cum"]
}

gwas_data$p_lrt_log <- (-1)*log10(gwas_data$p_lrt)

gwas_data$pt_size <- 0.9

#***** To reduce the size of the manhattan plot: use a subset of the data points (keep all the points with low p-values)
gwas_part1 <- subset(gwas_data,p_lrt_log < 0.75)
gwas_part2 <- subset(gwas_data,p_lrt_log >= 0.75)

gwas_part1_subset <- gwas_part1[seq(1,nrow(gwas_part1),by=5),]
gwas_part1_subset$pt_size <- 1.5

gwas_pruned <- rbind(gwas_part1_subset,gwas_part2)
gwas_order <- order(gwas_pruned$tot_pos)
gwas_pruned_sorted <- gwas_pruned[gwas_order,]

#***************

col1 <- "#A56066"
col2 <- "mistyrose3"

# Manhattan plot
gwas_plot <- ggplot(gwas_pruned_sorted,aes(tot_pos,p_lrt_log)) + 
  geom_point(aes(colour=as.factor(chr),shape="."),size=gwas_pruned_sorted$pt_size) +
  scale_color_manual(values=rep(c(col1,col2),num_contigs)) +
  scale_x_continuous(label=1:16, breaks=contig_lengths$cum_mid, expand=expansion(mult=0.02,add=0)) +
  xlab("Chromosome") +
  ylab("-log10(p)") +
  theme(axis.text.x=element_text(angle = 0), #45
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70",linewidth=0.005),
        panel.grid.minor.y=element_blank(), #element_line(colour="black",linewidth=0.05),
        legend.position='none',
        axis.text=element_text(size=13),
        axis.title=element_text(size=15),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),units="cm"),
        plot.title=element_text(margin=margin(t=10,r=0,b=20,l=0)),
        axis.title.y=element_text(margin=margin(t=0,r=20,b=0,l=10))) #+
#geom_vline(xintercept=contig_lengths$cum[-1],linetype="dashed",colour="black",linewidth=0.2)

gwas_plot


