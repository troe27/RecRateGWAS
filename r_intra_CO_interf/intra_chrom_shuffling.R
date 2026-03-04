rm(list=ls())
packages <- c(
  "dplyr", "purrr", "readr", "tidyr", "stringr", "ggplot2"
) ##https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed], dependencies = TRUE)
}
invisible(lapply(packages, library, character.only = TRUE))
co_data_unfilt <- list.files(path = "drones_GWAS/final_dataset/filt_input_YAPP/YAPP_output/", pattern = "^NC_.*_yapp_recombinations.txt_RIC$", full.names = TRUE) %>% map_dfr(~ read_table(.x, col_types = cols()) %>% mutate(source_file = basename(.x))) ##Read all recombination files into one table

# samples to filter out
bad_qual <- read.table("drones_GWAS/final_dataset/excluded_samp.txt",sep="\t",header=TRUE)
outlier_samp <- read.table("drones_GWAS/final_dataset/outliers_after_qual_filt.txt",sep="\t",header=TRUE)
co_data <- subset(co_data_unfilt,!(offspring %in% c(bad_qual$indv, outlier_samp$indv)))

co_data <- co_data %>% mutate(midpoint = (left + right) / 2) ##Calculate CO positions as the average of left and right coordinates of CO blocks
chr_cols <- c("NC_037638.1","NC_037639.1","NC_037640.1","NC_037641.1","NC_037642.1","NC_037643.1","NC_037644.1","NC_037645.1","NC_037646.1","NC_037647.1","NC_037648.1","NC_037649.1","NC_037650.1","NC_037651.1","NC_037652.1","NC_037653.1")
chrom_length <- c(27754200, 16089512, 13619445, 13404451, 13896941, 17789102, 14198698, 12717210, 12354651, 12360052, 16352600, 11514234, 11279722, 10670842, 9534514, 7238532)
chroms <- tibble(chrom = chr_cols, length = chrom_length) ##Make a table of chromosome name and chromosome length in bp
chroms <- chroms %>% mutate(L = length / sum(length)) ##Add a column for the fraction of chromosome length versus genome size
co_data <- co_data %>% left_join(chroms, by = "chrom") %>% mutate(frac_pos = midpoint / length) ##Add the position of CO as fraction of the total chromosome length
total_genome_length <- sum(chroms$length)
compute_pk <- function(midpoints, chrom_len) {
  pos_frac <- sort(midpoints) / chrom_len ##Calculate CO positions as fractions of the chromosome instead
  breaks <- c(0, pos_frac, 1) ##Record the position of the COs relatively to the start and stop of the chromosome as fractions
  segs <- diff(breaks) 
  p <- sum(segs[seq(1, length(segs), by = 2)]) ##alternate segments as coming from alternating homologs
  tibble(p = p, two_p_1mp = 2 * p * (1 - p))
} #Function to calculate the 2pk(1−pk) part
r_components <- co_data %>%
  group_by(parent, offspring, chrom) %>% ##Group CO by queen, drone name, and chromosome
  summarise(
    midpoints = list(midpoint),
    n_CO = n(),  # crossover count per chromosome
    .groups = "drop"
  ) %>% ##Make a list of all CO for that chromosome, and add a column for CO count per queen per drone per chromosome
  left_join(chroms, by = "chrom") %>% ##Adds the length (bp) and the relative fraction of genome length (L).
  mutate(pk_data = map2(midpoints, length, compute_pk)) %>% ##Use the 2pk(1−pk) function
  unnest(pk_data) %>%
  mutate(contrib = two_p_1mp * L^2) ##Calculate the rest of the equation
r_intra_per_gamete <- r_components %>%
  group_by(parent, offspring) %>%
  summarise(
    r_intra = sum(contrib),
    n_CO_drone = sum(n_CO),
    .groups = "drop") ##Intra-chromosomal shuffling per drone
r_intra_per_parent <- r_intra_per_gamete %>%
  group_by(parent) %>%
  summarise(
    mean_r_intra = mean(r_intra),
    sd_r_intra = sd(r_intra),
    n_gametes = n(),
    .groups = "drop"
  ) ##Intra-chromosomal shuffling per parent as mean and sd of intra-chromosomal shuffling per drone

corr_v1 <- cor.test(r_intra_per_gamete$n_CO_drone, r_intra_per_gamete$r_intra,method="spearman",alternative="two.sided") ##Correlation between CO count and intra-shuffling
ggplot(r_intra_per_gamete, aes(x=n_CO_drone, y = r_intra)) + geom_point(size=3, alpha=0.7) + geom_smooth(method="lm", se=TRUE) + labs(x="CO Count", y=expression(r[intra])) ##Correlation plotting
