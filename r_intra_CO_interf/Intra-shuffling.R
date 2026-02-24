rm(list=ls())
packages <- c(
  "dplyr", "purrr", "readr", "tidyr", "stringr", "ggplot2"
) ##https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed], dependencies = TRUE)
}
invisible(lapply(packages, library, character.only = TRUE))
co_data_unfilt <- list.files(path = "/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/filt_input_YAPP/YAPP_output/", pattern = "^NC_.*_yapp_recombinations.txt_RIC$", full.names = TRUE) %>% map_dfr(~ read_table(.x, col_types = cols()) %>% mutate(source_file = basename(.x))) ##Read all recombination files into one table
bad_qual <- read.table("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/excluded_samp.txt",sep="\t",header=TRUE)
outlier_samp <- read.table("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/outliers_after_qual_filt.txt",sep="\t",header=TRUE)
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
  summarise(r_intra = sum(contrib), .groups = "drop") ##Intra-chromosomal shuffling per drone
r_intra_per_parent <- r_intra_per_gamete %>%
  group_by(parent) %>%
  summarise(
    mean_r_intra = mean(r_intra),
    sd_r_intra = sd(r_intra),
    n_gametes = n(),
    .groups = "drop"
  ) ##Intra-chromosomal shuffling per parent as mean and sd of intra-chromosomal shuffling per drone
CO_corrected_bp_mean_var_per_colony_v2 <- read.delim("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/r_intra_CO_interference/CO_corrected_bp_mean_var_per_colony_v2.txt")
r_intra_per_gamete <- r_intra_per_gamete %>% mutate(colony_id = str_remove(parent, "_queen$")) ##Add a column of colony ID without the queen suffix
r_intra_filtered <- r_intra_per_gamete %>% semi_join(CO_corrected_bp_mean_var_per_colony_v2, by = c("colony_id" = "colony")) ##Remove any colony outlier based on existing filtered colonies
r_intra_filtered <- r_intra_filtered %>%
  left_join(
    CO_corrected_bp_mean_var_per_colony_v2 %>%
      select(colony, mean_num_CO_uncorr),
    by = c("colony_id" = "colony") 
  ) %>%
  mutate(total_CO = mean_num_CO_uncorr) %>%
  select(-mean_num_CO_uncorr)
r_intra_per_parent <- r_intra_filtered %>% group_by(parent) %>%
  summarise(
    mean_r_intra = mean(r_intra),
    sd_r_intra = sd(r_intra),
    n_gametes = n(),
    .groups = "drop"
  ) #Recalculate parents intra-shuffling frequency
cor.test(r_intra_filtered$total_CO, r_intra_filtered$r_intra) ##Correlation between CO count and intra-shuffling
ggplot(r_intra_filtered, aes(x=total_CO, y = r_intra)) + geom_point(size=3, alpha=0.7) + geom_smooth(method="lm", se=TRUE) + labs(x="CO Count", y=expression(r[intra])) ##Correlation plotting
