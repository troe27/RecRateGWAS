##Continue using data generated from Intra-shuffling R Script
library(xoi)
test_cM_Mb_per_SNP_all <- read.delim("C:/Users/Admin/Downloads/test_cM_Mb_per_SNP_all.txt") ##Read Demetris's map

test_cM_Mb_per_SNP_all <- test_cM_Mb_per_SNP_all %>% ##Make columns for the size of blocks between COs in bp and cM, and adding cumulative cM for blocks overtime
  group_by(chrom) %>%
  mutate(
    segment_bp = end_pos - start_pos,
    delta_cM = cM.Mb * (segment_bp / 1e6),
    cum_cM = cumsum(delta_cM)
  )

test_cM_Mb_per_SNP_all_condensed <- test_cM_Mb_per_SNP_all %>% ##Condense the table down for sparser data rows for easier interpolation of CO position downstream
  arrange(chrom, start_pos) %>%
  group_by(chrom) %>%
  mutate(group = cumsum(
    CO_per_bp != lag(CO_per_bp, default = first(CO_per_bp)) |
      cM.Mb != lag(cM.Mb, default = first(cM.Mb)) |
      start_pos != lag(end_pos, default = first(start_pos))
  )) %>%
  group_by(chrom, group, CO_per_bp, cM.Mb) %>%
  summarise(
    start_pos = min(start_pos),
    end_pos = max(end_pos),
    .groups = "drop"
  ) %>%
  select(chrom, start_pos, end_pos, CO_per_bp, cM.Mb)

test_cM_Mb_per_SNP_all_condensed <- test_cM_Mb_per_SNP_all_condensed %>% ##Finally calculate the chromosome size through sum of cumulative cM of blocks
  group_by(chrom) %>%
  mutate(
    segment_bp = end_pos - start_pos,
    delta_cM = cM.Mb * (segment_bp / 1e6),
    cum_cM = cumsum(delta_cM)
  )

test_length2 <- co_data %>% ##Interpolation of CO position in cM
  group_by(chrom) %>%
  mutate(
    cM_pos = {
      chr <- unique(chrom)
      chr_data <- test_cM_Mb_per_SNP_all_condensed %>%
        filter(chrom == chr) %>%
        arrange(start_pos)
      
      x <- chr_data$start_pos
      y <- chr_data$cum_cM
      
      fn <- splinefun(x, y, method = "monoH.FC")   ##Monotonic cubic (ensures that if the base-pair position increases, the position in cM also has to increase)
      
      vals <- fn(midpoint)
      
      vals <- pmax(min(y, na.rm = TRUE), pmin(max(y, na.rm = TRUE), vals))
      vals
    },
    L_cM = max(test_cM_Mb_per_SNP_all_condensed$cum_cM[test_cM_Mb_per_SNP_all_condensed$chrom == unique(chrom)], na.rm = TRUE)
  ) %>% ungroup()

ggplot(test_length2 %>% mutate(midpoint_Mb = midpoint / 1e6), aes(x = midpoint_Mb, y = cM_pos)) + ##Plot out position in basepair vs cM for gradient of quality (make sure the mapping is continuous)
  geom_point(alpha = 0.5, size = 0.8) +
  geom_line(alpha = 0.7) +  # optional: smoother trajectory
  facet_wrap(~ chrom, scales = "free_x", ncol = 4) +
  labs(
    x = "Physical position (Mb)",
    y = "Genetic position (cM)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

xoloc_results <- test_length2 %>% ##Make a table with an additional column called xoloc that contains lists of crossover locations (in cM), each component being a vector of locations for a different drone (necessary format for fitStahl).
  group_by(parent, chrom) %>%
  summarise(
    xoloc = list(
      group_by(cur_data(), offspring) %>%
        summarise(x = list(sort(unique(cM_pos))), .groups = "drop") %>%
        pull(x)
    ),
    L_cM = first(L_cM),
    .groups = "drop"
  )

new_xoloc_results <- xoloc_results %>%
  mutate(colony = sub("_queen$", "", parent)) %>%   ##Remove "_queen" due to slightly different naming between datasets
  semi_join(CO_corrected_bp_mean_var_per_colony_v2, 
            by = c("colony" = "colony"))

##Make new empty object
test_xoloc_results <- tibble()

all_parents <- unique(new_xoloc_results$parent)

##Define size for running the loop below
chunk_size <- 30

##Looping through parents based on the chunk size above to perform the interference analysis
for (i in seq(1, length(all_parents), by = chunk_size)) {
  
  parent_chunk <- all_parents[i:min(i + chunk_size - 1, length(all_parents))]
  
  cat("Processing parents", i, "to", min(i + chunk_size - 1, length(all_parents)), "\n")
  
  for (p in parent_chunk) {
    
    cat("Processing parent:", p, "\n")
    
    new_results <- new_xoloc_results %>% 
      filter(parent == p) %>%  
      mutate(
        fit = map2(xoloc, L_cM, ~ {
          if (is.na(.y) || length(.x) == 0) return(NA)
          attr(.x, "L") <- .y
          tryCatch(fitStahl(.x), error = function(e) NA)
        })
      )
    
    test_xoloc_results <- bind_rows(test_xoloc_results, new_results)
  }
  
  cat("Finished chunk", i, "to", min(i + chunk_size - 1, length(all_parents)), "\n")
  
  ##Pause for inspection; adjust chunk size if want to run more/all samples per run
  cat("PAUSE: Check output for this chunk before continuing...\n")
  readline(prompt = "Press [Enter] to continue to next chunk...")
}

test_xoloc_results <- test_xoloc_results %>% ##Extract the number of COs, and results from the interference analysis
  mutate(
    nCO = map_int(xoloc, ~ sum(lengths(.x))),
    nu = map_dbl(fit, ~ .x["nu"]),
    p  = map_dbl(fit, ~ .x["p"]),
    loglik = map_dbl(fit, ~ .x["loglik"]),
    nu0 = map_dbl(fit, ~ .x["nu0"]),
    loglik0 = map_dbl(fit, ~ .x["loglik0"]),
    ln_testing = map_dbl(fit, ~ .x["ln LR testing p=0"])
  )

lrt_results <- test_xoloc_results %>% ##Use a likelihood-ratio test to determine the validity of the nu and p values (loglik is maximixed log likelihood, and loglik0 is maximized log likelihood if p=0); Currently unsure if this is how it is supposed to be done. 
  mutate(
    LR_stat = 2 * (loglik - loglik0),
    LR_pvalue = pchisq(LR_stat, df = 1, lower.tail = FALSE),
    ##Classify significance
    p_signif = case_when(
      is.na(p) ~ "fit_failed",
      p == 0 & (LR_pvalue > 0.05 | is.na(LR_pvalue)) ~ "p=0 and significant - Support for all interference CO",
      p == 0 & LR_pvalue <= 0.05 ~ "p=0 and not significant - No support for all interference CO",
      p > 0 & LR_pvalue > 0.05 ~ "p>0, but not significant - No support for non-interference CO",
      p > 0 & LR_pvalue <= 0.05 ~ "p>0, but significant - Support for non-interference CO",
      TRUE ~ "other"
    )
  )

##Summary of the above table
lrt_summary <- lrt_results %>%
  count(p_signif, name = "n_chromosomes") %>%
  mutate(proportion = n_chromosomes / sum(n_chromosomes))
print(lrt_summary)

mean(test_xoloc_results$nu)
mean(test_xoloc_results$p)
cor.test(test_xoloc_results$nu, test_xoloc_results$nCO)
cor.test(test_xoloc_results$p, test_xoloc_results$nCO)

chromosome_summary <- test_xoloc_results %>% ##Summary of fitStahl results per chromosome
  group_by(chrom) %>%
  summarise(
    mean_nCO = mean(nCO, na.rm = TRUE),
    mean_L_cM = mean(L_cM, na.rm = TRUE),
    mean_nu = mean(nu, na.rm = TRUE),
    mean_p = mean(p, na.rm = TRUE),
  ) %>%
  rename(
    Chromosome = chrom,
    `Mean number of COs` = mean_nCO,
    `Genetic length (cM)` = mean_L_cM,
    `Mean interference index (ν)` = mean_nu,
    `Mean proportion of non-interfering crossovers (p)` = mean_p
  )

ggplot(test_xoloc_results, aes(x = nCO, y = nu)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  labs(
    title = "Interference Index (ν) vs. Total Number of Crossovers",
    x = "Total crossovers (per chromosome)",
    y = "Nu"
  ) +
  theme_minimal()

ggplot(test_xoloc_results, aes(x = nCO, y = p)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  labs(
    title = "Proportion of non-interference CO (p) vs. Total Number of Crossovers",
    x = "Total crossovers (per chromosome)",
    y = "P"
  ) +
  theme_minimal()
