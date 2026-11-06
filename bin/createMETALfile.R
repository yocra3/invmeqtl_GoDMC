arguments <- commandArgs(T)

results_Rdata <- arguments[1]
inv_sum <- arguments[2]
prefix <- arguments[3]

library(tidyverse)
## Load results
load(results_Rdata)
eqtls <-  me$all$eqtls

## Load sample sizes
inv_N <- read.table(inv_sum, header = TRUE)

metal_tab <- left_join(eqtls, inv_N[, c("inversion", "N")], by = c("snps" = "inversion")) %>%
    mutate(RefAllele = "N", 
           AltAllele = "I",
           Marker = paste(snps, gene, sep = "-")
    )


write.table(metal_tab, file = paste0("results/cohorts_METAL_tabs/", prefix, "_metal.txt"),
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
