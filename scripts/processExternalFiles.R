# Process external files
# 
# docker run --rm -it -v $PWD:$PWD -w $PWD godmc:1.2 R

# Combine the masks from the different arrays

## Load libraries
library(tidyverse)
## Load masks
mask_epic2 <- read_delim("data/EPICv2.hg38.mask.tsv.gz", delim = "\t")
mask_epic_plus <- read_delim("data/EPIC+.hg38.mask.tsv.gz", delim = "\t")
mask_epic <- read_delim("data/EPIC.hg38.mask.tsv.gz", delim = "\t")
mask_450 <- read_delim("data/HM450.hg38.mask.tsv.gz", delim = "\t")

all_masks <- full_join(mask_epic2 %>% mutate(Probe_ID = gsub("_[A-Z][A-Z][0-9]*$", "", Probe_ID)),
    mask_epic_plus %>% mutate(Probe_ID = gsub("_.*$", "", Probe_ID)), 
    by = "Probe_ID", suffix = c(".epic2", ".epic_plus")) %>%
    full_join(mask_epic %>% mutate(Probe_ID = probeID), by = "Probe_ID") %>%
    full_join(mask_450 %>% mutate(Probe_ID = probeID), by = "Probe_ID", 
    suffix = c(".epic", ".450k")) %>%
    select(Probe_ID, starts_with("M_general"), starts_with("MASK_general")) %>%
    filter(grepl("cg", Probe_ID))

probes <- all_masks$Probe_ID
table(nchar(probes))

table(all_masks$M_general.epic2, all_masks$M_general.epic_plus)
table(all_masks$M_general.epic2, all_masks$MASK_general.epic)
table(all_masks$M_general.epic2, all_masks$MASK_general.450k)
table(all_masks$M_general.epic_plus, all_masks$MASK_general.epic)
table(all_masks$M_general.epic_plus, all_masks$MASK_general.450k)
table(all_masks$MASK_general.epic, all_masks$MASK_general.450k)

all_masks_final <- all_masks %>%
    mutate(final_mask = ifelse(!is.na(M_general.epic2), M_general.epic2,
        ifelse(!is.na(M_general.epic_plus), M_general.epic_plus,
        ifelse(!is.na(MASK_general.epic), MASK_general.epic,
        ifelse(!is.na(MASK_general.450k), MASK_general.450k, NA)))))
save(all_masks_final, file = "results/preprocess/all_masks_final.Rdata")

## Select SNPs in trans and catalog only present in scoreInvHap 
library(scoreInvHap)

trans_all <- read_delim("data/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz", delim = "\t")
cis_all <- read_delim("data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz", delim = "\t")

inv_snps <- lapply(SNPsR2, function(x) names(x)[x > 0.6]) 
inv_snps_df <- tibble(
    SNP = unlist(inv_snps),
    inversion = rep(names(inv_snps), times = lengths(inv_snps))
)

trans_inv <- inner_join(trans_all, inv_snps_df, by = "SNP")
save(trans_inv, file = "results/preprocess/trans_eQTL_filtered.Rdata")

cis_inv <- inner_join(cis_all, inv_snps_df, by = "SNP")
save(cis_inv, file = "results/preprocess/cis_eQTL_filtered.Rdata")

## Combine eQTMs
framingham_eqtms <- read_delim("data/Framingham Heart Study_data.csv")
helix_eqtms <- read_delim("data/HELIX_ Whole Blood Autosomal_data.csv")
gtex_blood <- read_delim("data/GTEx_ WholeBlood_data.csv")
gtex_lung <- read_delim("data/GTEx_ Lung_data.csv")
gtex_colon <- read_delim("data/GTEx_ ColonTransverse_data.csv")
gtex_muscle <- read_delim("data/GTEx_ MuscleSkeletal_data.csv")

gtex_eqtms <- bind_rows(gtex_blood, gtex_lung, gtex_colon, gtex_muscle) %>%
    filter(!is.na(meth.chr) | !is.na(gene.chr)) %>%
    mutate(meth.chr = paste0("chr", meth.chr))
combined_eqtms <- bind_rows(framingham_eqtms, helix_eqtms, gtex_eqtms)
save(combined_eqtms, file = "results/preprocess/combined_eQTMs.Rdata")