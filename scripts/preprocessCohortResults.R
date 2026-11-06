# Create summary statistics for cohort results
# 
# docker run --rm -it -v $PWD:$PWD -w $PWD godmc:1.2 R

## Load libraries
library(tidyverse)
library(cowplot)

## Load results
load("results/preprocess/all_masks_final.Rdata")

## Create list of bad probes
bad_probes <- filter(all_masks_final, final_mask)$Probe_ID

folders <- list.dirs("cohort_results", recursive = FALSE)

## Get summary statistics per cohort
cohort_sum <- lapply(folders, function(folder) {
    print(folder)
    results_Rdata <- list.files(folder, pattern = "inversionsSummary.txt", full.names = TRUE, recursive = TRUE)
    if(length(results_Rdata) == 1) {
        inv_sum <- read_delim(results_Rdata, delim = "\t", col_types = cols())   
    }
    inv_sum
})
names(cohort_sum) <- basename(folders)

cohort_sum_df <- bind_rows(cohort_sum, .id = "cohort")

## Remove inversions genotyped in few cohorts
tab_inv <- table(cohort_sum_df$inversion)
tab_inv
#  inv1_004  inv1_008 inv10_005 inv11_001 inv11_004 inv12_004 inv12_006 inv14_005 
#         7        37        38        35        10        38         7         7 
# inv16_009 inv17_007  inv2_002  inv2_013 inv21_005  inv3_003  inv6_002  inv6_006 
#        38        37        37        36        34         1         8        37 
#  inv7_003  inv7_005  inv7_011  inv8_001 
#        38        38        36        38 
## Remove inversions genotyped in few cohorts
good_invs <- names(tab_inv)[tab_inv > 30]


cohort_sum_filt <- filter(cohort_sum_df, inversion %in% good_invs)

png("figures/preproc/inversion_freq_distr.png", width = 800, height = 600)
ggplot(cohort_sum_filt, aes(x = freq)) +
    geom_histogram(binwidth = 1,  fill = "black") +
    labs(title = "Inversion Allele Frequency Distribution Across Cohorts",
        x = "Inversion Allele Frequency",
        y = "Count") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~inversion, ncol = 4)
dev.off()

png("figures/preproc/Nsnp_distr.png", width = 800, height = 600)
ggplot(cohort_sum_filt, aes(x = n_SNPs)) +
    geom_histogram(binwidth = 0.1,  fill = "black") +
    labs(title = "Number of SNPs Distribution Across Cohorts",
        x = "Number of SNPs",
        y = "Count") +
    theme_bw() +
    scale_x_log10() +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~inversion, ncol = 4)
dev.off()

png("figures/preproc/Nsnp_invFreq_distr.png", width = 800, height = 600)
ggplot(cohort_sum_filt, aes(y = n_SNPs, x = freq)) +
    geom_point() +
    labs(title = "Number of SNPs vs Inversion Allele Frequency",
        y = "Number of SNPs",
        x = "Inversion Allele Frequency (%)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~inversion, ncol = 4, scales = "free")
dev.off()

cohort_sum_filt_stats <- cohort_sum_filt %>%
    filter(inversion != "inv2_002") %>%
    group_by(inversion) %>%
    mutate(median = median(freq),
        outlier = ifelse(freq - median < -20 | freq - median > 20, TRUE, FALSE))

good_invs <- good_invs[!good_invs == "inv2_002"]


bad_cohort_inv <- cohort_sum_filt_stats %>%
    filter(outlier) %>%
    select(cohort, inversion)
save(good_invs, bad_cohort_inv, file = "results/preprocess/good_invs.Rdata")

## Load cohort results
cohort_res <- lapply(folders, function(folder) {
    print(folder)
    results_Rdata <- list.files(folder, pattern = "invmeqtl.Rdata", full.names = TRUE, recursive = TRUE)
    if(length(results_Rdata) == 1) {
        load(results_Rdata)
        eqtls <-  me$all$eqtls
    }
    eqtls
})
names(cohort_res) <- basename(folders)


cohort_res_filt <- lapply(cohort_res, function(eqtls) {
    filter(eqtls, !gene %in% bad_probes & snps %in% good_invs)
})
save(cohort_res_filt, file = "results/preprocess/cohort_results.Rdata")
sapply(cohort_res, function(x) length(unique(x$gene)))
sapply(cohort_res_filt, function(x) length(unique(x$gene)))

## Create METAL files per cohort
lapply(names(cohort_res_filt), function(cohort_name) {
    print(cohort_name)
    eqtls <- cohort_res_filt[[cohort_name]]

    if (cohort_name %in% bad_cohort_inv$cohort) {
        bad_invs_cohort <- bad_cohort_inv %>%
            filter(cohort == cohort_name) %>%
            pull(inversion)
        eqtls <- filter(eqtls, !snps %in% bad_invs_cohort)
    }

    inv_sum <- cohort_sum_filt %>%
        filter(cohort == cohort_name) %>%
        select(inversion, N)

    metal_tab <- left_join(eqtls, inv_sum[, c("inversion", "N")], by = c("snps" = "inversion")) %>%
        mutate(RefAllele = "N", 
               AltAllele = "I",
           Marker = paste(snps, gene, sep = "-")
    )
 write.table(metal_tab, file = paste0("results/cohorts_METAL_tabs_filtered/", cohort_name, "_metal.txt"),
   sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

})



## Load sample size
cohort_samplesize_df <- cohort_sum_filt %>%
    group_by(cohort) %>%
    summarise(n = max(N))


## Compute lambda for each cohort
getLambda <- function(df, inv){
    tab <- subset(df, snps == inv)
    pvals <- tab$pvalue
    observed <- -log10(sort(pvals))
    expected <- -log10(ppoints(length(pvals)))

    chi2_obs <- qchisq(1 - sort(pvals), df = 1)
    median_expected <- qchisq(0.5, df = 1)
    lambda <- median(chi2_obs) / median_expected
    lambda
}

lambdas <- lapply(cohort_res_filt, function(df) {
    invs <- unique(df$snps)
    lambda <- sapply(invs, function(inv) getLambda(df, inv))
    tibble(inversion = invs, lambda = lambda)
}) 
lambdas_df <- bind_rows(lambdas, .id = "cohort")

png("figures/preproc/lambda_per_cohort.png", width = 800, height = 600)
ggplot(lambdas_df, aes(x = cohort, y = lambda)) +
    geom_boxplot() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title = "Lambda per cohort",
        x = "Cohort",
        y = expression(lambda)) +   
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

png("figures/preproc/lambda_sd_cohort.png", width = 800, height = 600)
lambdas_df %>%
    group_by(cohort) %>%
    summarise(median_lambda = median(lambda),
    sd_lambda = sd(lambda)) %>%
    left_join(cohort_samplesize_df) %>%
    ggplot(aes(x = n, y = sd_lambda)) +
    geom_point() +
    theme_bw() 
dev.off()

##########################################################################################################################################
## A partir de aquí no
##########################################################################################################################################

plotQQ <- function(df, inv){
    tab <- subset(df, snps == inv)
    pvals <- tab$pvalue
    observed <- -log10(sort(pvals))
    expected <- -log10(ppoints(length(pvals)))

    chi2_obs <- qchisq(1 - sort(pvals), df = 1)
    median_expected <- qchisq(0.5, df = 1)
    lambda <- median(chi2_obs) / median_expected

    qq <- ggplot(tab, aes(sample = -log10(pvalue))) +
        stat_qq(color = "blue") +
        stat_qq_line(color = "red") +
        labs(title = inv,
            x = "Theoretical Quantiles",
            y = "Sample Quantiles") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
         annotate("text", 
           x = min(expected) + 0.1, 
           y = max(observed) - 0.1, 
           label = paste("λ =", round(lambda, 2)),
           hjust = 1, color = "black", size = 5)
    qq
}
f5l <- cohort_res_filt[["F5L"]]
f5l_qq <- lapply(good_invs, function(inv) {
    print(inv)

    plots <- plotQQ(f5l, inv)
})
png("figures/preproc/qq_F5L.png", width = 800, height = 800)
plot_grid(plotlist = f5l_qq, ncol = 4)
dev.off()

bib <- cohort_res_filt[["BIB"]]
bib_qq <- lapply(good_invs, function(inv) {
    print(inv)

    plots <- plotQQ(bib, inv)
})
png("figures/preproc/qq_BIB.png", width = 800, height = 800)
plot_grid(plotlist = bib_qq, ncol = 4)
dev.off()

png("figures/preproc/BIB_inv1_008_betas.png", width = 400, height = 400)
hist(subset(bib, snps == "inv1_008")$beta, breaks = 100, main = "BIB inv1_008 betas", 
     xlab = "Effect size (beta)")
dev.off()

png("figures/preproc/BIB_inv1_008_pvals.png", width = 400, height = 400)
hist(subset(bib, snps == "inv1_008")$pvalue, breaks = 100, main = "BIB inv1_008 p-values",
     xlab = "P-value")
dev.off()

png("figures/preproc/BIB_inv2_002_betas.png", width = 400, height = 400)
hist(subset(bib, snps == "inv2_002")$beta, breaks = 100, main = "BIB inv2_002 betas",
     xlab = "Effect size (beta)")
dev.off()

png("figures/preproc/BIB_inv2_002_pvals.png", width = 400, height = 400)
hist(subset(bib, snps == "inv2_002")$pvalue, breaks = 100, main = "BIB inv2_002 p-values",
     xlab = "P-value")
dev.off()

png("figures/preproc/BIB_inv17_007_pvals.png", width = 400, height = 400)
hist(subset(bib, snps == "inv17_007")$pvalue, breaks = 100, main = "BIB inv17_007 p-values",
     xlab = "P-value")
dev.off()