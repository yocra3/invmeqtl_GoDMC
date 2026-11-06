# Explore results of the meta-analysis for all inversions
# 
# docker run -it -v $PWD:$PWD -w $PWD godmc:1.4 R

## Load libraries
library(tidyverse)
library(cowplot)
library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(ggmanh)
library(ggrepel)
library(ggbreak)
library(enrichR)
library(LOLA)
library(OmicCircos)

load("results/meta_analysis/invqtl_results.Rdata")
load("results/preprocess/cis_eQTL_filtered.Rdata")
load("results/preprocess/trans_eQTL_filtered.Rdata")
load("results/preprocess/combined_eQTMs.Rdata")

chain_file <- "data/hg19ToHg38.over.chain"
chain <- import.chain(chain_file)
# Distribution of sample sizes
png("figures/preproc/sample_size_distribution.png", res = 300, width = 1600, height = 1200)
ggplot(meta_cpg, aes(x = Weight)) +
    geom_histogram(binwidth = 500, fill ="black") +
    labs(title = "Meta-analysis sample sizes across meQTLs",
         x = "Sample size",
         y = "Count") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
dev.off()

## Filter CpGs present in at least 50% of samples
meta_cpg_filt <- subset(meta_cpg, Weight > 20e3)

## Make qqplots
plotQQ <- function(inv, df){
    tab <- subset(df, inversion == inv)
    pvals <- tab$P.Value
    observed <- -log10(sort(pvals))
    expected <- -log10(ppoints(length(pvals)))

    chi2_obs <- qchisq(1 - sort(pvals), df = 1)
    median_expected <- qchisq(0.5, df = 1)
    lambda <- median(chi2_obs) / median_expected

    qq <- ggplot(tab, aes(sample = -log10(P.Value))) +
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
inversions <- unique(meta_cpg$inversion)
for (inv in inversions){
    print(inv)
    png(paste0("figures/qq_plots/qq_", inv, ".png"))
    plotQQ(inv, meta_cpg) %>% print()
    dev.off()
    png(paste0("figures/qq_plots/qq_", inv, "filtered.png"))
    plotQQ(inv, meta_cpg_filt) %>% print()
    dev.off()
}
## Make Manhattan plot per inversion
n_tests <- nrow(meta_cpg)

inv_ranges37 <- read_table("data/inversion_ranges.txt", 
    col_names = c("chr", "start", "end", "inversion")) %>%
    mutate(chr = paste0("chr", chr))
inv_ranges37GR <- makeGRangesFromDataFrame(inv_ranges37)
names(inv_ranges37GR) <- inv_ranges37$inversion

inv_ranges38GR <- liftOver(inv_ranges37GR, chain) 
inv_ranges38GR <- lapply(inv_ranges38GR, range) %>% GRangesList() %>% unlist()

plotManhattan <- function(inv, df){

    tab <- subset(df, inversion == inv)
    tab$chr <- factor(tab$chr, levels = paste0("chr", c(1:22, "X")))
    tab$chrNum <- factor(gsub("chr", "", tab$chr), levels = c(1:22, "X"))
    tab$sel_CpG <- ifelse(tab$adj_pvalue < 0.05, tab$CpG, "")  
    tabGR <- makeGRangesFromDataFrame(tab, seqnames.field = "chr", 
        start.field = "pos", end.field = "pos", keep.extra.columns = TRUE)
    sel_CpGs <- subsetByOverlaps(tabGR, inv_ranges38GR[inv] + 100e3)
    highlight_colormap <- c("even" = "lightblue", "odd" = "lightgreen", inversion = "black")
    tab$color <- ifelse(tab$CpG %in% sel_CpGs$CpG, "inversion", 
        ifelse(as.numeric(tab$chrNum) %% 2 == 0, "even", "odd"))


    man_plot <- manhattan_plot(x = tab, pval.colname = "P.Value", 
                             signif = c(0.05/n_tests),
                             chr.colname = "chrNum", pos.colname = "pos",
                             rescale = TRUE,
                             plot.title = inv,
                             highlight.colname = "color",
                            highlight.col = highlight_colormap, 
                            color.by.highlight = TRUE) +
    theme(plot.title = element_text(hjust = 0.5)) 

    man_plot
}

inversions <- unique(meta_cpg$inversion)
for (inv in inversions){
    print(inv)
    png(paste0("figures/manhattans/manhattan_", inv, ".png"))
    plotManhattan(inv, meta_cpg_filt) %>% print()
    dev.off()
}

inv_order <- c("inv1_008", "inv2_013", "inv6_006", "inv7_003", "inv7_005", "inv7_011",
            "inv8_001", "inv10_005", "inv11_001", "inv12_004", "inv16_009", "inv17_007", "inv21_005")
manhattans <- lapply(inversions, function(inv){
   plotManhattan(inv, meta_cpg_filt) %>% print()
})
names(manhattans) <- inversions
png("figures/manhattans/manhattan_panel.png", res = 300, width = 3600, height = 2400)
plot_grid(plotlist = manhattans[inv_order], ncol = 4)
dev.off()

## Summary of inversion associations
inversionOverlap <- function(chr_inv, start_inv, end_inv, chr_cpg, cpg_pos,  window = 500e3){
    ifelse(chr_inv == chr_cpg & 
           (cpg_pos >= (start_inv - window) & cpg_pos <= (end_inv + window)),
           TRUE, FALSE)
}

inv_ranges_38_df <- as_tibble(inv_ranges38GR) %>%
    mutate(inversion = names(inv_ranges38GR),
           inv_chr = seqnames) %>%
    select(inversion, inv_chr, start, end)


getDistance <- function(pos, start, end){
    d_start <- abs(pos - start)
    d_end <- abs(pos - end)
    pmin(d_start, d_end)
}

meta_cpg_inv_annot <- meta_cpg_filt %>%
    left_join(inv_ranges_38_df, by = c("inversion" = "inversion")) %>%
    mutate(rel = ifelse(inversionOverlap(inv_chr, start, end, chr, pos, window = 0), "inv_region",
        ifelse(inversionOverlap(inv_chr, start, end, chr, pos, window = 1.5e6), "cis_extended",
            ifelse(chr == inv_chr, "trans_intra", "trans"))),
        distance = ifelse(rel %in% c("trans_intra", "cis"), getDistance(pos, start, end), NA))



cis_invmeqtl <- meta_cpg_inv_annot %>%
    mutate(inversion = factor(inversion, levels = inv_order),
            rel = factor(rel, levels = c("inv_region", "cis_extended", "trans_intra", "trans"))) %>%
    filter(adj_pvalue < 0.05) %>%
    group_by(inversion, rel) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    mutate(group = ifelse(rel %in% c("inv_region", "cis_extended"), "Cis", "Trans"))

table(meta_cpg_inv_annot$adj_pvalue < 0.05, meta_cpg_inv_annot$rel, meta_cpg_inv_annot$inversion)

sig_invmeqtl <- meta_cpg_inv_annot %>%
    filter(adj_pvalue < 0.05)



cis_eqtl <- cis_inv %>%
    filter(FDR < 0.05 & inversion %in% sig_invmeqtl$inversion) %>%
    group_by(inversion) %>%
    summarise(cis_eqtl = n_distinct(GeneSymbol)) %>%
    ungroup() %>%
    pivot_longer(cols = cis_eqtl, names_to = "rel", values_to = "n") %>%
    mutate(group = "Cis")


cis_eqtm_eqtl <- inner_join(select(sig_invmeqtl, inversion, CpG, Zscore, P.Value, adj_pvalue), 
    select(combined_eqtms, probeid, gene, beta, pvalue, FDR, tissue, cohort), by = c("CpG" = "probeid")) %>%
    inner_join(select(cis_inv,  GeneSymbol, Zscore, Pvalue, FDR, inversion), by = c("inversion" = "inversion", "gene" = "GeneSymbol")) %>%
    distinct(inversion, CpG, gene, .keep_all = TRUE) 

cis_eqtm_eqtl_sum  <- cis_eqtm_eqtl %>%
    select(inversion, gene) %>%
    distinct() %>%
    group_by(inversion) %>%
    summarise(eqtm_eqtl = n()) %>%
    ungroup() %>%
    pivot_longer(cols = eqtm_eqtl, names_to = "rel", values_to = "n") %>%
    mutate(group = "Cis")
    
trans_eqtl <- trans_inv %>%
    filter(FDR < 0.05) %>%
    group_by(inversion) %>%
    summarise(trans_eqtl = n_distinct(GeneSymbol)) %>%
    ungroup() %>%
    pivot_longer(cols = trans_eqtl, names_to = "rel", values_to = "n") %>%
    mutate(group = "Trans")

trans_eqtm_eqtl <- inner_join(select(sig_invmeqtl, inversion, CpG, P.Value, adj_pvalue), 
    select(combined_eqtms, probeid, gene, pvalue, FDR, tissue, cohort), by = c("CpG" = "probeid")) %>%
    inner_join(select(trans_inv, GeneSymbol, Pvalue, FDR, inversion), by = c("inversion" = "inversion", "gene" = "GeneSymbol")) %>%
    distinct(inversion, CpG, gene, .keep_all = TRUE)

trans_eqtm_eqtl_sum  <- trans_eqtm_eqtl %>%
    select(inversion, gene) %>%
    distinct() %>%
    group_by(inversion) %>%
    summarise(trans_eqtm_eqtl = n()) %>%
    ungroup() %>%
    pivot_longer(cols = trans_eqtm_eqtl, names_to = "rel", values_to = "n") %>%
    mutate(group = "Trans")

png("figures/general/invmeqtl_summary.png", res = 300, width = 2400, height = 1600)
bind_rows(cis_invmeqtl, cis_eqtl, cis_eqtm_eqtl_sum, trans_eqtl, trans_eqtm_eqtl_sum) %>%
    mutate(color = factor(rel, levels = c("inv_region", "cis_extended", "cis_eqtl", "eqtm_eqtl", "trans_intra", "trans", "trans_eqtl", "trans_eqtm_eqtl")),
    loc = recode(rel, inv_region = "Inv. Reg.", cis_extended = "Inv. Ext.", cis_eqtl = "inv-eQTL", eqtm_eqtl = "inv-eQTM-eQTL", 
    trans_intra = "Intrachrom.", trans = "Other Chrom.", trans_eqtl = "inv-eQTL", trans_eqtm_eqtl = "inv-eQTM-eQTL"),
    loc = factor(loc, levels = c("Inv. Reg.", "Inv. Ext.", "Intrachrom.", "Other Chrom.", "inv-eQTL", "inv-eQTM-eQTL")),
    inversion = factor(inversion, levels = rev(inv_order))) %>%
    ggplot(aes(x = loc, y = inversion, size = n, color = color)) +
    geom_point() +
    scale_color_manual(values = c("#CFE2F2", "#A6CEE3", "#1F78B4", "#08306B", "#FDD0A2", "#FDBF6F", "#FF7F00", "#7F2704")) +
    labs(title = "inversion meQTLs",
         x = "",
         y = "Inversion",
         color = "Location") +
    theme_bw() +
    guides(color = "none") +
    facet_grid(~group, scales = "free_x") +
    theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## Plot genomic context of inversion meQTLs
inversions <- unique(meta_cpg_inv_annot$inversion)
for(inv in inversions){

    print(inv)

    inv_range <- inv_ranges38GR[inv] 
    from <- start(inv_range) -500e3
    to <- end(inv_range) + 500e3
    chr <- as.character(seqnames(inv_range))

    cis_eqtls <- unique(subset(cis_inv, inversion == inv & FDR < 0.05)$GeneSymbol)
    axTrack <- GenomeAxisTrack()
    idxTrack <- IdeogramTrack(genome="hg38", chromosome=chr)

    knownGenes <- UcscTrack(
        genome = "hg38", chromosome = chr,
        track = "NCBI RefSeq", table = "ncbiRefSeqCurated", 
        from = from, to = to, trackType = "GeneRegionTrack", 
        rstarts = "exonStarts", rends = "exonEnds", 
        gene = "name2", symbol = "name2", transcript = "name",
        strand = "strand",
        collapseTranscripts = "meta", 
        transcriptAnnotation = "symbol",
        name = "RefSeq Genes")

    from <- min(start(ranges(knownGenes)))
    to <- max(end(ranges(knownGenes)))

    if(cis_eqtls %>% length() > 0){
        cis_eqtl_genesGR <- ranges(knownGenes)[ranges(knownGenes)$symbol %in% cis_eqtls]
        ranges(knownGenes) <- ranges(knownGenes)[!ranges(knownGenes)$symbol %in% cis_eqtls]
    } else{
        cis_eqtl_genesGR <- GRanges()
    }

    cis_eqtl_genes <- GeneRegionTrack(cis_eqtl_genesGR,
        genome = "hg38", chromosome = chr,
        name = "Cis eQTL genes",
        fill = "red",
        col = "red",
        transcriptAnnotation = "symbol",
        collapseTranscripts = "meta"
    )


    seqDups <- tryCatch(
        UcscTrack(genome = "hg38", chromosome = chr, 
                  track = "Segmental Dups", table = "genomicSuperDups", 
                  from = from, to = to,
                  trackType = "AnnotationTrack", 
                  start = "chromStart", end = "chromEnd",
                  name = "Seg. Dups", 
                  stacking = "dense"),
        error = function(e) {
            AnnotationTrack(name = "Seg. Dups", stacking = "dense", 
                            genome = "hg38", chromosome = chr)
        }
    )

    filt <- meta_cpg_inv_annot %>%
        filter(inversion ==  inv & pos > from & pos < to & chr == chr) %>%
        arrange(pos)

    dtrack <- DataTrack(data = filt$Zscore, start = filt$pos,
                    end = filt$pos, chromosome = chr, genome = "hg38", 
                    name = "Zscores", 
                    type = "p", 
                    legend = TRUE)
    displayPars(cis_eqtl_genes) <- list(rotation.title = 0)    
                
    ht <- HighlightTrack(trackList = list(knownGenes, cis_eqtl_genes, dtrack, seqDups),
                     start = start(inv_range), end = end(inv_range))

    png(paste0("figures/inversionContext/inversion_context_", inv, ".png"), res = 300, width = 1600, height = 1200)
    plotTracks(list(idxTrack, axTrack, ht), margin = 50)
    dev.off()


}

## LOLA enrichment for invmeQTL with trans effects
sel_invs <- c("inv8_001", "inv16_009", "inv17_007")
names(sel_invs) <- sel_invs

# Load LOLA database
regionDB <- loadRegionDB("data/LOLA_cache/t1/resources/regions/LOLACore/hg38")

## Define universe and query sets
### Add 100 bp from each position
universe <- makeGRangesFromDataFrame(meta_cpg_inv_annot %>%
    select(CpG, chr, pos) %>%
    distinct() %>%
    mutate(start = pos, end = pos), keep.extra.columns = TRUE) + 100
red_universe <- reduce(universe)

inv_list <- lapply(sel_invs, function(inv){
    gr <- makeGRangesFromDataFrame(sig_invmeqtl %>%
        filter(inversion == inv & rel %in% c("trans", "trans_intra")) %>%
        select(CpG, chr, pos) %>%
        mutate(start = pos, end = pos), keep.extra.columns = TRUE) + 100 
    reduce(gr)
}) %>% GRangesList()



lola_enrichments <- runLOLA(userSets = inv_list, 
            userUniverse = red_universe,
            regionDB = regionDB,
            cores = 10)

lola_sig_enrichments <- subset(lola_enrichments, pValueLog >= 7 & collection == "encode_tfbs") %>% select(-filename)
table(lola_sig_enrichments$userSet, lola_sig_enrichments$antibody)

inv17_lola <- subset(lola_enrichments, userSet == "inv17_007" & pValueLog >= 7) %>% select(-filename)
## Enriched in regions associated to promoters

inv8_lola <- subset(lola_enrichments, userSet == "inv8_001" & pValueLog >= 7) %>% select(-filename)
## Enriched in segmental dups

## Compute enrichments for inv8, inv17 and inv16
dbs <- listEnrichrDbs()

db.vec <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "DisGeNET", "Human_Phenotype_Ontology",
"OMIM_Expanded", "Jensen_DISEASES_Curated_2025", "TRRUST_Transcription_Factors_2019", "GO_Biological_Process_2025", "GWAS_Catalog_2023" )



enrich_tabs <- lapply(sel_invs, function(inv){
    
    trans_genes <- subset(trans_eqtm_eqtl, inversion == inv)$gene %>% unique()
    cis_genes <- subset(cis_eqtm_eqtl, inversion == inv)$gene %>% unique()
    all_genes <- c(cis_genes, trans_genes)

    trans_enrich <- enrichr(trans_genes, db.vec,  include_overlap = TRUE)
    cis_enrich <- enrichr(cis_genes, db.vec, include_overlap = TRUE)
    all_enrich <- enrichr(all_genes, db.vec, include_overlap = TRUE)

    names(trans_enrich) <- db.vec
    names(all_enrich) <- db.vec

    enrich_trans_df <-  Reduce(rbind,trans_enrich) %>%
        mutate(DB = rep(names(trans_enrich), sapply(trans_enrich, nrow)),
                n_genes = sapply(strsplit(Overlap, "/"), `[`, 1)) %>%
        filter(P.value < 0.05) 
    enrich_cis_df <-  Reduce(rbind,cis_enrich) %>%
        mutate(DB = rep(names(cis_enrich), sapply(cis_enrich, nrow)),
               n_genes = sapply(strsplit(Overlap, "/"), `[`, 1)) %>%
        filter(P.value < 0.05)
    enrich_all_df <-  Reduce(rbind,all_enrich) %>%
        mutate(DB = rep(names(all_enrich), sapply(all_enrich, nrow)),
               n_genes = sapply(strsplit(Overlap, "/"), `[`, 1)) %>%
        filter(P.value < 0.05) 
    list(trans = enrich_trans_df, cis = enrich_cis_df, all = enrich_all_df)
})
names(enrich_tabs) <- sel_invs

## Inv16 exploration
subset(cis_inv, inversion == "inv16_009" & GeneSymbol == "IL27")

subset(trans_inv, inversion == "inv16_009" & GeneSymbol %in% c("IRF1", "STAT1")) %>% data.frame()
inv16_snps <- subset(trans_inv, inversion == "inv16_009" & GeneSymbol %in% c("IRF1", "STAT1"))$SNP %>% unique()
Refs$inv16_009[inv16_snps]

subset(enrich_tabs$inv16_009$trans, n_genes > 2 & grepl("IRF1", Term))

## Inv17 exploration
subset(cis_inv, inversion == "inv17_007" & GeneSymbol == "HEXIM2")

subset(enrich_tabs$inv17_007$trans, n_genes > 2 & Adjusted.P.value < 0.05)
## No clear enrichment
subset(enrich_tabs$inv17_007$cis, n_genes > 2 & Adjusted.P.value < 0.05)

## Inv8 exploration
subset(enrich_tabs$inv8_001$trans, n_genes > 2 & Adjusted.P.value < 0.05)
subset(enrich_tabs$inv8_001$cis, n_genes > 2 & Adjusted.P.value < 0.05)
subset(enrich_tabs$inv8_001$cis, n_genes > 2 & DB != "GWAS_Catalog_2023")