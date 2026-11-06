# Explore results of the meta-analysis
# 
# docker run -it -v $PWD:$PWD -w $PWD godmc:1.2 R

# Load libraries
library(tidyverse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(rtracklayer)
library(GenomicRanges)


# Load annotation 
anno_450k <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
cg_450k <- rownames(anno_450k)

anno_v1 <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
cg_v1 <- rownames(anno_v1)

anno_v2 <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations
cg_38 <- gsub("_.*$", "", rownames(anno_v2))

## Create list of bad probes
load("results/preprocess/good_invs.Rdata")
load("results/preprocess/all_masks_final.Rdata")

bad_probes <- filter(all_masks_final, final_mask)$Probe_ID


# Load meta-analysis results
meta <- read_table("results/meta_analysis/invqtl_metanalysis1.tbl",
                  col_types = cols(
                     .default = col_guess(), 
                     `P-value` = col_character()
                   )) %>%
    mutate(inversion = gsub("-.*", "", MarkerName),
           CpG = gsub(".*-", "", MarkerName)) %>%
    filter(inversion %in% good_invs & !CpG %in% bad_probes)

meta <- mutate(meta, 
    P.Value = as.numeric(`P-value`),
    P.Value = ifelse(P.Value == 0, 1.482197e-323, P.Value),
    adj_pvalue = p.adjust(`P-value`, method = "bonferroni"))
all_cpgs <- unique(meta$CpG)

## Check CpGs are present in all annotation files
mean(all_cpgs %in% c(cg_38, cg_450k, cg_v1))

## Create Annotation for CpGs
anno_v2$CpG <- cg_38
anno <- subset(anno_v2, CpG %in% all_cpgs)

anno_v1$CpG <- rownames(anno_v1)
anno_450k$CpG <- rownames(anno_450k)

## Convert hg19 coordinates to hg38
anno_v1_sel <- subset(anno_v1, CpG %in% all_cpgs & !CpG %in% anno$CpG )
anno_450k_sel <- subset(anno_450k, CpG %in% all_cpgs & !CpG %in% anno$CpG & !CpG %in% anno_v1_sel$CpG )
anno_19 <- rbind(anno_v1_sel, anno_450k_sel)

chain_file <- "data/hg19ToHg38.over.chain"
chain <- import.chain(chain_file)

hg19_coords <- GRanges(seqnames = anno_19$chr,
                       ranges = IRanges(start = anno_19$pos, end = anno_19$pos))
hg38_coords <- liftOver(hg19_coords, chain) 

anno_19_lift <- anno_19 %>% 
    as_tibble() %>%
    mutate(chr_19 = chr,
           pos_19 = pos,
           chr = as.character(seqnames(hg38_coords)),
           pos = as.numeric(start(hg38_coords))) 

anno_38 <- rbind(select(as_tibble(anno), chr, pos, CpG), 
    select(anno_19_lift, chr, pos, CpG)) 

anno_38 <- anno_38 %>% 
    filter(!is.na(chr)) %>%
    distinct()
## Check CpGs are present in final annotation files
mean(all_cpgs %in% anno_38$CpG)

meta_cpg <- inner_join(meta, anno_38, by = "CpG", copy = TRUE)
save(meta_cpg, file = "results/meta_analysis/invqtl_results.Rdata")

