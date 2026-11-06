

ewas_res <- read_table("data/ewascatalog-results.txt.gz")
eqtm_res <- read_delim("data/eQTM_autosome_adj.cells_SIG.txt.gz", delim = "\t")
trans_res <- read_delim("data/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz", delim = "\t")
## Inv8 trans effects
trans_8 <- subset(meta_cpg, inversion == "inv8_001" & chr != "chr8" & adj_pvalue < 0.05)

trans_8_phenos <- left_join(trans_8, ewas_res, by = "CpG")
trans_8_eqtm <- left_join(trans_8, eqtm_res, by = "CpG")

## https://pubmed.ncbi.nlm.nih.gov/32490597/ a CpG found there with autism

trans_17 <- subset(meta_cpg, inversion == "inv17_007" & chr != "chr17" & adj_pvalue < 0.05)

trans_17_phenos <- left_join(trans_17, ewas_res, by = "CpG")
trans_17_eqtm <- left_join(trans_17, eqtm_res, by = "CpG")

# MAPK8IP1 associated with autism
map <- subset(trans_res, GeneSymbol == "MAP3K13")