# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Install R packages
RUN R -e 'BiocManager::install(c("tidyverse", "ggExtra", "cowplot", "ggmanh", "ggrepel", "ggbreak", "IlluminaHumanMethylationEPICv2anno.20a1.hg38", "IlluminaHumanMethylation450kanno.ilmn12.hg19"))'

# Install extra annotation
RUN R -e 'BiocManager::install(c("IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))'

# Install scoreInvHap  (v1.1)
RUN R -e 'BiocManager::install(c("scoreInvHap"))'

# Install Gviz for plotting (v1.2)
RUN R -e 'BiocManager::install(c("Gviz"))'

# Install enrichR for tests (v1.3)
RUN R -e 'BiocManager::install(c("enrichR"))'

# Install LOLA and OmicCircos (v1.4)
RUN R -e 'BiocManager::install(c("LOLA", "OmicCircos", "simpleCache"))'
