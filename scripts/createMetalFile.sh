#!/bin/bash

cat <<EOF > scripts/runMetal_filtered.txt

VERBOSE OFF

# Describe and process input files
MARKER   Marker
WEIGHT   N
ALLELE   RefAllele AltAllele
EFFECT   beta
STDERR   SE
PVAL     pvalue
EOF

for file in results/cohorts_METAL_tabs_filtered/*.txt
do
    echo "PROCESS $file" >> scripts/runMetal_filtered.txt
done

cat <<EOF >> scripts/runMetal_filtered.txt

# Execute meta-analysis
OUTFILE results/meta_analysis/invqtl_metanalysis .tbl
ANALYZE
EOF