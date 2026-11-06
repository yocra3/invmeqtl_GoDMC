#!/bin/bash
#===============================================================================
#
#          FILE:  extract_cohort_results.sh
#
#   DESCRIPTION:  Extract the .tgz result file from each cohort. Then, converts the 
#                 results from Rdata format to a tab file and copy it to the results
#                 folder.
#
#
#===============================================================================

for folder in cohort_results/*/
do
    tgz=$folder/*.tgz
    if [ -e $tgz ] && [ ! -d $folder/results ];
        then
        prefix=$(basename "$folder")
        echo "Uncompressing and creating file for " $prefix
        tar -xvzf $tgz -C $folder
        Rscript bin/createMETALfile.R $folder/results/08/invmeqtl.Rdata $folder/results/08/inversionsSummary.txt $prefix
    fi
done

