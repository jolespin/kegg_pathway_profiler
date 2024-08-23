#!/usr/bin/env bash
OUTPUT_DIRECTORY=${1:-"."}


# Pull http://rest.kegg.jp/list/module 
wget -O - http://rest.kegg.jp/list/module > ${OUTPUT_DIRECTORY}/pathway_names.tsv

# Pull definitions and classes
>${OUTPUT_DIRECTORY}/pathway_definitions.tsv
>${OUTPUT_DIRECTORY}/pathway_classes.tsv

for ID in $(cut -f1 ${OUTPUT_DIRECTORY}/pathway_names.tsv);
do
    wget -O - http://rest.kegg.jp/get/${ID} | while IFS= read -r line; do
        case "$line" in
            DEFINITION*)
                echo -e "${ID}\t${line:12}" >> ${OUTPUT_DIRECTORY}/pathway_definitions.tsv
                ;;
            CLASS*)
                echo -e "${ID}\t${line:12}" >> ${OUTPUT_DIRECTORY}/pathway_classes.tsv
                ;;
        esac
    done
done