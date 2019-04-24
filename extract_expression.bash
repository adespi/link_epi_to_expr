#!/bin/bash
if [ -e expression_temp.tsv ]; then
    rm expression_temp.tsv
fi
for list_of_genes in "some_genes";do #"some_genes_32_Tissue-specific_regulatory_networks_FANTOM5-v1" "all_genes"
    for TF in `tac "gene_info_"$list_of_genes".txt" |sed '$d'|cut -f1`;do
        read -r chromosome gene <<<$(less ../data/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |cut -f1-5|grep $TF|cut -f3-4)
        less ../data/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |grep -e $'\t'$chromosome$'\t'$gene$'\t'>>expression_temp.tsv
    done
    less ../data/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz|grep -e '^T'|cat - expression_temp.tsv > "expression_"$list_of_genes".tsv"
    rm expression_temp.tsv
done
