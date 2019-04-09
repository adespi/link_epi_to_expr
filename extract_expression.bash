#!/bin/bash
#gene_info_some_genes_32_Tissue-specific_regulatory_networks_FANTOM5-v1.txt
for list_of_genes in "gene_info_some_genes.txt" ;do #"gene_info_some_genes_32_Tissue-specific_regulatory_networks_FANTOM5-v1.txt";do
    for TF in `tac $list_of_genes |sed '$d'|cut -f1`;do
        read -r chromosome gene <<<$(less ../data/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |cut -f1-5|grep $TF|cut -f3-4)
        less ../data/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz |grep -e $'\t'$chromosome$'\t'$gene$'\t'>>expression_for_some_genes_temp.tsv
    done
done
less ../data/geuvadis_expression/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz|grep -e '^T'|cat - expression_for_some_genes_temp.tsv >expression_for_some_genes.tsv
rm expression_for_some_genes_temp.tsv