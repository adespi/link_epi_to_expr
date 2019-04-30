#!/home/antoine/miniconda3/envs/kipoi-DeepSEA__predict/bin/python

import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#gene_info = pd.read_csv("gene_info_some_genes.txt", sep='\t',header=0)
gene_info = pd.read_csv("gene_info_some_genes_32_Tissue-specific_regulatory_networks_FANTOM5-v1.txt", sep='\t',header=0)
len(gene_info["Gene name"])
 
 
 
nbr_genes=len(gene_info["Gene name"])
max_corr=pd.DataFrame(np.zeros([nbr_genes,nbr_genes]), index=gene_info["Gene name"], columns=gene_info["Gene name"])
nbr_signif_corr=pd.DataFrame(np.zeros([nbr_genes,nbr_genes]), index=gene_info["Gene name"], columns=gene_info["Gene name"])
sign_cutoff=0.1550258
gene_index=0
gene_info["Gene name"]
for gene in gene_info["Gene name"]:
    print(gene_index)
    gene_folder=glob.glob("correlations/correlations_*_*_"+gene+"/")
    #print(gene_folder)
    if len(gene_folder) ==1 :
        sub_gene_index=0
        for sub_gene in gene_info["Gene name"]:
            sub_gene_file=glob.glob(gene_folder[0]+sub_gene+".csv.gz")
            #print(sub_gene_file)
            if len(sub_gene_file) ==1:
                if len(glob.glob(gene_folder[0]+sub_gene+".csv.gz.info"))==0 :
                    my_data = pd.read_csv(sub_gene_file[0], sep=',',header=0,index_col=0)
                    max_corr.values[gene_index,sub_gene_index]=my_data.abs().max().max()
                    nbr_signif_corr.values[gene_index,sub_gene_index]=np.sum(((my_data.abs())>sign_cutoff).values)
                    with open(sub_gene_file[0]+".info","w+") as f:
                        f.write(str(max_corr.values[gene_index,sub_gene_index])+"\n"+str(nbr_signif_corr.values[gene_index,sub_gene_index]))
                        f.close()
                else:
                    with open(sub_gene_file[0]+".info","r") as f:
                        max_corr.values[gene_index,sub_gene_index] = float(f.readline().strip())
                        nbr_signif_corr.values[gene_index,sub_gene_index] = (f.readline().strip())
            sub_gene_index+=1
    gene_index+=1