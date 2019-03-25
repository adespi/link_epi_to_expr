#!/home/antoine/miniconda3/envs/kipoi-DeepSEA__predict/bin/python

# change above line to point to local 
# python executable

#import time

with open('list_gene_positions_chr1') as f:
    listg = f.read().splitlines()

with open('list_snp_positions_chr1') as f:
    lists = f.read().splitlines()

for g in listg:
    for s in lists:
        if int(s)>(int(g)-250000) and int(s)<(int(g)+250000):
            print (int(g))
