#!/home/antoine/miniconda3/envs/kipoi-DeepSEA__predict/bin/python
#sys.argv=["program_name","chr_pos","nbr_intervals","gene_name",1,100,500000,"batch_size"]
#sys.argv=[" ","5_158526788",5000,"EBF1",1,100,500000,20]
#sys.argv=[" ","9_37034476",5000,"PAX5",1,100,500000,20]
#sys.argv=[" ","2_127864931",5000,"BIN1",1,100,500000,20]
#python
import glob
import re
import os
import sys
import kipoi
import numpy as np
import datetime
import warnings
import pandas as pd

#go to folder where script is. Enables lauch of script from any folder
os.chdir(os.path.dirname(os.path.realpath(__file__)))

start = datetime.datetime.now()

#import DeepSEA model without printing stder
old_stdout = sys.stdout
sys.stdout = open(os.devnull, "w")
model = kipoi.get_model('DeepSEA/predict')
sys.stdout = old_stdout

#format input args strings->ints
sys.argv[2]=int(sys.argv[2])
sys.argv[4]=int(sys.argv[4])
sys.argv[5]=int(sys.argv[5])
sys.argv[6]=int(sys.argv[6])
batch_size=int(sys.argv[7])

#fetch expression data
dfexpr = pd.read_csv("expression_for_some_genes.tsv", sep="\t",header = 0)

#create empty files for predictions and expression, will be filled in loop
answ=np.empty([sys.argv[2],445,919])
expr=np.empty([dfexpr.shape[0],445])
i=0

#predict genomic marks and load expression for each patient
for g in sorted(glob.glob("temp/"+sys.argv[1]+"/fa_output/out"+sys.argv[1]+"_*.fa.gz")):
    print (i)
    dl_kwargs = {'intervals_file': 'temp/'+sys.argv[1]+'/intervals/'+os.path.splitext(os.path.basename(g))[0], 'fasta_file': g, "num_chr_fasta": "False"}
    #filter warnings to clear output
    warnings.filterwarnings('ignore',category=FutureWarning)
    #initiate the list of intervals
    dl = model.default_dataloader(**dl_kwargs)
    #set batch size for the predictions
    it = dl.batch_iter(batch_size=batch_size)
    #run all the predictions on all the batches
    for b in range(0,sys.argv[2],batch_size):
        batch = next(it)
        warnings.filterwarnings('default',category=FutureWarning)
        if (sys.argv[2]<=b+batch_size):
            answ[b:sys.argv[2],i,:]=model.predict_on_batch(batch['inputs'])
        else:
            answ[b:b+batch_size,i,:]=model.predict_on_batch(batch['inputs'])
    #get expression for the patient
    expr[:,i]=dfexpr[(re.split('_|\.',g)[-3])]
    i+=1


start_cor = datetime.datetime.now()

#informations for output file
gene_info = pd.read_csv("gene_info_some_genes.txt", sep="\t",header = 0,index_col=0)
column_names = list(pd.read_csv("correlations_small/correlations_"+sys.argv[1]+"_"+sys.argv[3]+".csv.gz",header = 0))
with open('deepsea_postprocessing/predictor.names') as f:
    row_names = f.read().splitlines()

gene_n=0
for gene in dfexpr["TargetID"]:
    gene=gene.split(".")[0]
    #create empty correlation file
    correlations=np.empty([919,sys.argv[2]])
    #ignore warnings when std == 0
    np.seterr(divide='ignore', invalid='ignore')
    #compute all the correlations
    for i in range(sys.argv[2]):
        correlations[:,i]=np.corrcoef(np.transpose(answ[i,:,:]),expr[gene_n])[-1][:-1]
    np.seterr(divide='warn', invalid='warn')
    #print(i)
    #convert np_array to pd DataFrame to add column and row names and then save correlations to .csv.gz
    df = pd.DataFrame(correlations, columns=column_names, index=row_names)
    df.to_csv("correlations/correlations_"+sys.argv[1]+"_"+sys.argv[3]+"/"+gene_info.loc[gene][4]+".csv.gz")
    gene_n+=1

test_time = open("correlations/correlations_"+sys.argv[1]+"_"+sys.argv[3]+"/time", 'w')
test_time.write(str(start)+"\n"+str(start_cor)+"\n"+str(datetime.datetime.now()))
test_time.close()
