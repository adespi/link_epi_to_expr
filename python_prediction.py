#!/home/antoine/miniconda3/envs/kipoi-DeepSEA__predict/bin/python
#sys.argv=["program_name","chr_pos","nbr_intervals","gene_name",1,100,500000]
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
start = datetime.datetime.now()

old_stdout = sys.stdout
sys.stdout = open(os.devnull, "w")
model = kipoi.get_model('DeepSEA/predict')
sys.stdout = old_stdout
#sys.argv=[" ","5_158526788",5000,"gene_name",1,100,500000,100]
#sys.argv=[" ","9_37034476",5000,"PAX5",1,100,500000,100]
sys.argv[2]=int(sys.argv[2])
sys.argv[4]=int(sys.argv[4])
sys.argv[5]=int(sys.argv[5])
sys.argv[6]=int(sys.argv[6])
batch_size=int(sys.argv[7])
#print(sys.argv[1])
with open("temp/"+sys.argv[1]+"/expression/"+sys.argv[1]+".tsv") as f:
    content = f.readlines()

content = [x.strip() for x in content]
content[0] = content[0].split()
content[1] = content[1].split()
answ=np.empty([sys.argv[2],445,919])
expr=np.empty([445])
i=0
for g in sorted(glob.glob("temp/"+sys.argv[1]+"/fa_output/out"+sys.argv[1]+"_*.fa.gz")):
    dl_kwargs = {'intervals_file': 'temp/'+sys.argv[1]+'/intervals/'+os.path.splitext(os.path.basename(g))[0], 'fasta_file': g, "num_chr_fasta": "False"}
    warnings.filterwarnings('ignore',category=FutureWarning)
    dl = model.default_dataloader(**dl_kwargs)
    #it = dl.batch_iter(batch_size=sys.argv[2])
    it = dl.batch_iter(batch_size=batch_size)
    for b in range(0,sys.argv[2],batch_size):
        batch = next(it)
        warnings.filterwarnings('default',category=FutureWarning)
        if (sys.argv[2]<=b+batch_size):
            answ[b:sys.argv[2],i,:]=model.predict_on_batch(batch['inputs'])
        else:
            answ[b:b+batch_size,i,:]=model.predict_on_batch(batch['inputs'])
    expr[i]=content[1][content[0].index(re.split('_|\.',g)[-3])]
    i+=1

start_cor = datetime.datetime.now()
correlations=np.empty([919,sys.argv[2]])
np.seterr(invalid='ignore')
for i in range(sys.argv[2]):
    correlations[:,i]=np.corrcoef(np.transpose(answ[i,:,:]),expr)[-1][:-1]
np.seterr(invalid='warn')
#print(i)

column_names = np.arange(sys.argv[4],sys.argv[6],sys.argv[5])-1-sys.argv[6]/2
with open('deepsea_postprocessing/predictor.names') as f:
    row_names = f.read().splitlines()

df = pd.DataFrame(correlations, columns=column_names, index=row_names)
df.to_csv("correlations_para/correlations_"+sys.argv[1]+"_"+sys.argv[3]+".csv.gz")
#np.savez_compressed(/"+sys.argv[1]+"/fa_output/out"+sys.argv[1]+"_*.fa.gz',answ) #53 34
#np.savetxt("correlations_para/correlations_"+sys.argv[1]+"_"+sys.argv[3]+".csv", correlations, delimiter=",")

test_time = open("correlations_para/correlations_"+sys.argv[1]+"_"+sys.argv[3]+".time", 'w')
test_time.write(str(start)+"\n"+str(start_cor)+"\n"+str(datetime.datetime.now()))
test_time.close()
