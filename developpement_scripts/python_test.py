python
import glob
import re
import os
import sys
import kipoi
import numpy as np
import datetime
import warnings
start = datetime.datetime.now()

model = kipoi.get_model('DeepSEA/predict')
sys.argv=[" ","9_37034476",5000,"PAX5",1,100,500000]
batch_size=100
sys.argv[2]=int(sys.argv[2])
with open("temp/"+sys.argv[1]+"/expression/"+sys.argv[1]+".tsv") as f:
    content = f.readlines()

content = [x.strip() for x in content]
content[0] = content[0].split()
content[1] = content[1].split()
answ=np.empty([sys.argv[2],445,919])
expr=np.empty([445])
i=0

g=sorted(glob.glob("temp/"+sys.argv[1]+"/fa_output/out"+sys.argv[1]+"_*.fa.gz"))[0]
#g=sorted(glob.glob("temp/"+sys.argv[1]+"/fa_output/out"+sys.argv[1]+"_*.fa"))[0]
dl_kwargs = {'intervals_file': 'temp/'+sys.argv[1]+'/intervals/'+os.path.splitext(os.path.basename(g))[0], 'fasta_file': g, "num_chr_fasta": "False"}
#dl_kwargs = {'intervals_file': 'temp/'+sys.argv[1]+'/intervals/'+os.path.basename(g), 'fasta_file': g, "num_chr_fasta": "False"}

dl = model.default_dataloader(**dl_kwargs)

for b in range(0,sys.argv[2],batch_size):
    it = dl.batch_iter(batch_size=batch_size)
    batch = next(it)
    #warnings.filterwarnings('default',category=FutureWarning)
    if (sys.argv[2]<=b+batch_size):
        answ[b:sys.argv[2],i,:]=model.predict_on_batch(batch['inputs'])
    else:
        answ[b:b+batch_size,i,:]=model.predict_on_batch(batch['inputs'])

expr[i]=content[1][content[0].index(re.split('_|\.',g)[-3])]
correlations=np.empty([919,sys.argv[2]])
correlations[:,i]=np.corrcoef(np.transpose(answ[i,:,:]),expr)[-1][:-1]

import pandas as pd
column_names = np.arange(sys.argv[4],sys.argv[5],sys.argv[6])-1
with open('deepsea_postprocessing/predictor.names') as f:
    row_names = f.read().splitlines()

df = pd.DataFrame(my_data, columns=column_names, index=row_names)

###

b=range(0,sys.argv[2],batch_size)[0]
dl = model.default_dataloader(**dl_kwargs)
it = dl.batch_iter(batch_size=100)
batch = next(it)
answ[b:b+batch_size,i,:]=model.predict_on_batch(batch['inputs'])

###

it = dl.batch_iter(batch_size=sys.argv[2])
batch = next(it)
answ[:batch_size,i,:]=model.predict_on_batch(batch['inputs'])
expr[i]=content[1][content[0].index(re.split('_|\.',g)[-3])]
