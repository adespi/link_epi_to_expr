#!/home/antoine/miniconda3/envs/kipoi-DeepSEA__predict/bin/python
#sys.argv=["program_name","chr_pos","nbr_intervals","gene_name",1,100,500000,"batch_size"]
#sys.argv=[" ","5_158526788",5000,"EBF1",1,100,500000,100]
#sys.argv=[" ","9_37034476",5000,"PAX5",1,100,500000,100]
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
import matplotlib.pyplot as plt
from sklearn.preprocessing import quantile_transform

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
with open("temp/"+sys.argv[1]+"/expression/"+sys.argv[1]+".tsv") as f:
    content = f.readlines()

content = [x.strip() for x in content]
content[0] = content[0].split()
content[1] = content[1].split()

#create empty files for predictions and expression, will be filled in loop
answ=np.empty([sys.argv[2],445,919])
expr=np.empty([445])
i=0

#predict genomic marks and load expression for each patient
for g in sorted(glob.glob("temp/"+sys.argv[1]+"/fa_output/out"+sys.argv[1]+"_*.fa.gz")):
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
    expr[i]=content[1][content[0].index(re.split('_|\.',g)[-3])]
    i+=1


start_cor = datetime.datetime.now()
#create empty correlation file
correlations=np.empty([919,sys.argv[2]])
#ignore warnings when std == 0
np.seterr(divide='ignore', invalid='ignore')
#compute all the correlations
expr=quantile_transform(expr.reshape(-1,1), output_distribution="normal", copy=True)[:,0]
for i in range(sys.argv[2]):
    correlations[:,i]=np.corrcoef(np.transpose(quantile_transform(answ[i,:,:], output_distribution="normal", copy=True)),expr)[-1][:-1]
np.seterr(divide='warn', invalid='warn')
#print(i)

#convert np_array to pd DataFrame to add column and row names and then save correlations to .csv.gz
column_names = np.arange(sys.argv[4],sys.argv[6],sys.argv[5])-1-sys.argv[6]/2
with open('deepsea_postprocessing/predictor.names') as f:
    row_names = f.read().splitlines()

x=np.log(np.max(np.max(answ, axis=1)/np.min(answ, axis=1),axis=1))
y=np.max(correlations,axis=0)
plt.plot(x,y,".")
plt.title("prediction magnitude vs correlation")
plt.xlabel("log(magnitude)")
plt.ylabel("correlation")

plt.savefig("correlations/prediction_magnitude_vs_correlation_{}_{}".format(sys.argv[1],sys.argv[3]))
plt.close()
plt.hist2d(x, y, (100, 100), cmap=plt.cm.jet)
plt.colorbar()
plt.title("prediction magnitude vs correlation heatmap")
plt.xlabel("log(magnitude)")
plt.ylabel("correlation")

plt.savefig("correlations/prediction_magnitude_vs_correlation_heatmap_{}_{}".format(sys.argv[1],sys.argv[3]))
plt.close()

nbr_points=np.zeros([2,2])
nbr_points[0,0]=sum(np.logical_and(y>0.18,x<1.5))
nbr_points[0,1]=sum(np.logical_and(y>0.18,x>1.5))
nbr_points[1,0]=sum(np.logical_and(y<0.18,x<1.5))
nbr_points[1,1]=sum(np.logical_and(y<0.18,x>1.5))

np.savetxt("correlations/prediction_range_vs_correlation_values_{}_{}".format(sys.argv[1],sys.argv[3]),nbr_points)


df = pd.DataFrame(correlations, columns=column_names, index=row_names)

position=df.abs().max().idxmax()
position_i=df.columns.get_loc(position)
mark=df[position].abs().idxmax()
mark_i=df.index.get_loc(mark)
X=np.copy(answ[position_i,:,mark_i])
Y=np.copy(expr)
X=X.reshape(-1, 1)
Y=Y.reshape(-1, 1)
X=quantile_transform(X, output_distribution="normal", copy=True)
Y=quantile_transform(Y, output_distribution="normal", copy=True)
for x in range(-5,6):
   if len(Y[np.round(X)==x])>0:
      a = plt.violinplot(Y[np.round(X)==x],[x])
      for pc in a['bodies']:
         pc.set_color('C0')
      a['cbars'].set_color('C0')
      a['cmins'].set_color('C0')
      a['cmaxes'].set_color('C0')
#plt.plot(X, Y,'.')
plt.xlabel(mark)
plt.ylabel(sys.argv[3]+" expression")
plt.title("Best correlation normalized for {} ({})".format(sys.argv[3],np.corrcoef(X.reshape(-1), Y.reshape(-1))[0,1]))
#plt.show()
plt.savefig("correlations/best_correlation_normalized_{}_{}.png".format(sys.argv[1],sys.argv[3]))





df.to_csv("correlations/correlations_"+sys.argv[1]+"_"+sys.argv[3]+".csv.gz")
#np.savez_compressed(/"+sys.argv[1]+"/fa_output/out"+sys.argv[1]+"_*.fa.gz',answ) #53 34
#np.savetxt("correlations/correlations_"+sys.argv[1]+"_"+sys.argv[3]+".csv", correlations, delimiter=",")

test_time = open("correlations/correlations_"+sys.argv[1]+"_"+sys.argv[3]+".time", 'w')
test_time.write(str(start)+"\n"+str(start_cor)+"\n"+str(datetime.datetime.now()))
test_time.close()
