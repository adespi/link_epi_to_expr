#!/home/antoine/miniconda3/envs/kipoi-DeepSEA__predict/bin/python
import sys
import pandas as pd
df= pd.read_csv(sys.argv[1], header=0)
print(df.idxmin(axis=1)[0])
