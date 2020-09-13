from os import listdir
from os.path import join
import pandas as pd
from operator import add
import time


mypath = "ARC/work-od-pair"
files = [f for f in listdir(mypath) if ".csv" in f]
df_all = None
for file in files[0:3]:
    t0 = time.time()
    df = pd.read_csv(join(mypath,file))
    df = df.groupby(["orig-block", "dest-block"]).apply(len)
    df = pd.DataFrame(df,columns=['Count'])
    print('test 1')
    print(time.time()-t0)
    if df_all is None:
        df_all = df
    else:
        df_all = df_all.combine(df, add, fill_value=0)
    print('test 2')
    print(time.time()-t0)
df_all.to_csv("arc_od_pairs.csv")