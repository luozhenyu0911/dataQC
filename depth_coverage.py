from __future__ import print_function
from __future__ import division
import sys
import pandas as pd
import numpy as np
import glob
import subprocess
input1=open(sys.argv[1],'r')
outf1=sys.argv[2]
for i in input1:
    i=i.strip()
    prefix=i.strip().split('.')[0]
    outfile_d=open(prefix+".depth",'w')
    outfile_c=open(prefix+".bedgraph",'w')
    subprocess.call(["samtools","depth","-a",i],stdout=outfile_d)
    subprocess.call(["bedtools","genomecov","-ibam",i,'-bga'],stdout=outfile_c)
    outfile_d.close()
    outfile_c.close()
    
def get_depth(x):
    depth=x[2].sum()/len(x)
    return depth
    
def coverage_1(x):
    value1=0
    value2=0
    for i in range(len(x)):
        if x.iloc[i,3]==0:
            value1+=int(x.iloc[i,2]-x.iloc[i,1])
        else:
            value2+=int(x.iloc[i,2]-x.iloc[i,1])
    return value2/(value1+value2)


depth_dict={}
for i in glob.glob('*.depth'):
    df=pd.read_csv(i.strip(),sep='\t',header=None)
    depth_dict[i.strip().split('.')[0]]=get_depth(df)
depth=pd.DataFrame.from_dict(depth_dict,orient='index',columns=['depth'])
depth=depth.reset_index().rename(columns={'index':'sample'})

coverage_dict={}
for i in glob.glob('*.bedgraph'):
    df2=pd.read_csv(i.strip(),sep='\t',header=None)
    coverage_dict[i.strip().split('.')[0]]=coverage_1(df2)
coverage=pd.DataFrame.from_dict(coverage_dict,orient='index',columns=['coverage'])
coverage=coverage.reset_index().rename(columns={'index':'sample'})

df1=pd.merge(depth, coverage, left_on="sample", right_on="sample", how="inner")
df1.to_csv(outf1,sep='\t',index=None)
