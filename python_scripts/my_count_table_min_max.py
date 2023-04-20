#!/usr/bin/env python3
import pandas as pd
import subprocess
import sys
import warnings
import os
  
if os.path.exists("table.tsv"):
    uncompress_table='table.tsv'            
else:
    uncompress_table='results/qiime2/abundance_tables/feature-table.tsv'

# adapted from count_table_minmax_reads.py @author Daniel Straub
# collected from nf-core/ampliseq
# read tsv and skip first two rows
data = pd.read_csv(uncompress_table, sep="\t", skiprows=[0, 1], header=None)  # count table

# drop feature ids
df = data.drop(data.columns[0], axis=1)

# make sums
sums = df.sum()

# we want minimum values
mindepth = int(sums.min())
maxdepth = int(sums.max())

if mindepth > 10000:
    print("Use the sampling depth of " +str(mindepth)+" for rarefaction")
elif mindepth < 10000 and maxdepth > 5000: 
    print("WARNING The sampling depth of "+str(mindepth)+" is quite small for rarefaction")
elif mindepth < 5000 and mindepth > 1000:
    print("WARNING The sampling depth of "+str(mindepth)+" is very small for rarefaction")
elif mindepth < 1000:
    print("WARNING The sampling depth of "+str(mindepth)+" seems too small for rarefaction")
else:
    print("ERROR this shouldn't happen")
    exit(1)

        
#check values
        
if maxdepth > 75000:
    maxdepth = 75000
        
if maxdepth > 5000:
    maxsteps=250
else:
    maxsteps=(maxdepth/20)

file = open("srs_curve_val.txt", "w")
file.write(str(maxdepth))
file.close