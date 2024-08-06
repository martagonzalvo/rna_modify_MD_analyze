
# PLOT ANALYSIS 

# Usually have several replicates of same simulation with same name and different numbers, e.g. ogdup1, ogdup2, ... and long2dup1, long2dup2 etc. 
# This file plots all replicates of each type of simulation together

# python visualiz_analysis.py resultsfolder


import pandas as pd
import matplotlib.pyplot as plt
import sys, os, time

start = time.time()

folder = sys.argv[1]
cwd = os.getcwd()

path = cwd+'/'+folder

toanalyze = []
for file in os.listdir(path):
    if 'csv' in file:
        toanalyze.append(file)
toanalyze=sorted(toanalyze)


sample = toanalyze[0]
dfsample = pd.read_csv('{}/{}'.format(path,toanalyze[0]))
colns = dfsample.columns


typesfiles = [ 
    'md',
    ]

# plot each type of file separately
for coln in colns[2:]:
    for typefile in typesfiles:
        for f, file in enumerate(toanalyze):
            if typefile in file:
                print('Plotting',coln, file)
                df = pd.read_csv(path+'/'+file)

                prename = file.split('_')[1]
                name = prename.split('.')[0]
                plt.plot(df[coln], label=name)
            else:
                continue
        plt.title(coln)
        plt.xlabel('Simulation snapshots')
        plt.ylabel('Measure')
        plt.legend()
        plt.savefig("{}/{}{}.png".format(folder,coln,typefile), bbox_inches="tight")
        plt.clf()