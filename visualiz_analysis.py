
# PLOT ANALYSIS 

# Usually have several replicates of same simulation with same name and different numbers, e.g. ogdup1, ogdup2, ... and long2dup1, long2dup2, ... etc. 
# This file plots all replicates of each type of simulation together

# python visualiz_analysis.py resultsfolder


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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


# putting all results in joint DataFrame for easier later analysis
for f, file in enumerate(toanalyze):
    df = pd.read_csv(path+'/'+file)

    if f == 0:
        avgdf = df
    else:
        avgdf = pd.concat([df,avgdf], axis=0)
    print('adding', file)

avgdf.to_csv('{}/alldatadf.csv'.format(path))


# Plotting averages of all columns by type

for coln in avgdf.columns[4:]:
    fig, ax = plt.subplots()
    sns.lineplot(ax=ax, data=avgdf, x='Simulation time (ns)', y=coln, hue='type')
    plt.ylabel('{}'.format(coln))
    plt.legend(bbox_to_anchor=(1, 1))
    plt.xlabel('Time (ns)')
    plt.title('Comparison {} T1 structures'.format(coln))
    plt.show()
