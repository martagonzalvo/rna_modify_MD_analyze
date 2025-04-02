# Calculates Helix parameters for each snapshot of the trajectory, for all trajectories (pdb or gro) in folder directory, using DSSR program
# Then gathers all parameters calculated into helixparam.csv

# Usage:

#   python dssrx3dna_folder.py folder dssrdirectory

    # dssrdirectory is location of DSSR program: https://inventions.techventures.columbia.edu/technologies/dssr-an-integrated--CU20391



import os, subprocess, sys, json 
import pandas as pd

folder = sys.argv[1]
dssrdirectory = sys.argv[2]
cwd = os.getcwd()


def getparamsjson(jsonfile):
    '''Defines helical parameters to save in helixparam.csv'''
    with open(jsonfile) as f:
        jsondata = json.load(f)

    params = pd.DataFrame(columns=['shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening','shift', 'slide', 'rise', 'tilt', 'roll', 'twist'])
    pair=0
    for h in jsondata['helices']:
        for p in h['pairs']:
            if 'bp1_params' in p:
                params.loc[pair,['shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening']]=p['bp1_params']
                params.loc[pair,['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']]=p['step_params']
                pair+=1
    return params


# Calculating parameters for all files in folder using DSSR
for file in os.listdir(cwd+'/'+folder):

    print(file)
    if 'pdb' not in file and 'gro' not in file:
        print('skipping')
        continue
    if 'dssr' in file:
        continue

    name = file.split('.')[0]

    if 'gro' in file:
        subprocess.call('''gmx trjconv -f {}.gro -s {}.gro -o {}.pdb<< EOF
        1
        EOF'''.format(name, name, name), cwd=cwd+'/'+folder, shell=True)
    

    subprocess.call('''{} -i={}.pdb -o={}.json --more --json'''.format(dssrdirectory, name, name), cwd=cwd+'/'+folder, shell=True)



    subprocess.call('''gmx trjconv -f {}.pdb -s {}.pdb -o {}_sep.pdb -sep<< EOF
        0
        EOF'''.format(name, name, name), cwd=cwd+'/'+folder, shell=True)
    

    total = [ts for ts in os.listdir(cwd+'/'+folder) if name in ts and 'pdb' in ts and 'sep' in ts]

    for i, subfile in enumerate(total):
        
        if 'sep' not in subfile:
            continue

        newname = subfile.split('.')[0]
        subprocess.call('''{} -i={}.pdb -o={}.json --more --json'''.format(dssrdirectory, newname, newname), cwd=cwd+'/'+folder, shell=True)

        # deleting all but last one
        if i != len(total):
            subprocess.call('''rm {}'''.format(subfile), cwd=cwd+'/'+folder, shell=True)

    
    print('Did ',file)




listparams = []
for jsonfile in os.listdir(cwd+'/'+folder):
    if 'sep' not in jsonfile:
        continue
    if '.json' not in jsonfile:
        continue
    name =  jsonfile.split('.')[0] # These conventions are based on specific naming system
    time = name.split('sep')[-1] 

    params = getparamsjson(cwd+'/'+folder+'/'+jsonfile)
    fileparams = pd.to_numeric(params.mean())
    fileparams['name'] = name
    fileparams['time'] = time
    listparams.append(fileparams)
    print('Did',jsonfile)


allparams = pd.concat(listparams, axis=1) 
allparams= allparams.T

# Reordering columns to desired order
# Time is converted into ns based on the specific dump frequency of snapshots
allparams=pd.concat([allparams['name'],allparams['time'].astype(int)/2, allparams['buckle'],allparams['rise'],allparams['twist'],allparams['opening'],allparams['propeller'],allparams['roll'],allparams['shear'],allparams['shift'],allparams['slide'],allparams['stagger'],allparams['stretch'],allparams['tilt']],axis=1)


allparams.to_csv('{}/helixparam.csv'.format(cwd+'/'+folder))