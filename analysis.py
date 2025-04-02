
# Usage:
#   python analysis_sirna.py struct.pdb templatefolder

# struct.pdb is trajectory pdb with all timesteps in a single pdb file
# templatefolder is folder with template for pdb with correct chain information and other topology used to define the indices of the measures of interest (e.g. line 120 of this file). The file that is output of the gromacs command "gmx pdb2gmx -f duplex.pdb -o duplexgmx.pdb ..." is usually a good template, e.g. in this case duplexgmx.pdb

# Outputs csv: summary_namepdb.csv

# Distances of interest are described in line 120-151 of this file. For a different system, they would change.


import pandas as pd
import mdtraj as md
import numpy as np
from sklearn.decomposition import PCA
import sys, os, time

start = time.time()
typesfiles = [ 
   'md',
   ]

folder = sys.argv[1]
templatefolder = sys.argv[2]

def calc_coord(posit_name, listnames, full_traj, indeces,  summary_df, nth=1):
    ''''same as posit but saving coords, not doing linalg
    '''
    traj = md.Trajectory.atom_slice(full_traj,indeces).xyz

    for i, data  in enumerate(zip(traj[::nth],listnames[::nth])):
        group, name = data
        if 'noallos' not in name: 
            positval = np.mean(group, axis=0)
            summary_df.loc[summary_df.name==name, '{}_x'.format(posit_name)] = positval[0]
            summary_df.loc[summary_df.name==name, '{}_y'.format(posit_name)] = positval[1]
            summary_df.loc[summary_df.name==name, '{}_z'.format(posit_name)] = positval[2]
    
    return summary_df


def calc_dist(distname, listnames, traj, indexg1, indexg2, summary_df, nth=1):
    group1_traj = md.Trajectory.atom_slice(traj,indexg1).xyz
    group2_traj = md.Trajectory.atom_slice(traj,indexg2).xyz
    
    # if only 1 atom, use mdtraj distance to account for boundary conditions correctly
    if len(indexg1) == 1 and len(indexg2) ==1:
        distsmd = md.compute_distances(traj, np.array([indexg1, indexg2]).T, periodic=True, opt=True)
        summary_df['{}mdtraj'.format(distname)] = distsmd
    else:
        for i, data  in enumerate(zip(group1_traj[::nth], group2_traj[::nth],listnames[::nth])):
                group1, group2, name = data

                distval = np.linalg.norm(np.mean(group1, axis=0)-np.mean(group2, axis=0))
                summary_df.loc[summary_df.name==name, '{}'.format(distname)] = distval
    
    return summary_df

def calc_vector_pca(posit):
    if len(posit)==1:
        vector = posit[0]
    elif len(posit)==2:
        vector = posit[1]-posit[0]
    else:
        pca = PCA(n_components=3)
        pca.fit(posit)
        vector=pca.singular_values_
    return vector

def angle(posit1, posit2, posit3):

    vector1 = calc_vector_pca(posit1)
    vector2 = calc_vector_pca(posit2)

    if posit3 is not None: ### CHECK THIS WORKS
        vector3 = calc_vector_pca(posit3)  

        vector1 = np.mean(vector1, axis=0)-np.mean(vector2, axis=0)
        vector2 = np.mean(vector3, axis=0)-np.mean(vector2, axis=0)

    unitv1 = vector1 / np.linalg.norm(vector1)
    unitv2 = vector2 / np.linalg.norm(vector2)
    angle =  np.degrees(np.arccos(np.clip(np.dot(unitv1, unitv2), -1.0, 1.0)))

    return angle

def calc_angle(angle_name, listnames, traj, indexg1, indexg2, indexg3, summary_df, nth=1):

    group1_traj = md.Trajectory.atom_slice(traj,indexg1).xyz
    group2_traj = md.Trajectory.atom_slice(traj,indexg2).xyz
    coord1 = group1_traj - np.mean(group1_traj,axis=0)
    coord2 = group2_traj - np.mean(group2_traj,axis=0)

    if  indexg3 is not None: 
        group3_traj = md.Trajectory.atom_slice(traj,indexg3).xyz
        coord3 = group3_traj - np.mean(group3_traj,axis=0)
    else:
        coord3 = [None]*len(group2_traj)    
    
    for i, data  in enumerate(zip(coord1[::nth], coord2[::nth], coord3[::nth],listnames[::nth])):
        coord1, coord2, coord3, name = data
        angleval = angle(coord1, coord2, coord3)
        summary_df.loc[summary_df.name==name, '{}'.format(angle_name)] = angleval
    return summary_df


cwd = os.getcwd()

path = cwd+'/'+folder


for trajpdb in os.listdir(path):
    print(trajpdb)
    if 'pdb' not in trajpdb:
        continue
    for type in typesfiles:
        if type not in trajpdb:
            print(type, 'cont')
            continue
        elif "copy" in trajpdb:
            print(type, 'cont')
            continue
        elif type in trajpdb:
            print('Loading trajectory')
            print(trajpdb)

            traj = md.load(path+'/{}'.format(trajpdb))
            trajtempl = md.load('{}'.format(templatefolder+'/'+type+'.pdb'))

            topo, _ = trajtempl.topology.to_dataframe()
            print('Loaded trajectory')
            namesyst = trajpdb.split('.')[0]
            listnames = ['{}_'.format(namesyst)+str(i) for i in range(len(traj))]

            res_endoverhang = topo['resSeq'].max()

            # List of indeces involved in each measure- depends on structure analyzing!
            ind_overhang = topo[(topo['name']=="C3'") & (topo['resSeq'] > 21) &  (topo['chainID']==0)].index
            ind_endoverhang = topo[(topo['name']=="O3'") & (topo['resSeq'] == res_endoverhang) &  (topo['chainID']==0)].index
            ind_topduplex = topo[(topo['name']=="C1'") & ((topo['resSeq'] == -22)) & (topo['chainID']==1)].index

            ind_b_12 = topo[(topo['name']=="C3'") & (topo['resSeq'] == 12) & (topo['chainID']==1)].index
            ind_b_11 = topo[(topo['name']=="C3'") & (topo['resSeq'] == 11) & (topo['chainID']==1)].index

            ind_allnucleic = topo[(topo['chainID']==0) | (topo['chainID']==1) | (topo['chainID']==2)].index
            ind_helixlong = topo[((topo['chainID']==0)& (topo['resSeq'] < 22)) | ((topo['chainID']==1)& (topo['resSeq'] < 22)& (topo['resSeq'] > 1))].index
            ind_helixshort = topo[((topo['chainID']==1)& (topo['resSeq'] < 0)) | ((topo['chainID']==2))].index 


            list_distances = [
                    ['overhang-topduplex',ind_topduplex,ind_overhang],
                    ['topoverhang-topduplex', ind_topduplex,ind_endoverhang],
                    ['gapcore_dist1112', ind_b_11,ind_b_12],
                    ['helices_dist', ind_helixlong,ind_helixshort],
            ]

            list_coords= [        
                    ['endoverhang', ind_endoverhang],
                    ['ind_b_12',ind_b_12],
                    ['ind_b_11',ind_b_11],
            ]

            list_angles = [
                ['helices_ang', ind_helixlong,ind_helixshort, None],
                ['overhang-helix_ang',ind_helixlong, ind_overhang, None],
            ]

            print('Analyzing and saving data')

            summary_df =pd.DataFrame(listnames, columns=['name'])


            for i, data in enumerate(list_distances): 
                distname, indexg1, indexg2 = data
                summary_df = calc_dist(distname, listnames, traj, indexg1, indexg2, summary_df, nth=1)

            for i, data in enumerate(list_angles):
                positname, indexg1, indexg2, indexg3 = data
                summary_df = calc_angle(positname, listnames, traj, indexg1, indexg2, indexg3,  summary_df, nth=1)


            summary_df['rmsd'] = md.rmsd(traj, traj, frame=0, atom_indices=ind_allnucleic)
            summary_df['rmsd_overhang'] =  md.rmsd(traj, traj,frame=0, atom_indices=ind_overhang)


            summary_df.to_csv('{}/summary_{}.csv'.format(path,namesyst))
            
            print('Done, file in summary_{}.csv'.format(namesyst))



print('Done all in {} mins'.format((time.time()-start)/60))

