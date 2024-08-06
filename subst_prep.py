
# Creates 1 folder per simulation in larger folder, substitutes seq into pymol duplex, prepares all gromacs files, adds H-bond constraints to .mdp files, 
# Input json with all necessary data per simulation, pointing to location pymol duplex and seq substitute
# folderfiles is files with starting pdb and txt files describing subsitutions


# Usage:
#       python subst_prep.py simulation.json folderfiles

import sys, subprocess, os, time, json

listsimulat = sys.argv[1]
folderfiles = sys.argv[2]

path = os.getcwd()+'/'+folderfiles

with open('{}'.format(listsimulat)) as f:
    datasimulat = json.load(f)

starttime = time.time()


if os.path.exists(path+'/prepared')==False:
    os.mkdir(path+'/prepared')
else:
    # Delete WHOLE FOLDER in case of previous iterations same folder
    print('Deleting everything in folder prepared')
    subprocess.call('rm -rf prepared/*', cwd=path, shell=True)



def prepare_gmx(cwd, mgions):
    '''Run all steps before equilibration in gromacs: solvation, add ions'''
    # prep file gromacs
    subprocess.call('''gmx pdb2gmx -f substituted.pdb -o duplexgmxmer.pdb -water tip3p -ignh -ff amber14sblnaintorna -missing
    gmx editconf -f duplexgmxmer.pdb -o dup_boxmer.pdb -c -d 1 -bt cubic 
    gmx solvate -cp dup_boxmer.pdb -o dup_solvmer.gro -p topol.top 
    gmx grompp -f ions.mdp -c dup_solvmer.gro -p topol.top -o ions.tpr -maxwarn 2''', cwd=cwd, shell=True)

    # add charge
    subprocess.call('''gmx genion -s ions.tpr -o dup_ca.gro -p topol.top -pname CA -np {} << EOF  
    3
    EOF'''.format(mgions), cwd=cwd, shell=True)

    if mgions%2 == 0:
        na=0
        SOLlast = 4
    else:
        na=1
        SOLlast = 5

    subprocess.call('''gmx grompp -f ions.mdp -c dup_ca.gro -p topol.top -o ions.tpr -maxwarn 2
    gmx genion -s ions.tpr -o dup_na.gro -p topol.top -pname NA -np {} << EOF  
    {}
    EOF'''.format(na, SOLlast-1), cwd=cwd, shell=True)

    subprocess.call('''gmx grompp -f ions.mdp -c dup_na.gro -p topol.top -o ions.tpr -maxwarn 2
    gmx genion -s ions.tpr -o dup_ions.gro -p topol.top -pname NA -nname CL -conc 0.154 << EOF
    {}
    EOF'''.format(SOLlast), cwd=cwd, shell=True)

    subprocess.call("rm -rf \#*", cwd=cwd, shell=True)


filesupdate = ['nvt', 'md']
def makectmdp(cwd, basestrucct, bpcts, additcts, backbonects, repulsects):
    '''Add h-bond and other constraints to each .mdp file depending on structure'''
    
    # MAKING PDB TO READ TO GET COORD TO ADD CONSTRAINTS LATER
    subprocess.call('''gmx editconf -f dup_ions.gro -o dup_ions.pdb''', cwd=cwd, shell=True)

    # ADDING CONSTRAINTS
    # only em change index
    subprocess.call('''python ~/rna_md_analysis/biashbond_rna.py duplexgmxmer.pdb em.mdp yes {} {} {} {} {}'''.format(basestrucct, bpcts, additcts, backbonects, repulsects), cwd=cwd, shell=True)

    # rest don't change index

    for file in filesupdate:
        subprocess.call('''python ~/rna_md_analysis/biashbond_rna.py duplexgmxmer.pdb {}.mdp no {} {} {} {} {}'''.format(file, basestrucct, bpcts, additcts, backbonects, repulsects), cwd=cwd, shell=True)


for sim in datasimulat:
    
    print(sim)
    name = datasimulat[sim]['name']
    basestrucct = datasimulat[sim]['basestrucct']
    initduplex = datasimulat[sim]['initduplex']
    mgions = datasimulat[sim]['mgions']
    seqsubst = datasimulat[sim]['seqsubst']
    bpcts = datasimulat[sim]['bpcts']
    additcts = datasimulat[sim]['additcts']
    backbonects = datasimulat[sim]['backbonects']
    repulsects =datasimulat[sim]['repulsects']
    namehpc = datasimulat[sim]['namehpc']

    
    # Making 1 folder/simulation 
    cwd = path+'/prepared/'+name
    os.mkdir(cwd)


    subprocess.call('''cp -r ~/rna_md_analysis/samplefiles/amber14sblnaintorna.ff .
    cp ~/rna_md_analysis/samplefiles/*mdp .
    cp ~/rna_md_analysis/samplefiles/*sh .''', cwd=cwd, shell=True)

    # Substituting unnat residues into pymol duplex
    subprocess.call('''python ~/rna_md_analysis/substnucl_wholefile.py ~/rna_md_analysis/{}/{} ~/rna_md_analysis/{}/{} ~/rna_md_analysis/{}/monomer_struct substituted.pdb'''.format(folderfiles,initduplex, folderfiles,seqsubst,folderfiles), cwd=cwd, shell=True)


    # preparing gmx
    prepare_gmx(cwd, mgions)

    # make add constants to mdp files
    makectmdp(cwd, basestrucct, bpcts, additcts, backbonects, repulsects)


    print('Added constraints to all mdp files, files to run are the modified following files with {} _ct.mdp'.format(filesupdate))


timeend = (time.time()-starttime)/60
print('It took {} mins'.format(timeend))

print('Prepared folders are in folder prepared')
print('Will delete everything in prepared if script run again!')


