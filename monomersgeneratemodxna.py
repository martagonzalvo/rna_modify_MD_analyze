# Uses modxna (https://modxna.chpc.utah.edu/about/, https://pubs.acs.org/doi/10.1021/acs.jctc.4c01164) to generate topology and template monomer structures of all modified nucleotides to use in future simulations

# usage: in folder where want to generate all files
#    python monomersgeneratemodxna.py

# Takes almost 1 minute (0.8 min) for 28 residues below


import os, subprocess, time
import parmed as pmd

start = time.time()

monomers = [ ["DPO LNA RAA", "LA"], # L is locked nucleic acid 
["DPO LNA RGG", "LG"], 
["DPO LNA RCC", "LC"],
["DPO LNA DTT", "LT"],
["DPO LNA M5C", "LH"], # H is 5-Methylcytosine
["DPO OME RAA", "MA"], # M is 2’-O-methyl ribose
["DPO OME RCC", "MC"],
["DPO OME RGG", "NG"], # NG is 2’-O-methyl ribose guanine, to differentiate from ion Mg
["DPO OME RUU", "MU"],
["DPO AF2 RAA", "FA"], # F is 2’-fluororibose
["DPO AF2 RCC", "FC"],
["DPO AF2 RGG", "FG"],
["DPO AF2 RUU", "FU"],
["PS1 AF2 RGG", "FGS"], # S is phosphorothioate
["PS1 LNA DTT", "LTS"],
["PS1 OME RAA", "MAS"],
["PS1 OME RGG", "NGS"],
["PS1 OME DTT", "MTS"],
["PS1 OME M5C", "MHS"],
]


atomstodel = ["P", "OP1", "OP2", "S1"]

cwd = os.getcwd()

for data in monomers:
    monolist, name = data
    with open("{}/{}.modxna".format(cwd,name), 'w') as f:
        f.writelines(monolist+'\n')
    
    print("echo modxna.1.6.sh -i {}.modxna -m {}".format(name, name))
    subprocess.call("bash modxna.1.6.sh -i {}.modxna -m {}".format(name, name), cwd=cwd, shell=True)

    restemplate = pmd.load_file('{}.lib'.format(name))
    restemplate['{}'.format(name)].save('{}.pdb'.format(name), overwrite=True)

    # also deleting Phosphate including O5' and O3' due to how substitute monomer into structure

    with open('{}.pdb'.format(name), 'r') as f:
        pdblines = f.readlines()

    newlines = []
    for line in pdblines:
        save = True
        for at in atomstodel:
            if at in line:
                save = False
        if save==True:
            if 'HETATM' in line:
                spaces = (3-len(name))*' '
                line = line[:17]+spaces+'{} A   '.format(name)+line[25:]

            newlines.append(line)

    with open('{}.pdb'.format(name), 'w') as f:
        f.writelines(newlines)


print('Done in this many minutes, 18 unique residues ',str((time.time()-start)/60))


