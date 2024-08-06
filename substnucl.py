
# Functions needed to substitute nucleic acid monomers into a given pdb, using MDAnalysis to align the structures and pandas to parse the data

# Used in files like align_subst_mergechains_sirna.py

import os
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align


cwd = os.getcwd()

def get_indeces(atoms_align, structure, residue, chain):
    '''Gets index number for a given atom name in a residue and chain'''
    atom_indeces = {}
    chainid = ord(str(chain))%32 - 1
    topo, _ = structure.topology.to_dataframe()

    for at in atoms_align:
        atom_indeces[at] = topo[(topo['resSeq']==int(residue)) & (topo['chainID']==chainid) & (topo['name']==at)].index[0]

    return atom_indeces

def make_pdb_df(df):
    '''Manually writing lines in pdb format from pd dataframe'''

    lines_file = []

    for i, row in df.iterrows():
        if row['ATOM']=='TER':
            lines_file.append('TER')
            return lines_file
                
        # skipping H atoms
        if 'H' in row['attype']:
            continue
        spac_atnum = (7-len(str(row['atnum'])))*' '
        spac_attype = (4-len(row['attype']))*' ' + (3-len(row['res']))*' '
        space_res = (2-len(row['chain']))*' '

        space_resnum = (4-len(str(row['resnum'])))*' '
        sp_beforex = (12-len(row['x']))*' '
        sp_beforey = (8 - len(row['y']))*' '
        sp_beforez = (8 - len(row['z']))*' '
        if row['elem']==None:
            row['elem']=row['attype'][0]
            
        line_spaced = row['ATOM']+spac_atnum+str(row['atnum'])+'  '+row['attype']+spac_attype+row['res']+space_res+row['chain']+space_resnum+str(row['resnum'])+sp_beforex+row['x']+sp_beforey+row['y']+sp_beforez+row['z']+'  '+row['1']+'  '+row['0']+'           '+row['elem']+' \n'
        lines_file.append(line_spaced)
    return lines_file

def check_delete(string, todelete):
    for word in todelete:
        if word in string:
            return 0
    return 1

dict_rna = {
    'U5':'RU5',
    'A5':'RA5',
    'C5':'RC5',
    'G5':'RG5',
    'U3':'RU3',
    'A3':'RA3',
    'C3':'RC3',
    'G3':'RG3',
    'U':'RU',
    'A':'RA',
    'C':'RC',
    'G':'RG',
    'lA':'LA',
    'lU':'LU',
    'lG':'LG',
    'lC':'LC',
    'lT':'LT',
    'mA':'MA',
    'mU':'MU', 
    'mG':'NG',  ### NG BECAUSE MG is magnesium
    'mC':'MC',
}
def rename_res_rna(res):
    '''Renames residue from name file to names in dictionary'''
    if res in dict_rna:
        return dict_rna[res]
    else:
        return res

def make_df(file, todelete):
    '''Parses pdb into df for ease of modification'''
    file = [line for line in file if check_delete(line, todelete) ]
    file = [line.split() for line in file]


    for n, line in enumerate(file):
        # deleting extra chain added by mdanalysis at before last field
        if len(line) >= 13:
            line.pop(11)
        for elem in line:
            if len(elem) > 8:
                line.pop(10)
        
        # splitting coordinates when there are no spaces between them
        if len(line) < 12:
            for i, element in enumerate(line):
                if element.count('.') > 1:
                    toadd = []
                    for el in element:
                        if el == '.':
                            coord = element[:element.index(el)+4]
                            element = element.replace(coord, '')
                            toadd.append(coord)
                    file[n] = line[:i]+toadd+line[i+1:]

    df = pd.DataFrame(file, columns=['ATOM', 'atnum', 'attype', 'res','chain','resnum','x','y','z','1','0','elem'])
    df['res'] = df['res'].apply(lambda x: rename_res_rna(x))
    return df


def align_subst_single(basestruc, monostruc, atoms_align, ind_temp, chain, keep, todelete, nameoutput, basename, mononame, cwd, ind_ref=1, chain_ref='A', multipleres=False, prime5=False,substnucleotide=False, thiophosphate=False):
    '''Aligns base starting pdb structure with monomer to substitute in and changes in place, returning clean pdb'''

    indeces_base = get_indeces(atoms_align, basestruc, ind_temp, chain)

    indeces_mono = get_indeces(atoms_align, monostruc, ind_ref, chain_ref)

    mobile = mda.Universe(basename)
    ref = mda.Universe(mononame)

    mobile_at = mda.core.groups.AtomGroup(list(indeces_base.values()), mobile)
    ref_at = mda.core.groups.AtomGroup(list(indeces_mono.values()), ref)

    mobileats = 'name'
    for at in atoms_align:
        mobileats = mobileats+' {} or name'.format(at)
    mobileats = mobileats[:-7]

    mobile_string = (mobileats+"and resid {} and segid {}".format( ind_temp, chain))

    refats = 'name'
    if substnucleotide:
        for at in atoms_align:
            refats = refats+' {} or name'.format(at)
    else:
        for at in keep:
            refats = refats+' {} or name'.format(at)
    refats = refats[:-7]
    ref_string = (refats+"and resid {} and segid {}".format(ind_ref, chain_ref))


    alignment = align.AlignTraj(mobile, ref,  filename='{}/aligned_temp.pdb'.format(cwd), prefix=None, weights=None, match_atoms=False,force=True,in_memory=False, reference_atoms = ref_at, mobile_atoms = mobile_at, verbose=True, select={'mobile':mobile_string, 'reference': ref_string})
    alignment.run()


    monostruc.save('mono_temp.pdb', force_overwrite=True)
    # Substituting coordinates 
    with open('aligned_temp.pdb') as f:
        template =  f.readlines()   
    templf = make_df(template[:-1], todelete)

    with open('mono_temp.pdb') as f:
        mono =  f.readlines() 
    monof = make_df(mono, todelete)

    monof.loc[:, 'chain'] = str(chain)

#   from phosphate 2ble struct merger
    if multipleres:
        monof.loc[:, 'resnum'] = monof['resnum'].apply(lambda x: str(int(x)+int(ind_temp)-1))
    else:
        monof.loc[:, 'resnum'] = ind_temp

    #print(monof)

    onlyres = templf[(templf['resnum']==str(ind_temp)) & (templf['chain']==str(chain))]
    if substnucleotide:
        keeprows = onlyres.loc[onlyres['attype'].isin(keep)]
    else:
        keeprows = onlyres.loc[~onlyres['attype'].isin(keep)]
    keeprows.loc[:, 'res'] = monof['res'].unique()[0]

    if thiophosphate:
        if thiophosphate =='right':
            srow = keeprows.loc[keeprows['attype']=='O1P']
            oprow = keeprows.loc[keeprows['attype']=='O2P']
            opind = oprow.index.tolist()[0]
            keeprows.loc[opind, 'attype'] = 'O1P'
            
        else:
            srow = keeprows.loc[keeprows['attype']=='O2P']
        sind= srow.index.tolist()[0]
        keeprows.loc[sind, 'attype'] = 'S'
        keeprows.loc[sind, 'elem'] = 'S'
        monof.loc[:, 'res'] = monof['res'].unique()[0]+'S'
        keeprows.loc[:, 'res'] = monof['res'].unique()[0]

#   from phosphate 2ble struct merger
    if multipleres:
        if prime5:
            keeprows.loc[:, 'resnum'] = str(int(ind_temp) +1) 
    if prime5:
        substituted = pd.concat([monof, keeprows], ignore_index=True)
    else:
        substituted = pd.concat([keeprows, monof], ignore_index=True)

    before = templf[templf['atnum'].astype('int') < int(list(onlyres['atnum'])[0])]

    after = templf[templf['atnum'].astype('int') > int(list(onlyres['atnum'])[-1])]


    if multipleres:
         #modifying number residues after added residue
         after.loc[:, 'resnum'] = (after['resnum'].astype(int)+1).astype('str')

    full = pd.concat([before, substituted, after], ignore_index=True)
    full['atnum'] = full.index.values +1 
    newlines = make_pdb_df(full)

    with open('{}.pdb'.format(nameoutput), "w") as f:
      f.writelines(newlines)

    return '{}.pdb'.format(nameoutput)