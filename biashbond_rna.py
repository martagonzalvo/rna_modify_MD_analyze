# Adds biases to modified mdp file for gromacs simulation given .gro or .pdb file

# Adds:
# - H-bond biases each base pair, keeping 2/3 h-bonds in place and bases pairing correctly

# Can also add:
# - additional constraints: fix gaps, merge ends closer
# - backbonects: O3'-P all backbones
# - repulsive cts: beyond distance e.g 10A, between 2 atoms

# Returns:
# -  new {md}_const.mdp and index.ndx
# - (optional) index being made off dup_ions.gro (from genion ions.tpr, bc groups with water and ions), and atom numbers selected from duplexgmxer.pdb (from pdb2gmx, bc chains) - both need to exist


# Usage:
#   python biashbond_rna.py structbase.pdb template.mdp (no) structype True True True False

#   python biashbond_rna.py duplexgmxmer.pdb em.mdp (no) ev3 True True True False

# Inputs:
#   - duplexgmxer.pdb is structure to base indexes and atoms off, preferably a pdb. Can use before adding waters, since they are added after the protein indeces
#   - em.mdp is mdp to use as template
#   - (no) is for not modifying index.ndx, after creating em.mdp, just prints mdps without modifying index.ndx
#   - structype, e.g. ev3 is type of structure: for each type of structure, the base pair and atom constraint indices are going to change. need to choose from list below or create new indeces
#   - True1 is for base pair constraints
#   - True2 is add additional cts: fix gaps, merge ends closer etc
#   - True3 is add backbone cts: O3'-P
#   - False is add repulsive cts, for overhang

import mdtraj as md, numpy as np
import sys, subprocess, os

# structure as reference
struct = sys.argv[1]
# mdp file to add constraints: em, nvt, annealing etc
mdpfile = sys.argv[2]
# create and modify index file with atom indeces
modindex = True
if sys.argv[3] =='no' or sys.argv[3] =='NO' or sys.argv[3] =='No':
     modindex = False

# type of structure
structype = sys.argv[4]

# constraints to apply
bpcts = False
additcts = False
backbonect = False
repulsct = False


if sys.argv[5]=='True':
    bpcts = True
if sys.argv[6] =='True':
    additcts = True
if sys.argv[7] =='True':
    backbonect = True
if sys.argv[8] =='True':
    repulsct = True


structuretypes=[
'ogdup',
'long2dup',
'long6dup',]


if structype not in structuretypes:
    print('''Please choose one of the following structure types or add new atom indexes: {}
    
          Usage example: 
          
          python biashbond_rna.py duplexgmxmer.pdb em.mdp no ev3

          EXITING, NOT GENERATING NEW MDPS
'''.format(structuretypes))
    exit(0)

# -1 bc rest are 1-indexed, but I already accounted for that


additcts_list = []
repulsivects_list = []


if structype == 'ogdup':
    strands = [
        ['A', np.arange(1,27)],
        ['B', np.concatenate([np.arange(11,-1),np.arange(-22,0),np.arange(22,11)])],
        ['C', np.arange(1,21)],
        ]
    additcts_list = [[[0.16, 1915-1, 1907-1]], # C3 RL 24 and O3' RC 21, chain B
        [[0.16, 1256-1, 1257-1]], # RU 0 O3'- P RC 1, chain B
        [[0.16, 2220-1, 2219-1]], #  RU -12 P - RU -13 O3', chain B
        [[0.25, 2249-1 ,899-1]], #  RU -12 O3' - RU -11 O5', chain B
        ]

hbond_pairs = [    
    [('C', -1), ('B', 1)],
    [('C', -2), ('B', 2)],
    [('C', -3), ('B', 3)],
    [('C', -4), ('B', 4)],
    [('C', -5), ('B', 5)],
    [('C', -6), ('B', 6)],
    [('C', -7), ('B', 7)],
    [('C', -8), ('B', 8)],
    [('C', -9), ('B', 9)],
    [('C', -10), ('B',10)],
    [('C', -11), ('B',11)],
    [('C', -12), ('B',12)],
    [('C', -13), ('B',13)],
    [('C', -14), ('B',14)],
    [('C', -15), ('B',15)],
    [('C', -16), ('B',16)],
    [('C', -17), ('B',17)],
    [('C', -18), ('B',18)],
    [('C', -19), ('B',19)],
    [('C', -20), ('B',20)],
    [('C', -21), ('B',21)],
    [('A', 1), ('B', -1)],
    [('A', 2), ('B', -2)],
    [('A', 3), ('B', -3)],
    [('A', 4), ('B', -4)],
    [('A', 5), ('B', -5)],
    [('A', 6), ('B', -6)],
    [('A', 7), ('B', -7)],
    [('A', 8), ('B', -8)],
    [('A', 9), ('B', -9)],
    [('A', 10), ('B',-10)],
    [('A', 11), ('B',-11)],
    [('A', 12), ('B',-12)],
    [('A', 13), ('B',-13)],
    [('A', 14), ('B',-14)],
    [('A', 15), ('B',-15)],
    [('A', 16), ('B',-16)],
    [('A', 17), ('B',-17)],
    [('A', 18), ('B',-18)],
    [('A', 19), ('B',-19)],
    [('A', 20), ('B',-20)],
    [('A', 21), ('B',-21)],
    [('A', 22), ('B',-22)],
    ]


CV_ATOMS = {
    'CG': [
        (1/6, 'H41', 'O6'),
        (1/6, 'H42', 'O6'),
        (1/6, 'O2', 'H21'),
        (1/6, 'O2', 'H22'),
        (1/3, 'N3', 'H1')
    ],
    'GU': [
        (1/2, 'H1', 'O2'),
        (1/2, 'O6', 'H3')
    ],
    'AU': [
        (1/4, 'H61', 'O4'),
        (1/4, 'H62', 'O4'),
        (1/2, 'N1', 'H3')
    ]
}

MD = '''
title                   = RNA NPT equilibration 
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 500000    ; 2 * 500000 = 1000 ps (1 ns)
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 0         ; suppress bulky .trr file by specifying 
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps
compressed-x-grps       = System    ; save the whole system
; Bond parameters
continuation            = yes       ; Restarting after NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighborsearching
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = RNA Water_and_ions   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Velocity generation
gen_vel                 = no        ; Velocity generation is off 

'''

PULL = '''
pull                    = yes
pull_ngroups            = {groups}  ;
pull_ncoords            = {coords}  ;
{names}
{defs}
pull-nstxout            = {save_cv}
pull-nstfout            = 0 ; dont think we need the force, its implied by the position

'''

PULL_DEF = '''
pull-coord{n}-type         = umbrella
pull-coord{n}-geometry     = distance
pull-coord{n}-groups       = {i} {j}
pull-coord{n}-dim          = Y Y Y
pull-coord{n}-k            = 0      ; avoid forces working directly on this distance

'''

PULL_REPULS = '''
pull-coord{n}-type         = umbrella
pull-coord{n}-geometry     = transformation
pull-coord{n}-expression   = (1 - 1/(1+exp(50*({xnum}-{x_nm})))  ) ; x1 refers to the value of coord1 
pull-coord{n}-init         = 1 ; restrains to beyond distance
pull-coord{n}-k            = 10000

'''

PULL_DIST = '''
pull-coord{n}-type         = umbrella
pull-coord{n}-geometry     = distance
pull-coord{n}-groups       = {i} {j}
pull-coord{n}-dim          = Y Y Y
pull-coord{n}-k            = 200000
pull-coord{n}-init         = {d}    

'''

PULL_FINAL = '''
pull-coord{n}-type         = umbrella
pull-coord{n}-geometry     = transformation
pull-coord{n}-expression   = {expr}  ; x1 and x2 refer to the value of coord1 and coord2
pull-coord{n}-init         = {x_nm}    ; restrains the average distance
pull-coord{n}-k            = {k_kj_nm2}

'''



def get_index_pairs(hbond_pairs, topo):
    pair_list = []
    for pair in hbond_pairs:
        # getting names bases
        bases = ''
        bp = []
        chains = []
        for base in pair:
            chain, res = base
            chainid = ord(str(chain))%32 - 1 # from letter to number 
            resname = topo[(topo['resSeq']==res)& (topo['chainID']==chainid)]['resName'].unique()[0]

            bp.append(res)
            chains.append(chainid)

            if 'A' in resname:
                bases = bases+'A'
            elif 'U' in resname or 'T' in resname:
                bases = bases+'U'
            elif 'C' in resname:
                bases = bases+'C'
            elif 'G' in resname:
                bases = bases+'G'
        # ordering base pair alphabetically
        if bases[1] < bases[0]:
            bp, chains, bases = bp[::-1], chains[::-1], bases[::-1]
        # getting atom indices for each h-bonding atom pair
        atomnames = CV_ATOMS[bases]

        pairs = []
        for atompairs in atomnames:
            atoms1 = topo[(topo['resSeq']==bp[0]) & (topo['chainID']==chains[0]) & (topo['name']==atompairs[1])].index[0]
            atoms2 = topo[(topo['resSeq']==bp[1]) & (topo['chainID']==chains[1]) & (topo['name']==atompairs[2])].index[0]
            coeff = [atompairs[0], atoms1, atoms2]
            pairs.append(coeff)
        
        pair_list.append(pairs)
    return pair_list

def get_index_pairs_backbone(strands, topo):
    pair_list = []
    for namestrand, residuenums in strands: 
        for r, res in enumerate(residuenums):
            if r==0:
                continue
            chainid = ord(str(namestrand))%32 - 1 # from letter to number 
            pind = topo[(topo['chainID']==chainid) & (topo['resSeq']==res) & (topo['name']=="P")]
            o3ind = topo[(topo['chainID']==chainid) & (topo['resSeq']==(res-1)) & (topo['name']=="O3'")]
            if not pind.empty and not o3ind.empty:
                pair_list.append([[0.16, pind.index[0], o3ind.index[0]]])
    return pair_list


def generate_pull(pair_lists, save_cv, additcts_list,backbone_pairlists,repulsivects_list, bp_cts=True, additcts=True, backbonect=True, repulsct=False):
    '''pairs is a list of (c, i, j) with c the coefficient of the distance between atom indices i and j'''
    indices = [k+1 for pl in pair_lists for (c, i, j) in pl for k in (i, j)] # gromacs indices are one based!
    z = 1
    groups = 1
    defs = ''
    keys = []

    # CREATING H-BOND BASE PAIR CONSTRAINTS
    if bp_cts:
        for w, pairs in enumerate(pair_lists):
            if len(pairs) == 2:
                expr=' + '.join(f'{c} * x{z+n}' for n, (c, _, _) in enumerate(pairs))
            elif len(pairs) == 3:
                expr = f'0.5 * min(x{z+0}, x{z+1}) + 0.5 * x{z+2}'
            else:
                expr = f'{1/3} * min(x{z+0}, x{z+1}) + {1/3} * min(x{z+2}, x{z+3}) + {1/3} * x{z+4}'

            defs += ''.join(PULL_DEF.format(n=z+n, i=groups+2*n, j=groups+2*n+1) for n in range(len(pairs)))
            #defs += ''.join(PULL_DEF.format(n=z+n, i=z+2*n, j=z+2*n+1) for n in range(len(pairs)))

            defs += PULL_FINAL.format(n=z+len(pairs), k_kj_nm2=int(2e5), x_nm=0.35, expr=f'max(0.35, {expr})')
            keys += [(w, (i, j)) for _, i, j in pairs] + [(w, None)]

            z += len(pairs) + 1
            groups += len(pairs)*2

        coords = sum(len(pairs)+1 for pairs in pair_lists)
    else:
        coords = 0
        w = 0
    
    # ADDITIONAL DISTANCES - GAP CHAINS, C3', MA
    if additcts:
        for v, pairs in enumerate(additcts_list):
            d, i, j= pairs[0]
            defs+=''.join(PULL_DIST.format(n=z+n, i=groups+2*n, j=groups+2*n+1, d=d) for n in range(len(pairs)))
            keys += [(v+w, (i, j)) for _, i, j in pairs] + [(v+w, None)]
            z += len(pairs)
            groups += len(pairs)*2

        # changing pairs for correct number of coord!        
        indices = indices+[k+1 for pl in additcts_list for (_,i,j) in pl for k in (i,j)]
    
        coords = coords + len(additcts_list)
    
    if backbonect:
        for v, pairs in enumerate(backbone_pairlists):
            d, i, j= pairs[0]
            defs+=''.join(PULL_DIST.format(n=z+n, i=groups+2*n, j=groups+2*n+1, d=d) for n in range(len(pairs)))
            keys += [(v+w, (i, j)) for _, i, j in pairs] + [(v+w, None)]
            z += len(pairs)
            groups += len(pairs)*2

        # changing pairs for correct number of coord!        
        indices = indices+[k+1 for pl in backbone_pairlists for (_,i,j) in pl for k in (i,j)]
    
        coords = coords + len(backbone_pairlists)
    
    if repulsct:
        for v, pairs in enumerate(repulsivects_list):
            d, i, j=pairs[0]
            defs += ''.join(PULL_DEF.format(n=z+n, i=groups+2*n, j=groups+2*n+1) for n in range(len(pairs)))

            defs += PULL_REPULS.format(n=z+len(pairs), xnum=z , x_nm=d)
            keys += [(v+w, (i, j)) for _, i, j in pairs] + [(v+w, None)]
            z += len(pairs)
            groups += len(pairs)*2


        # changing pairs for correct number of coord!        
        indices = indices+[k+1 for pl in repulsivects_list for (_,i,j) in pl for k in (i,j)]
    
        coords = coords + len(repulsivects_list)*2

    return indices, keys, PULL.format(groups=len(indices), 
        coords=coords, save_cv=save_cv,
        names='\n'.join('pull_group{}_name = a_{}'.format(1+n, k) for n, k in enumerate(indices)),
        defs=defs)


### RUNNING HERE

cwd = os.getcwd()

gro = md.load(struct)
topo, _ = gro.topology.to_dataframe()


# pairs of atoms in H-bonds in base pairs
pair_lists_hbond = get_index_pairs(hbond_pairs, topo)
print('Got h-bond pair lists')

# backbone constraints
pairs_backbone = get_index_pairs_backbone(strands, topo)
print('Got backbone pair lists')

# generate all forces and constraints
indices, keys, pull = generate_pull(pair_lists_hbond, 20, additcts_list, pairs_backbone,repulsivects_list,bp_cts=bpcts, additcts=additcts, backbonect=backbonect, repulsct=repulsct)
print('made pull mdp section')

# 10 is index number, should be higher than existing number of groups - ESSENTIAL!!
stringnames = '\n'.join(f'a {a}\nname {10+n} a_{a}' for n, a in enumerate(indices))

if modindex:
    subprocess.call("""gmx make_ndx -f dup_ions.pdb -o index.ndx << EOF
    {}
    q
    EOF""".format(stringnames), cwd=cwd, shell=True)

if mdpfile:
    with open(mdpfile, 'r') as f:
        mdplines = f.read()
    md_constraints = mdplines +pull
    namefile = mdpfile.split('.')[0]
else:
    md_constraints = MD+pull
    namefile = 'md'



with open('{}_ct.mdp'.format(namefile), 'w') as f:
    f.writelines(md_constraints)

print('Made file with h-bond and closing duplex constants, {}_const.mdp'.format(namefile))
