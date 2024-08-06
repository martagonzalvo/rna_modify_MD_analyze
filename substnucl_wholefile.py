# insert ALL OF THE modified residues on rosetta-generated double-stranded rna

#   python substnucl_wholefile.py template.pdb to_subst.txt monomerstructfolder nameoutput.pdb


#   - template is  pdb
#   - to_subst.txt is file formatted: N C modX
#     separated by single spaces!
#     N, e.g. 3 is number of nucleotide in template where modified structure will be placed into
#     C, e.g. A is chain where found
#     modX, e.g. lC is name of modified nucleotide
#       Example:
#           3 A lC
#           7 C mU
#           NO SPACES AT THE END of to_subst.txt!
#   - monomerstructfolder is folder with pdb structures of monomers to substitute in
#   - nameoutput.pdb is desired final name



import mdtraj as md
import os, sys
from substnucl import *

# set variables - check same as substnucl.py!
atoms_align = ["C3'", "C4'", "C5'"] # , "C1'"
todelete = ['HEADER', 'EXPDTA', 'REMARK', 'MODEL', 'END', 'TER', '~', 'N_BS', 'N_NWC', 'N_WC', 'score', 'TITLE', 'CRYST1']
keep = ['P', 'O2P', 'O1P', "O3'", "O5'"]

template = sys.argv[1]
mono_file = sys.argv[2]
monostructfolder = sys.argv[3]
nameoutput = sys.argv[4]

listmonomers = os.listdir(monostructfolder)
basestruc = md.load(template)

cwd = os.getcwd()

intermoutput = 'subst_interm'

with open(mono_file) as f:
    modify_positions = f.readlines()

erred = []

print('starting, takes about 1-2 mins for ~96 residues')

for line in modify_positions:
    data = line.split('\n')[0]
    position, chain, monomer =  data.split(' ')
    thiophosphate = False

    if 'sr' in monomer or 'Sr' in monomer or 'SR' in monomer:
        thiophosphate = 'right'
        monomer = monomer[:-2]
    if 'sl' in monomer or 'Sl' in monomer or 'Sl' in monomer:
        thiophosphate = 'left'
        monomer = monomer[:-2]

    if monomer+'.pdb' not in listmonomers:
        print("Can't substitute {}, don't have monomer structure".format(data))
        continue

    else:
        index_template = int(position)
        try: 
            monostruc = md.load(monostructfolder+'/'+monomer+'.pdb')
            newtemplate = align_subst_single(basestruc, monostruc, atoms_align, index_template, chain, keep, todelete, intermoutput, template,monostructfolder+'/'+monomer+'.pdb', cwd, ind_ref=1, chain_ref='A', multipleres=False, prime5=False,substnucleotide=True, thiophosphate=thiophosphate) 
            # UPDATE NAME SO ITERATE AND SAVE OVER IT
            template = 'subst_interm.pdb'
            basestruc = md.load(newtemplate)
            print('Substituted {}'.format(data))

        except Exception as e:
            print("Couldn't do substitution {}".format(data))
            erred.append(data)
            print(e)
            break



os.rename(intermoutput+'.pdb', nameoutput)

print('Substituted all monomers defined in {} into {}, output is file {}'.format(mono_file, template, nameoutput))
print("Couldn't do these substitutions: {}".format(erred))



