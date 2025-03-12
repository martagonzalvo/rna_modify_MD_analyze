

# Usage:

# python parameterize_gmx_modxna.py substituted

    # substituted is namefile



import subprocess, os, sys
import parmed as pmd


path =  os.getcwd()

namefile = sys.argv[1].split('.')[0]


subprocess.call('''cp ~/*.lib  .''', cwd=path, shell=True)
subprocess.call('''cp -r ~/dat  .''', cwd=path, shell=True)


tleapfile = '''loadamberparams dat/frcmod.modxna
source leaprc.DNA.OL15
source leaprc.water.tip3p
loadoff LG.lib
loadoff LC.lib
loadoff LT.lib
loadoff LH.lib
loadoff EG.lib
loadoff MA.lib
loadoff MC.lib
loadoff NG.lib
loadoff MU.lib
loadoff FA.lib
loadoff FC.lib
loadoff FG.lib
loadoff FU.lib
loadoff LAS.lib
loadoff LGS.lib
loadoff LCS.lib
loadoff LTS.lib
loadoff LHS.lib
loadoff EAS.lib
loadoff ECS.lib
loadoff EGS.lib
loadoff ETS.lib
loadoff EHS.lib
loadoff MAS.lib
loadoff MCS.lib
loadoff NGS.lib
loadoff MUS.lib
loadoff FUS.lib

{name} = loadpdb {name}.pdb



saveamberparm {name} {name}.topo {name}.coords
savepdb {name} amber{name}.pdb
quit
    '''.format(name=namefile)

with open('tleap.in', 'w') as f:
    f.writelines(tleapfile)



subprocess.call('''tleap -s -f tleap.in > testtleap.out''', cwd=path, shell=True)



res = pmd.load_file('{}.topo'.format(namefile))
res.save('modrna.top',overwrite=True)


# Parses topology generated and splitting into as many itp files as molecules
with open('modrna.top', 'r') as f:
    topfile = f.readlines()

itplines = ""
namemols = ""
savemol = False
prevline = None
for i,line in enumerate(topfile):

    if '[ moleculetype ]' in line:
        molsave = []
        savemol = True

    if savemol==True:
    
        if 'system' in line and '[ system ]' not in line:
            name = line.split(' ')[0]

        # if two consecutive lines empty, molecule is done
        if (line == "\n" and prevline =="\n") or '[ system ]' in line:
            with open('{}.itp'.format(name), 'w') as f:
                f.writelines(molsave)
            savemol=False
            itplines = itplines + """#include "{}.itp"\n""".format(name)
            namemols = namemols + "{}         1\n".format(name)

        molsave.append(line)
    prevline = line


itplines = itplines[:-1] # no space at end


topology = '''; Include forcefield parameters
#include "./amber14sb.ff/forcefield.itp"

; Include chain topologies
{itps}

; Include water topology
#include "./amber14sb.ff/tip3p.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include topology for ions
#include "./amber14sb.ff/ions.itp


[ system ]
; Name
ModsiRNA

[ molecules ]
; Compound        #mols
{namemols}'''.format(itps=itplines, namemols=namemols)


with open('topol.top', 'w') as f:
    f.write(topology)





