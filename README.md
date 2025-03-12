Author: Marta Gonzalvo-Ulla, Caltech

Scripts to build RNA structures with modified and unnatural residues (2'-O-methyl, locked nucleic acids, 2’-fluororibose, 2’-O-methoxyethyl, 5’-methylcytosine, phosphorothioate), run molecular dynamics simulations, and analyze the impact of different sequences on structure and stability.

The sample siRNA is comprised of one long core strand, and two shorter sensor and guide strands that both base pair with core strand, as reported by Han et. al (1). They form 2 double helices. Image included below. 

Protocol:

1. 2 A-form DNA double helices are created using PyMOL's (2) builder tool. Sample structures: helix1.pdb, helix2.pdb

2. The two helices are manually joined in single file and pulled close together with the right alignment for merging them. The linker C3 is also added manually. Then, the chains are manually reordered in pdb to have each A,B,C chain sequentially. The RL residue 2 is then renamed to become part of the adjacent residue (-22). Intermediate structure: mergedstructures.pdb. Sample final structure: duplex.pdb

3. The modifications are substituted into the structures and folders are created for each simulation with subst_prep.py based on parameters in json file. Sample json file: simulations.json. Sample final folder: foldersimulation.

    This involves a few major steps:

    3a) Substituting monomers in place with substnucl_wholefile.py (depends on functions in substnucl.py, needs file describing modifications+RNA residues, sample: substitute.txt, monomers in monomer_struct)

    3b) Parameterize, create topol.top and other .itp files using monomers from modXNA library and custom fragments parameterized with REDS using parameterize_gmx_modxna.py (hardcoded for monomers in this work). A detailed note on parameterization is included in parameterization.txt (11-19).

    3c) Create box, run solvation, add ions using Gromacs (3-10).

    3d) Creating .mdp files with correct and desired constraints using biashbond_rna.py, (customized hard-coded file)

    -- Command: python subst_prep.py samplefiles/simulations.json samplefiles (script is also hard-coded)

4. Equilibration and production molecular dynamics are run using Gromacs. Sample script: runsimulation.sh. 

    - Will generate xtc, other standard gromacs files. Then pdb of non-water atoms is generated from xtc with a command similar to: 
        gmx trjconv -f md.xtc -s md.gro -o md.pdb -pbc nojump
        Sample pdb: md.pdb

5. Results are analyzed from trajectory pdbs:

    5a) First, a csv is generated with the values of structural measures of interest at every timestep in the output pdbs using analysis.py. Outputs certain coordinates, distances and angles, hardcoded into the file. Sample output: summary_result1.csv

    --- Command: python analysis.py samplefiles/results samplefiles/templatefolder

    5b) The csv results are visualized for one or more simulations using visualiz_analysis.py. Types of files is hardcoded. Sample plot 5 simulations: rmsd5simulations.png

    --- Command: python visualiz_analysis.py samplefiles/results


This github repository is a work in progress and accompanies a publication in preparation. No  The project has been funded by Switch Therapeutics (https://www.switchthera.com).

Software used: pandas, MDTraj, MDAnalysis, Matplotlib (21-25)

(1-25): All references can be found in references.txt file.

![3-strand double helix construct image](https://github.com/martagu/rna_modify_MD_analyze/blob/main/samplefiles/md1pdb.png?raw=true) XXXX REDO WITH CORRECT T2/T1 STRUCTURE
