Author: Marta Gonzalvo

Project goal: understand siRNA structures through simulation

Approach: Build rna structures with modified and unnatural residues, run molecular dynamics simulations and analyze impact of different sequences on structure and stability


siRNA of interest is constituted of one long core strand and two shorter sensor and guide strands that both base pair with core strand. They form 2 double helices.

Protocol:

1. 2 A-form DNA double helices are created using PyMOL's builder tool. Other options including ROSETTA were explored, PyMOL was found to work best to start helices in regular conformation and stable in MD since later need to join them. Sample structures: helix1.pdb, helix2.pdb

2. The two helices are manually joined in single file and pulled close together with the right alignment for merger in later steps. Linker C3 is also added manually (RL residue 2 manually changed name to be part of previous residue). Chains are manually reordered in pdb to have each A,B,C chain sequentially.  Intermediate structure: mergedstructures.pdb. Sample structure: duplex.pdb


3. Structures are modified and gromacs folders are created and run with subst_prep.py based on parameters in json file. Sample json file: simulations.json. Sample final folder: foldersimulation.

This involves a few major steps:

- Substituting monomers in place with substnucl_wholefile.py (depends on functions in substnucl.py, needs file describing modifications+RNA residues, sample: substitute.txt, monomers in monomer_struct)
- Run solvation, adding ions using Gromacs: Using modified force field for modified and non-standard residues*
- Creating .mdp files with correct and desired constraints (depends on biashbond_rna.py, very customized hard-coded file)

-- Command: python subst_prep.py samplefiles/simulations.json samplefiles (script is also hard-coded)

4. Equilibration and production molecular dynamics are run using gromacs. Sample script: runsimulation.sh in foldersimulation. 

- Will generate xtc, other standard gromacs files. Then pdb of non-water atoms is generated from xtc with a command similar to: 
    gmx trjconv -f md.xtc -s md.gro -o md.pdb -pbc nojump
    Sample pdb: md.pdb

5. Results are analyzed from trajectory pdbs:

5a) First, a csv is generated with the values of structural measures of interest at every timestep in the output pdbs using analysis.py. Outputs certain coordinates, distances and angles, hardcoded into file. Sample output: summary_result1.csv
--- Command: python analysis.py samplefiles/results samplefiles/templatefolder

5b) The csv results are visualized for one or more simulations using visualiz_analysis.py. Types of files is hardcoded. Sample plot 5 simulations: rmsd5simulations.png
--- Command: python visualiz_analysis.py samplefiles/results


This github repository accompanies a publication in preparation.



*The force field (FF) AMBER14SB was used (1). Parameters for the modified nucleotides were added to make a modified FF from the literature: 2'OMe (2), LNA (3) and PS (4). For the unreported LNA thymidine, the parameters were derived from (3) for the sugar ring and from AMBER03FF for the base (5). Its charges were calculated with the RESP ESP charge Derive (RED) server (6). The parameters for the C3 linker are from the GAFF FF (7). 

References:

1. Maier, J. A. et al. ff14SB: improving the accuracy of protein side chain and backbone parameters from ff99SB. J. Chem. Theory Comput. 11, 3696–3713 (2015) 
2. Aduri, R., Psciuk, B.T., Saro, P., Taniga, H., Schlegel, H.B., and SantaLucia, J. (2007). AMBER force field parameters for the naturally occurring modified nucleosides in RNA. J. Chem. Theor. Comput. 3, 1464–1475.
3. Condon, D.E., Yildirim, I., Kennedy, S.D., Mort, B.C., Kierzek, R., and Turner, D.H. (2014). Optimization of an AMBER force field for the artificial nucleic acid, LNA, and benchmarking with NMR of L(CAAU). J. Phys. Chem. B 118, 1216–1228.
4. Lind, K.E., Sherlin, L.D., Mohan, V., Griffey, R.H., and Ferguson, D.M. (1997). Parameterization and simulation of the physical properties of phosphorothioate nucleic acids. In Molecular Modeling of Nucleic Acids, 682 (American Chemical Society), pp. 41–54.
5. Duan, Y., Wu, C., Chowdhury, S., Lee, M.C., Xiong, G., Zhang, W., Yang, R., Cieplak, P., Luo, R., Lee, T., et al. (2003). A point-charge force field for molecular mechanics simulations of proteins based on condensed-phase quantum mechanical calculations. J. Comput. Chem. 24, 1999–2012.
6.(http://q4md- forcefieldtools.org/REDServer/).
7. Wang, J., Wolf, R.M., Caldwell, J.W., Kollman, P.A., and Case, D.A. (2004). Development and testing of a general amber force field. J. Comput. Chem. 25, 1157–1174.

![3-strand double helix construct image](https://github.com/martagu/rna_modify_MD_analyze/blob/main/samplefiles/md1pdb.png?raw=true)
