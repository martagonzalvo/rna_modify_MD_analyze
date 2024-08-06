Author: Marta Gonzalvo-Ulla, Caltech

Scripts to build RNA structures with modified and unnatural residues (2'-O-methyl [2'OMe], locked nucleic acids [LNA], phosphorothioate [PS]), run molecular dynamics simulations, and analyze the impact of different sequences on structure and stability.

The sample siRNA is comprised of one long core strand, and two shorter sensor and guide strands that both base pair with core strand. They form 2 double helices. Image included below.

Protocol:

1. 2 A-form DNA double helices are created using PyMOL's* builder tool. Sample structures: helix1.pdb, helix2.pdb

2. The two helices are manually joined in single file and pulled close together with the right alignment for merging them. The linker C3 is also added manually. Then, the chains are manually reordered in pdb to have each A,B,C chain sequentially. The RL residue 2 is then renamed to become part of the adjacent residue (-22). Intermediate structure: mergedstructures.pdb. Sample final structure: duplex.pdb


3. The modifications are substituted into the structures and folders are created for each simulation with subst_prep.py based on parameters in json file. Sample json file: simulations.json. Sample final folder: foldersimulation.

        This involves a few major steps:

            - Substituting monomers in place with substnucl_wholefile.py (depends on functions in substnucl.py, needs file describing modifications+RNA residues, sample: substitute.txt, monomers in monomer_struct)
            - Run solvation, adding ions using Gromacs**: Using modified force field for modified and non-standard residues***
            - Creating .mdp files with correct and desired constraints (depends on biashbond_rna.py, customized hard-coded file)

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


This github repository accompanies a publication in preparation.

Software used: pandas, MDTraj, MDAnalysis, Matplotlib****

*PyMOL Reference: The PyMOL Molecular Graphics System, Version 2.5.8 Schrödinger, LLC.

** Gromacs References:
1. H. Bekker, H.J.C. Berendsen, E.J. Dijkstra, S. Achterop, R. van Drunen, D. van der Spoel, A. Sijbers, and H. Keegstra et al., “Gromacs: A parallel computer for molecular dynamics simulations”; pp. 252–256 in Physics computing 92. Edited by R.A. de Groot and J. Nadrchal. World Scientific, Singapore, 1993.
2. H.J.C. Berendsen, D. van der Spoel, and R. van Drunen, “GROMACS: A message-passing parallel molecular dynamics implementation,” Comp. Phys. Comm., 91 43–56 (1995).
3. E. Lindahl, B. Hess, and D. van der Spoel, “GROMACS 3.0: A package for molecular simulation and trajectory analysis,” J. Mol. Mod., 7 306–317 (2001).
4. D. van der Spoel, E. Lindahl, B. Hess, G. Groenhof, A.E. Mark, and H.J.C. Berendsen, “GROMACS: Fast, Flexible and Free,” J. Comp. Chem., 26 1701–1718 (2005).
5. B. Hess, C. Kutzner, D. van der Spoel, and E. Lindahl, “GROMACS 4: Algorithms for Highly Efficient, Load-Balanced, and Scalable Molecular Simulation,” J. Chem. Theory Comput., 4 [3] 435–447 (2008).
6. S. Pronk, S. Páll, R. Schulz, P. Larsson, P. Bjelkmar, R. Apostolov, M.R. Shirts, and J.C. Smith et al., “GROMACS 4.5: A high-throughput and highly parallel open source molecular simulation toolkit,” Bioinformatics, 29 [7] 845–854 (2013).
7. S. Páll, M.J. Abraham, C. Kutzner, B. Hess, and E. Lindahl, “Tackling exascale software challenges in molecular dynamics simulations with GROMACS”; pp. 3–27 in Solving software challenges for exascale. Edited by S. Markidis and E. Laure. Springer International Publishing Switzerland, London, 2015.
8. M.J. Abraham, T. Murtola, R. Schulz, S. Páll, J.C. Smith, B. Hess, and E. Lindahl, “GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers,” SoftwareX, 1–2 19–25 (2015).

***The force field (FF) AMBER14SB was used (1). Parameters for the modified nucleotides were added to make a modified FF from the literature: 2'OMe (2), LNA (3) and PS (4). For the unreported LNA thymidine, the parameters were derived from (3) for the sugar ring and from AMBER03FF for the base (5). Its charges were calculated with the RESP ESP charge Derive (RED) server (6). The parameters for the C3 linker are from the GAFF FF (7). FF References:

1. Maier, J. A. et al. ff14SB: improving the accuracy of protein side chain and backbone parameters from ff99SB. J. Chem. Theory Comput. 11, 3696–3713 (2015) 
2. Aduri, R., Psciuk, B.T., Saro, P., Taniga, H., Schlegel, H.B., and SantaLucia, J. (2007). AMBER force field parameters for the naturally occurring modified nucleosides in RNA. J. Chem. Theor. Comput. 3, 1464–1475.
3. Condon, D.E., Yildirim, I., Kennedy, S.D., Mort, B.C., Kierzek, R., and Turner, D.H. (2014). Optimization of an AMBER force field for the artificial nucleic acid, LNA, and benchmarking with NMR of L(CAAU). J. Phys. Chem. B 118, 1216–1228.
4. Lind, K.E., Sherlin, L.D., Mohan, V., Griffey, R.H., and Ferguson, D.M. (1997). Parameterization and simulation of the physical properties of phosphorothioate nucleic acids. In Molecular Modeling of Nucleic Acids, 682 (American Chemical Society), pp. 41–54.
5. Duan, Y., Wu, C., Chowdhury, S., Lee, M.C., Xiong, G., Zhang, W., Yang, R., Cieplak, P., Luo, R., Lee, T., et al. (2003). A point-charge force field for molecular mechanics simulations of proteins based on condensed-phase quantum mechanical calculations. J. Comput. Chem. 24, 1999–2012.
6.(http://q4md- forcefieldtools.org/REDServer/).
7. Wang, J., Wolf, R.M., Caldwell, J.W., Kollman, P.A., and Case, D.A. (2004). Development and testing of a general amber force field. J. Comput. Chem. 25, 1157–1174.

****Citations software:
1. The pandas development team. pandas-dev/pandas: Pandas (2020) Zenodo, https://doi.org/10.5281/zenodo.3509134
2. McKinney, W. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference (2010), 56-61
3. McGibbon, R. T., Beauchamp, K. A., Harrigan, M. P., Klein, C., Swails, J. M., Hernández, C. X., Schwantes, C. R., Wang, L.-P., Lane, T. J., and Pande, V. S. (2015). MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories. Biophysical Journal. 109, 8, 1528-1532.
4. R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler, D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein. MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations. In S. Benthall and S. Rostrup, editors, Proceedings of the 15th Python in Science Conference, pages 98-105, Austin, TX, 2016. SciPy, doi:10.25080/majora-629e541a-00e.
5. N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319-2327, doi:10.1002/jcc.21787. PMCID:PMC3144279
6. J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.



![3-strand double helix construct image](https://github.com/martagu/rna_modify_MD_analyze/blob/main/samplefiles/md1pdb.png?raw=true)
