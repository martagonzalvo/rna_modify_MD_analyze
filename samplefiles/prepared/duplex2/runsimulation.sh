gmx grompp -f em_ct.mdp -c dup_ions.gro -p topol.top -o em.tpr -n index.ndx -maxwarn 2
gmx mdrun -deffnm em -v

gmx grompp -f nvt_ct.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx -maxwarn 1
gmx mdrun -deffnm nvt -v

gmx grompp -f md.mdp -c nvt.gro -r nvt.gro -p topol.top -o md.tpr -n index.ndx -maxwarn 4
gmx mdrun -deffnm md -v
