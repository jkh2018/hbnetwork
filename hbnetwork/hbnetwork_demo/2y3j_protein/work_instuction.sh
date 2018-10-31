1) source /scratch/kjong/.hbnet/bin/GMXRC

2) gmx trjconv -f 2y3j.xtc -s 2y3j.tpr -o sample.pdb -dump 0 

3) echo "1 2 3 4 5 6" > 2y3j_res.txt

4) /scratch/kjong/hbnetwork/ksp_xtc_cpp_all_ends/hb_net_ndx.sh 

#add following info

sample.pdb

2y3j_res.txt

no

2y3j.ndx

5) /scratch/kjong/hbnetwork/ksp_xtc_cpp_all_ends/bin/hb_network_hydration.x < input_2y3j
