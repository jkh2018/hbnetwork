1) source /scratch/kjong/.hbnet/bin/GMXRC

2) gmx trjconv -f tip4p2005.xtc -s tip4p2005.tpr -o sample.pdb -dump 2000 

3) echo "0" > tip4p2005_res.txt 

4) /scratch/kjong/hbnetwork/ksp_xtc_cpp_all_ends/hb_net_ndx.sh 

#add following info

sample.pdb

tip4p2005_res.txt

no

tip4p2005.ndx

5) /scratch/kjong/hbnetwork/ksp_xtc_cpp_all_ends/bin/hb_network_hydration.x < input_tip4p2005
