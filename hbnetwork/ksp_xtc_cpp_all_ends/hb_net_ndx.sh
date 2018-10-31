###### Let us set the initial input as a pdb file
echo "#####################################################################################"
echo "#### !!!! Please be careful of version of gromacs with MPI ##########################"
echo "#### If you find error with tmp.ndx please check the script file to version of gromacs"
echo "#### gmx, gmx_d, gmx_mpi, gmx_mpi_d			     ######################"
echo "#####################################################################################"
echo "#### Example to generate *.pdb file     "
echo "#gmx trjconv -f traj.xtc -s topol.tpr -o sample.pdb "
echo "#### Example to generate residue file   "
echo "#for i in {1..432}; do echo -n $i" "; done > res.txt"
echo "#####################################################################################"

echo "Please input the name of the pdb file generated from Gromacs 'trjconv' funtion"
read pdb_file_in

if [ ${pdb_file_in##*.} != "pdb" ]
then
	pdb_file=$pdb_file_in".pdb"
else
	pdb_file=$pdb_file_in
fi

echo "pdb file to be used file index file generation >>> " $pdb_file

#echo "Please input the order of residue number in column for the pdb file starting from 1 (ATOM)"
#read res_col_ndx
#echo "Index of residue number in the pdb file is set to 5" $res_col_ndx
#echo "Warning !!! Index of residue number in the pdb file is set to 5"

echo "Please input the File name to include all residue indices (for example "res.txt")"

read res_file

echo "File including the indicies of the residues " $res_file

echo "Please input \"yes\" for MPI version of groamcs, otherwise input \"no\" to generate default index file from Gromacs"

read mpi_option

if [ ${mpi_option} == "yes" ]
then
	echo "Version of gromacs installed in your computer is with MPI : gmx_mpi (_d)"
else
	echo "Version of gromacs installed in your computer is without MPI : gmx"
fi

echo "Output filename for the index file"
read out_file_in
if [ ${out_file_in##*.} != "ndx" ]
then
        out_file=$out_file_in".ndx"
else
        out_file=$out_file_in
fi
echo "Output file >>> " $out_file

#Initialization of index data from the "*.pdb" file
if [ -f file_D_H.txt ]
then
	rm file_D_H.txt
fi

if [ -f file_A.txt ]
then
	rm file_A.txt
fi

if [ -f file_C.txt ]
then
	rm file_C.txt
fi

if [ -f file_ALL.txt ]
then
	rm file_ALL.txt
fi


cat $res_file | awk '{for(i=1;i<=NF;++i)printf "%d\n",$i}' > res_file_col

res_array=()
res_ndx=-1
while read -r res
do
	echo "input residue "$res
	
	#cat $pdb_file |	awk -v rr1=$res  -v col=$res_col_ndx '{if($1=="ATOM" && $col==rr1) print }' | awk 'col ~ /^[0-9]+$/ {print}' > tmp_jkh
	cat $pdb_file |	awk -v rr1=$res '{if(NF==12){if($1== "ATOM" && $4 != "SOL" && $6==rr1) print } else if($1== "ATOM" && $4 != "SOL" && NF==11) {if($5==rr1) print }}' > tmp_jkh

	chk_o_n="false"
	chk_h="false"

	num_h=0
	ndx_o_n_old=0
	ndx_o_n_vec_ndx=-1
	ndx_o_n_vec=()
	name_o_n_vec=()
	num_h_vec=()
	ndx_h_vec_ndx=-1
	ndx_h_vec=()
	ndx_C_vec_ndx=-1
	ndx_C_vec=()
	ndx_ALL_vec_ndx=-1
	ndx_ALL_vec=()

	while read -r line
	do
		name="$line"
		echo "$name" > tmp_line_jkh
		name_atom=`awk '{print $NF}' tmp_line_jkh`
		ndx_atom=`awk '{print $2}' tmp_line_jkh`
		#echo $name_atom" "$ndx_atom
		ndx_ALL_vec_ndx=$(($ndx_ALL_vec_ndx+1))
		ndx_ALL_vec[$ndx_ALL_vec_ndx]=$ndx_atom
		if [ "$name_atom" == "C" ]
		then

			#saving the atom number of C
			ndx_C_vec_ndx=$(($ndx_C_vec_ndx+1))
			ndx_C_vec[$ndx_C_vec_ndx]=$ndx_atom
		fi
		if [ "$name_atom" == "N" ] || [ "$name_atom" == "O" ] || [ "$name_atom" == "S" ]
		then
			chk_o_n="true"
			
			#saving the atom number of D or A
			ndx_o_n_vec_ndx=$(($ndx_o_n_vec_ndx+1))
			ndx_o_n_vec[$ndx_o_n_vec_ndx]=$ndx_atom
			name_o_n_vec[$ndx_o_n_vec_ndx]=$name_atom
			#echo ">>>D or A >> " $ndx_o_n_vec_ndx " ::: name " ${name_o_n_vec[$ndx_o_n_vec_ndx]} " ::: atom_num >> " ${ndx_o_n_vec[$ndx_o_n_vec_ndx]}
			
			#counting the total number of Hydrogen atoms attached to each D or A
			num_h_vec[ndx_o_n_vec_ndx]=0
		fi
		#echo $ndx_o_n" "$chk_o_n
		if [ "$chk_o_n" == "true" ]
		then
			if [ "$name_atom" == "H" ]
			then
				#Increasing number of hydrogen atoms attached to each D or A
				num_h_vec[$ndx_o_n_vec_ndx]=$((${num_h_vec[$ndx_o_n_vec_ndx]}+1))
				#echo "number of hydrogens attached  >> " ${num_h_vec[$ndx_o_n_vec_ndx]}
				
				#Indexing hydrogen atoms in serial way while continuously counting
				ndx_h_vec_ndx=$(($ndx_h_vec_ndx+1))
				ndx_h_vec[$ndx_h_vec_ndx]=$ndx_atom
				#echo ">>>H >> " $ndx_h_vec_ndx " :: atom_num >> " ${ndx_h_vec[$ndx_h_vec_ndx]}
			fi
			if [ "$name_atom" == "C" ]
			then
				chk_o_n="false"
			fi
		fi
		ndx_o_n_old=$ndx_o_n
	done < tmp_jkh

	#echo "Total number of donor and acceptor groups ---  " $(($ndx_o_n_vec_ndx+1))
	#echo "Total number of hydrogens attached to     ---  " $(($ndx_h_vec_ndx+1))

	echo "Printing index file for H-bond network analysis......"
	h_out_ndx=0
	for ((i=0; i<=$ndx_o_n_vec_ndx; i++))
	do
		#echo $i" "${name_o_n_vec[$i]}" "${ndx_o_n_vec[$i]}
		#echo "# of H :"${num_h_vec[$i]}
		if [[ "${num_h_vec[$i]}" -eq 0 ]]
		then
			echo ${ndx_o_n_vec[$i]} >> file_A.txt
		elif [[ "${num_h_vec[$i]}" -eq 1 && "${name_h_vec[$i]}" -eq "N" ]]
		then
			echo ${ndx_o_n_vec[$i]}" "${ndx_h_vec[$h_out_ndx]} >> file_D_H.txt
			h_out_ndx=$(($h_out_ndx+1))
		elif [[ "${num_h_vec[$i]}" -eq 1 && "${name_h_vec[$i]}" -eq "O" ]]
		then
			echo ${ndx_o_n_vec[$i]} >> file_A.txt
			echo ${ndx_o_n_vec[$i]}" "${ndx_h_vec[$h_out_ndx]} >> file_D_H.txt
			h_out_ndx=$(($h_out_ndx+1))
		else
			for ((j=0; j<${num_h_vec[$i]}; j++))
			do
				echo  ${ndx_o_n_vec[$i]}" "${ndx_h_vec[$h_out_ndx]} >> file_D_H.txt
				h_out_ndx=$(($h_out_ndx+1))
			done	
		fi
	done

	echo "Printing index file for Carbon atoms for hydration shell anaylsis......"
	for ((i=0; i<=$ndx_C_vec_ndx; i++))
	do
		echo ${ndx_C_vec[$i]} >> file_C.txt
	done
	for ((i=0; i<=$ndx_ALL_vec_ndx; i++))
	do
		echo ${ndx_ALL_vec[$i]} >> file_ALL.txt
	done
done < res_file_col

if [ ${mpi_option} == "yes" ]
then
	echo "q" | gmx_mpi make_ndx -f $pdb_file -o tmp.ndx
else
	echo "q" | gmx make_ndx -f $pdb_file -o tmp.ndx
fi
#echo "q" | gmx_mpi_d make_ndx -f sample.pdb -o tmp.ndx
#echo "q" | mpirun -np 1 gmx_mpi make_ndx -f sample.pdb -o tmp.ndx
#echo "q" | mpirun -np 1 gmx_mpi_d make_ndx -f sample.pdb -o tmp.ndx

echo "[ N-HBOND ]" > $out_file
awk '{if(NR%16==15){printf"%d\n", $1}else{printf"%d ", $1}}' file_D_H.txt >> $out_file
echo " " >> $out_file
echo "[ H-HBOND ]" >> $out_file 
awk '{if(NR%16==15){printf"%d\n", $2}else{printf"%d ", $2}}' file_D_H.txt >> $out_file
echo " " >> $out_file
echo "[ O-HBOND ]" >> $out_file
awk '{if(NR%16==15){printf"%d\n", $1}else{printf"%d ", $1}}' file_A.txt >> $out_file
echo " " >> $out_file
echo "[ HEAVY_ACTIVE_GROUP ]" >> $out_file
awk '{if(NR%16==15){printf"%d\n", $1}else{printf"%d ", $1}}' file_C.txt >> $out_file
echo " " >> $out_file
echo "[ ALL_ACTIVE_GROUP ]" >> $out_file
awk '{if(NR%16==15){printf"%d\n", $1}else{printf"%d ", $1}}' file_ALL.txt >> $out_file
echo " " >> $out_file
cat tmp.ndx >> $out_file

#rm tmp.ndx
#rm tmp*
rm *#

