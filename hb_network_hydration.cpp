#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <map>
#include <queue>
#include <vector>
#include <algorithm>
#include <time.h>
#include <math.h>
#include "Graph.h"
#include <omp.h>

#define INF 9999999

//using namespace std;
using namespace Graph_all;

//-----defintion of vector type--------------------//
typedef vector<double> double_vec_t;
typedef vector<int> int_vec_t;


class HBNetwork
{

private:
void read_input(FILE* fp, 
	int& first_frame, 
	int& last_frame, 
	string& xtcfile, 
	string& tprfile, 
	string& ndxfile, 
	string& outfile, 
	int& k_max, 
	string& direct_id, 
	string& weight_id, 
	string& ring_id, 
	string& print_path_id, 
	string& shell_id, 
	double& pop_shell, 
	string& scramble_id, 
	int& last_pro_res_ndx, 
	float& hb_dist_cut, 
	float& hb_ang_cut, 
	float& hb_weight_dist_cut, 
	double& dist_DA_D0, 
	double& dist_DA_R0, 
	double& dist_HA_D0, 
	double& dist_HA_R0, 
	double& ang_D0, 
	double& ang_R0, 
	int& max_lag, 
	double& r_shell_1_N_O, 
	double& r_shell_1_C, 
	double& r_shell_bulk, 
	double& scrambling_ti, 
	double& scrambling_tf, 
	int& WDIM, 
	string& n_hb_grp, 
	string& h_hb_grp, 
	string& o_hb_grp, 
	string& ow_grp, 
	string& c_grp, 
	string& active_grp, 
	string& print_adjmatrix_id, 
	string& print_defect_id, 
	string& fpt_id, 
	string& fpt_interface_grp, 
	string& fpt_jump_id, 
	double& fpt_box_epsilon, 
	double& fpt_box_width, 
	float& fpt_time_width)
{
	first_frame=0;
	last_frame=1;
	xtcfile="traj.xtc";
	tprfile="topol.tpr";
	ndxfile="index.ndx";
	outfile="test";
	k_max=1;
	direct_id="directed";
	weight_id="weight";
	ring_id="DA";
	print_path_id="no";
	shell_id="no";
	//pop_shell=0.3;
	scramble_id="no";
	last_pro_res_ndx = 6;
	//hb_dist_cut=0.35;
	//hb_ang_cut=double(30)/180*M_PI;
	//hb_weight_dist_cut = 0.45;
	//dist_DA_D0 = 0.0;
	//dist_DA_R0 = 0.35;
	//dist_HA_D0 = 0.0;
	//dist_HA_R0 = 0.25;
	//ang_D0 = 0; //0 degrees
	//ang_R0 = double(30)/180*M_PI; //30 degrees
	max_lag = 5000;
	r_shell_1_N_O = 0.35;
	r_shell_1_C = 0.5;
	r_shell_bulk = 1.5;
	scrambling_ti = 0; //ps
	scrambling_tf = 1000; //ps
	WDIM = 3;
	print_adjmatrix_id = "no";
	print_defect_id = "no";
  
	n_hb_grp = "N-HBOND";
	h_hb_grp = "H-HBOND";
	o_hb_grp = "O-HBOND";
	ow_grp   = "SOL";
        c_grp    = "HEAVY_ACTIVE_GROUP";
	active_grp = "ALL_ACTIVE_GROUP";
	fpt_id = "no";
	fpt_interface_grp = "Protein";
	fpt_jump_id = "x";
	fpt_box_epsilon = 0.1;
	fpt_box_width = 0.5;
	fpt_time_width = 100;
	
	string line;
  	line.resize(1000);
  	char buffer[1000];
  	char buffer1[1000];

  	while(fgets(buffer,256,fp)){
   		line=buffer;
    		for(int i=0;i<line.length();++i) if(line[i]=='#' || line[i]=='\n') line.erase(i);
    		for(int i=line.length()-1;i>=0;--i){
      			if(line[i]!=' ')break;
      			line.erase(i);
    		}
		buffer[0]=0;
    		sscanf(line.c_str(),"%s",buffer);
	    	if(strlen(buffer)==0) continue;
	    	string keyword=buffer;
	    	if(keyword=="first_frame")
		{
	      		sscanf(line.c_str(),"%s %d",buffer,&first_frame);
		}
		else if(keyword=="last_frame")
		{
	      		sscanf(line.c_str(),"%s %d",buffer,&last_frame);
		}
	    	else if(keyword=="xtcfile")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			xtcfile=buffer1;
	    	}
	    	else if(keyword=="tprfile")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			tprfile=buffer1;
	    	}
	    	else if(keyword=="ndxfile")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			ndxfile=buffer1;
	    	}
	    	else if(keyword=="outfile")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			outfile=buffer1;
	    	}
	    	else if(keyword=="k_max")
		{
	      		sscanf(line.c_str(),"%s %d",buffer,&k_max);
		}
	    	else if(keyword=="hbond_dist")
		{
	      		sscanf(line.c_str(),"%s %f",buffer,&hb_dist_cut);
		}
	    	else if(keyword=="hbond_ang")
		{
	      		sscanf(line.c_str(),"%s %f",buffer,&hb_ang_cut);
			hb_ang_cut = hb_ang_cut/180*M_PI;
		}
	    	else if(keyword=="direct_id")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			direct_id=buffer1;
	    	}
	    	else if(keyword=="weight_id")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			weight_id=buffer1;
	    	}
	    	else if(keyword=="hbond_weight_dist")
		{
	      		sscanf(line.c_str(),"%s %f", buffer, &hb_weight_dist_cut);
		}
	    	else if(keyword=="ring_id")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			ring_id=buffer1;
	    	}
	    	else if(keyword=="print_path")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
	      		if(buffer1[0]=='T' || buffer1[0]=='t') print_path_id="yes";
	    	}
	    	else if(keyword=="print_adjmatrix")
	    	{
			sscanf(line.c_str(),"%s %s",buffer,buffer1);
                        if(buffer1[0]=='T' || buffer1[0]=='t') print_adjmatrix_id="yes";
	    	}
	    	else if(keyword=="print_defect")
	    	{
			sscanf(line.c_str(),"%s %s",buffer,buffer1);
                        if(buffer1[0]=='T' || buffer1[0]=='t') print_defect_id="yes";
	    	}
	    	else if(keyword=="shell_analysis")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
	      		if(buffer1[0]=='T' || buffer1[0]=='t') shell_id="yes";
	    	}
	    	else if(keyword=="pop_shell")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&pop_shell);
		}
	    	else if(keyword=="scrambling_analysis")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
	      		if(buffer1[0]=='T' || buffer1[0]=='t') scramble_id="yes";
	    	}
	    	else if(keyword=="last_protein_residue")
		{
	      		sscanf(line.c_str(),"%s %d",buffer,&last_pro_res_ndx);
		}
	    	else if(keyword=="dist_DA_D0")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&dist_DA_D0);
		}
	    	else if(keyword=="dist_DA_R0")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&dist_DA_R0);
		}
	    	else if(keyword=="dist_HA_D0")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&dist_HA_D0);
		}
	    	else if(keyword=="dist_HA_R0")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&dist_HA_R0);
		}
	    	else if(keyword=="ang_D0")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&ang_D0);
		}
	    	else if(keyword=="ang_R0")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&ang_R0);
			ang_R0 = ang_R0/180*M_PI;
		}
	    	else if(keyword=="max_lag")
		{
	      		sscanf(line.c_str(),"%s %d",buffer,&max_lag);
		}
	    	else if(keyword=="r_shell_1_N_O")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&r_shell_1_N_O);
		}
	    	else if(keyword=="r_shell_1_C")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&r_shell_1_C);
		}
	    	else if(keyword=="r_shell_bulk")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&r_shell_bulk);
		}
	    	else if(keyword=="scrambling_ti")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&scrambling_ti);
		}
	    	else if(keyword=="scrambling_tf")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&scrambling_tf);
		}
	    	else if(keyword=="WDIM")
		{
	      		sscanf(line.c_str(),"%s %d",buffer,&WDIM);
		}
	    	else if(keyword=="Hbond_donor_grp")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			n_hb_grp=buffer1;
	    	}
	    	else if(keyword=="Hbond_H_grp")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			h_hb_grp=buffer1;
	    	}
	    	else if(keyword=="Hbond_acceptor_grp")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			o_hb_grp=buffer1;
	    	}
	    	else if(keyword=="water_grp")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			ow_grp=buffer1;
	    	}
	    	else if(keyword=="active_heavy_grp")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			c_grp=buffer1;
	    	}
	    	else if(keyword=="active_all_grp")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			active_grp=buffer1;
	    	}
	    	else if(keyword=="fpt_analysis")
	    	{
			sscanf(line.c_str(),"%s %s",buffer,buffer1);
                        if(buffer1[0]=='T' || buffer1[0]=='t') fpt_id="yes";
	    	}
	    	else if(keyword=="fpt_interface_grp")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			fpt_interface_grp=buffer1;
	    	}
	    	else if(keyword=="fpt_jump_direction")
	    	{
	      		sscanf(line.c_str(),"%s %s",buffer,buffer1);
			fpt_jump_id=buffer1;
	    	}
	    	else if(keyword=="fpt_box_epsilon")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&fpt_box_epsilon);
		}
	    	else if(keyword=="fpt_box_width")
		{
	      		sscanf(line.c_str(),"%s %lf",buffer,&fpt_box_width);
		}

	    	else if(keyword=="fpt_time_width")
		{
	      		sscanf(line.c_str(),"%s %f",buffer,&fpt_time_width);
		}
	    	else{
	      		fprintf(stderr,"Unknown keywords :%s\n",keyword.c_str());
	      		exit(1);
	    	}
	  }

}

void adjcent_matrix_calc(string direct_id, triclinicbox box, vector <coordinates> xyz_nNH, vector <coordinates> xyz_hNH, vector <coordinates> xyz_CO, vector <coordinates> xyz_OW, double** adj_matrix, int** connection_matrix, int* defect_in, int* defect_out, float hb_dist_cut, float hb_ang_cut, int WDIM)
{
	unsigned i, j, k, m;

	int n_CO=xyz_CO.size();
	int n_NH=xyz_nNH.size();
	int n_OW=xyz_OW.size()/WDIM;
	int N_matrix=n_CO+n_NH+n_OW;
	
	#pragma omp parallel private(i, j)
	{
		#pragma omp for 
		for(i =0; i<N_matrix; ++i)
		{
			for(j =0; j<N_matrix; ++j){
				adj_matrix[i][j] = INF;
			}
		}
		#pragma omp for
		for(i =0; i<N_matrix; ++i)
		{
			for(j =0; j<N_matrix; ++j){
				connection_matrix[i][j] = 0;
			}
		}
	}
				
      	double hb_dist, hb_ang;
	
      	/*Construction of Adjcent Matrix*/
	#pragma omp parallel private(j, k, hb_dist, hb_ang)
	{
		#pragma omp for 
		for(j =0;j<n_CO;++j){   //Oxygens
			for(k =0;k<n_NH;k++){ // NH groups as Donor and CO of Protein as Acceptor
				hb_dist = distance(xyz_nNH.at(k),xyz_CO.at(j),box);
				if(hb_dist < hb_dist_cut){
					hb_ang = bond_angle(xyz_hNH.at(k), xyz_nNH.at(k), xyz_CO.at(j), box); //!(H,Donor,Acceptor)->(x,y,z)
					if(hb_ang < hb_ang_cut){
						adj_matrix[k][n_OW+n_NH+j] = 1;
						connection_matrix[k][n_OW+n_NH+j] = 1;
						if(direct_id.compare( "undirected") == 0){
							adj_matrix[n_OW+n_NH+j][k] = 1;
						}
					}
				}
			}
		}
	}
	#pragma omp parallel private(j, k, m, hb_dist, hb_ang)
	{
		#pragma omp for
		for(j =0;j<n_OW;++j){
			for(k =0;k<n_NH;k++){ 
				hb_dist = distance(xyz_nNH.at(k),xyz_OW.at(j*WDIM),box);
				if(hb_dist < hb_dist_cut){
					hb_ang = bond_angle(xyz_hNH.at(k), xyz_nNH.at(k), xyz_OW.at(j*WDIM), box); // !(H,Donor,Acceptor)->(x,y,z)
					if(hb_ang < hb_ang_cut) {
						adj_matrix[k][n_NH+j] = 1;
						connection_matrix[k][n_NH+j] = 1;
						if(direct_id.compare( "undirected") == 0){
							adj_matrix[n_NH+j][k] = 1;
						}
					}
				}
			}     
			//Hbonds inside water : first Oxygen is Donor and second one is Acceptor 
			for(k =(j+1);k<n_OW;k++){
				hb_dist = distance(xyz_OW.at(j*WDIM),xyz_OW.at(k*WDIM),box);
				if(hb_dist < hb_dist_cut) {
					//calcuation of hbonds from j to k
					for(m =1;m<3;m++){
						hb_ang = bond_angle(xyz_OW.at(j*WDIM+m), xyz_OW.at(j*WDIM), xyz_OW.at(k*WDIM), box); //!(H,Donor,Acceptor)->(x,y,z)
						if(hb_ang < hb_ang_cut) {
							adj_matrix[n_NH+j][n_NH+k] = 1;
							connection_matrix[n_NH+j][n_NH+k] = 1;
							if(direct_id.compare( "undirected") == 0){
								adj_matrix[n_NH+k][n_NH+j]  = 1;
							}
							break;
						}
					}	
					//Calculation of hbonds from k to j
					for(m =1;m<3;m++){
						hb_ang = bond_angle(xyz_OW.at(k*WDIM+m), xyz_OW.at(k*WDIM), xyz_OW.at(j*WDIM), box); //!(H,Donor,Acceptor)->(x,y,z)
						if(hb_ang < hb_ang_cut) {
							adj_matrix[n_NH+k][n_NH+j] = 1;
							connection_matrix[n_NH+k][n_NH+j] = 1;
							if(direct_id.compare( "undirected") == 0){
								adj_matrix[n_NH+j][n_NH+k]  = 1;
							}
							break;
						}
					}
				}
			}  
			//Hbonds between CO group and water (oxygen) => CO is only Acceptor, and O is Donor
			for(k =0;k<n_CO;k++){
				hb_dist = distance(xyz_OW.at(j*WDIM),xyz_CO.at(k),box);
				if(hb_dist < hb_dist_cut) {
					for(m =1;m<3;m++){
						hb_ang = bond_angle(xyz_OW.at(j*WDIM+m), xyz_OW.at(j*WDIM), xyz_CO.at(k), box); //!(H,Donor,Acceptor)->(x,y,z)
						if(hb_ang < hb_ang_cut) {
							adj_matrix[n_NH+j][n_OW+n_NH+k] = 1;
							connection_matrix[n_NH+j][n_OW+n_NH+k] = 1;
							if(direct_id.compare( "undirected") == 0){
								adj_matrix[n_OW+n_NH+k][n_NH+j]  = 1;
							}
							break;
						}
					}
				}
			}
		}
	}
	// #pragma omp parallel private(i, j)
	{
		// #pragma omp for
		for(i=0; i < N_matrix; ++i){
			defect_in[i] = 0;
			defect_out[i] = 0;
			for(j=0; j < N_matrix; ++j){
				defect_in[i] += connection_matrix[j][i];
				defect_out[i] += connection_matrix[i][j];
			}
		}
	}
	return;
}

void adjcent_matrix_weight_calc(string direct_id, triclinicbox box, vector <coordinates> xyz_nNH, vector <coordinates> xyz_hNH, vector <coordinates> xyz_CO, vector <coordinates> xyz_OW, double** adj_matrix, int** connection_matrix, int* defect_in, int* defect_out, float hb_weight_dist_cut, double dist_DA_D0, double dist_DA_R0, double dist_HA_D0, double dist_HA_R0, double ang_D0, double ang_R0, int WDIM)
{
	unsigned i, j, k, m;

	int n_CO=xyz_CO.size();
	int n_NH=xyz_nNH.size();
	int n_OW=xyz_OW.size()/WDIM;
	int N_matrix=n_CO+n_NH+n_OW;
	
	#pragma omp parallel private(i, j)
	{
		#pragma omp for
		for(i =0; i<N_matrix; ++i)
		{
			for(j =0; j<N_matrix; ++j){
				adj_matrix[i][j] = INF;
				connection_matrix[i][j] = 0;
			}
		}
	}
    	double hb_dist, hb_oh_dist, hb_ang;
	
    	/*Construction of Adjcent Matrix*/
	#pragma omp parallel private(j, k, hb_dist, hb_ang, hb_oh_dist)
	{
		#pragma omp for
		for(j =0;j<n_CO;++j){ 
			for(k =0;k<n_NH;k++){ 
				hb_dist = distance(xyz_nNH.at(k),xyz_CO.at(j),box);
				if(hb_dist < hb_weight_dist_cut){
					hb_ang = bond_angle(xyz_hNH.at(k), xyz_nNH.at(k), xyz_CO.at(j), box); //!(H,Donor,Acceptor)->(x,y,z)
					hb_oh_dist = distance(xyz_hNH.at(k), xyz_CO.at(j), box);
					adj_matrix[k][n_OW+n_NH+j] = inv_switch_func(hb_dist, dist_DA_D0, dist_DA_R0)*inv_switch_func(hb_oh_dist, dist_HA_D0, dist_HA_R0)*inv_switch_func(hb_ang, ang_D0, ang_R0);
				}
			}
		}
	}
	
	#pragma omp parallel private(j, k, m, hb_dist, hb_ang, hb_oh_dist)
	{
		#pragma omp for
		for(j =0;j<n_OW;++j){
			for(k =0;k<n_NH;k++){ 
				hb_dist = distance(xyz_nNH.at(k),xyz_OW.at(j*WDIM),box);
				if(hb_dist < hb_weight_dist_cut){
					hb_ang = bond_angle(xyz_hNH.at(k), xyz_nNH.at(k), xyz_OW.at(j*WDIM), box); // !(H,Donor,Acceptor)->(x,y,z)
					hb_oh_dist=distance(xyz_hNH.at(k), xyz_OW.at(j*WDIM), box);
					adj_matrix[k][n_NH+j] = inv_switch_func(hb_dist, dist_DA_D0, dist_DA_R0)*inv_switch_func(hb_oh_dist, dist_HA_D0, dist_HA_R0)*inv_switch_func(hb_ang, ang_D0, ang_R0);
				}
			}
			//Hbonds inside water : first Oxygen is Donor and second one is Acceptor 
			for(k =(j+1);k<n_OW;k++){
				hb_dist = distance(xyz_OW.at(j*WDIM),xyz_OW.at(k*WDIM),box);
				if(hb_dist < hb_weight_dist_cut){
					//calcuation of hbonds from j to k
					double tt = 0.;
					for(m =1;m<3;m++){
						hb_ang = bond_angle(xyz_OW.at(j*WDIM+m), xyz_OW.at(j*WDIM), xyz_OW.at(k*WDIM), box); //!(H,Donor,Acceptor)->(x,y,z)
						hb_oh_dist = distance(xyz_OW.at(j*WDIM+m), xyz_OW.at(k*WDIM), box);
						tt +=switch_func(hb_oh_dist, dist_HA_D0, dist_HA_R0)*switch_func(hb_ang, ang_D0, ang_R0);
					}
					adj_matrix[n_NH+j][n_NH+k] = inv_switch_func(hb_dist, dist_DA_D0, dist_DA_R0)/tt;
					//Calculation of hbonds from k to j
					tt=0.;
					for(m =1;m<3;m++){
						hb_ang = bond_angle(xyz_OW.at(k*WDIM+m), xyz_OW.at(k*WDIM), xyz_OW.at(j*WDIM), box); //!(H,Donor,Acceptor)->(x,y,z)
						hb_oh_dist = distance(xyz_OW.at(k*WDIM+m), xyz_OW.at(j*WDIM), box);
						tt += switch_func(hb_oh_dist, dist_HA_D0, dist_HA_R0)*switch_func(hb_ang, ang_D0, ang_R0);
					}
					adj_matrix[n_NH+k][n_NH+j] = inv_switch_func(hb_dist, dist_DA_D0, dist_DA_R0)/tt;
				}
			}
			//Hbonds between CO group and water (oxygen) => CO is only Acceptor, and O is Donor
			for(k =0;k<n_CO;k++){
				hb_dist = distance(xyz_OW.at(j*WDIM),xyz_CO.at(k),box);
				if(hb_dist < hb_weight_dist_cut){
					double tt=0.;
					for(m =1; m<3; m++){
						hb_ang = bond_angle(xyz_OW.at(j*WDIM+m), xyz_OW.at(j*WDIM), xyz_CO.at(k), box); //!(H,Donor,Acceptor)->(x,y,z)
						hb_oh_dist = distance(xyz_OW.at(j*WDIM+m), xyz_CO.at(k), box);
						tt += switch_func(hb_oh_dist, dist_HA_D0, dist_HA_R0)*switch_func(hb_ang, ang_D0, ang_R0);
					}
					adj_matrix[n_NH+j][n_OW+n_NH+k] = inv_switch_func(hb_dist, dist_DA_D0, dist_DA_R0)/tt;
				}
			}
		}
	}
	//Choosing the smallest weight direction as D-A direction between water molecules for the connection matrix
	#pragma omp parallel private(i, j)
	{
		#pragma omp for
		for(i=0; i < N_matrix; ++i){
			for(j =i+1; j < N_matrix; ++j){
				if( adj_matrix[i][j] > INF ) { adj_matrix[i][j] = INF; }
				if( adj_matrix[j][i] > INF ) { adj_matrix[j][i] = INF; }
				if( adj_matrix[i][j] < adj_matrix[j][i] )
				{ 
					connection_matrix[i][j] = 1;
				}else if(adj_matrix[i][j] > adj_matrix[j][i] ){
					connection_matrix[j][i] = 1;
				}else{
					connection_matrix[i][j] = 1;
					connection_matrix[j][i] = 1;
				}
			}
		}
	}
	//For the undirected adj_matrix, we made the symmetric matrix by assigning the edge weight (i,j) as minimum weight between i and j
	if(direct_id.compare( "undirected") == 0){
		#pragma omp parallel private(i, j)
		{
			#pragma omp for
			for(i=0; i < N_matrix; ++i){
				for(j =i+1; j < N_matrix; ++j){
					if( adj_matrix[i][j] < adj_matrix[j][i] )
					{ 
						adj_matrix[j][i] = adj_matrix[i][j];
					}else{
						adj_matrix[i][j] = adj_matrix[j][i];
					}
				}
			}
		}
	}
	// #pragma omp parallel private(i, j)
	{
		// #pragma omp for
		for(i=0; i < N_matrix; ++i){
			defect_in[i] = 0;
			defect_out[i] = 0;
			for(j=0; j < N_matrix; ++j){
				defect_in[i] += connection_matrix[j][i];
				defect_out[i] += connection_matrix[i][j];
			}
		}
	}
	return;
}

void floyd_calc(unsigned int N_floyd, double** adj_matrix_spanning, double** dist, int** through)
{	
	unsigned int m, n, k;
	//Initializing the distance and connection medium points
	#pragma omp parallel private(m, n)
	{
		#pragma omp for
		for (m = 0; m < N_floyd; ++m){
			for (n = 0; n < N_floyd; ++n){
				dist[m][n] =  adj_matrix_spanning[m][n];   //dist[m,n] is distance between node m and node n
				through[m][n]  =  -1;                //through[m,n] is intermediate (first nearest neighbour) of node m connecting to node n
			}
		}
	}		
	//Reinitializing diagnal elements as zeros
	#pragma omp parallel private(n)
	{
		#pragma omp for
		for(n = 0; n < N_floyd; ++n){
			dist[n][n]=0;                         //distance of itself is zero
		}
	}
    	//Running Floyd Warshall Algorithm for the spanning network
	for(k = 0; k < N_floyd; ++k)
    	{
		// #pragma omp parallel for private(m)
		for(m = 0; m < N_floyd; ++m)
		{
			for(n = 0; n < N_floyd; ++n)
			{
				// If The k node or target node is infinity (INF), then we are just going to overflow the numbers. Ignore them.
				// We can also safely ignore nodes between to themselves
				if (dist[m][k] == INF || dist[k][n] == INF || m == n){
				continue;
			}
			// Check if there is a faster path through node k
			int new_dist = dist[m][k] + dist[k][n];
			if(dist[m][n] <= new_dist){
				continue;
			}
			// This way is faster! Update the min distance and keep track of the node
			dist[m][n] = new_dist;
			through[m][n] = k;
			}
		}
    	}
}

double switch_func(double dist, double dist_D0, double dist_R0)
{
	int n=12; 
	int m=24;
 	return ((dist - dist_D0) == dist_R0) ? double(n)/m : (1-pow((dist-dist_D0)/dist_R0, n))/(1-(pow((dist-dist_D0)/dist_R0, m)));
}

double inv_switch_func(double dist, double dist_D0, double dist_R0)
{
	int n=12; 
	int m=24;
 	return ((dist - dist_D0) == dist_R0) ? double(m)/n : (1-pow((dist-dist_D0)/dist_R0, m))/(1-(pow((dist-dist_D0)/dist_R0, n)));
}

/* Calculate the first Legendre polynomial autocorrelation for vectors <cos(theta(t))> */
double_vec_t rot_first_ACF_func(vector<coordinates> x, int max_lag)
{
    	int size_x;
    	size_x = x.size();
    	double_vec_t R;
    	R.resize(max_lag);
    	double sum;
    	int i,j;
    	for (i=0;i<max_lag;++i) {
		sum=0;
		for (j=0;j<size_x-i;++j) {
			sum += dot(x[j],x[j+i])/(magnitude(x[j])*magnitude(x[j+i]));
		}
		R[i]=sum/(size_x-i);
    	}
    	if (R[0] != 0){
		for (i=1;i<max_lag;++i) {
			R[i] = R[i]/R[0];
		}
		R[0]=1;
    	}
    	return R;
}

/* Calculate the second Legendre polynomial autocorrelation for vectors <3/2*cos(theta(t))^2-1/2> */
double_vec_t rot_second_ACF_func(vector<coordinates> x, int max_lag)
{
    	int size_x=x.size();
    	double_vec_t R(max_lag);
    	double sum;
    	int i,j;

    	for (i=0;i<max_lag;++i) {
		sum=0;
		for (j=0;j<size_x-i;++j) {
			sum += 1.5* pow(dot(x[j],x[j+i])/(magnitude(x[j])*magnitude(x[j+i])), 2)-0.5;
		}
		R[i]=sum/(size_x-i);
    	}
    	if (R[0] != 0){
		for (i=1;i<max_lag;++i) {
			R[i] = R[i]/R[0];
		}
		R[0]=1;
    	}
    	return R;
}
/* Calculate the autocorrelation for hbond of each water molecule <h(0)h(t)> */
double_vec_t gen_ACF_func(double_vec_t x, int max_lag)
{
    	int size_x=x.size();
    	double_vec_t R(max_lag);
    	double sum;
    	int i,j;

    	for (i=0;i<max_lag;++i) {
		sum=0;
		for (j=0;j<size_x-i;++j) {
			sum+=x[j]*x[j+i];    
		}
		R[i]=sum/(size_x-i);
    	}
    	if (R[0] != 0){
		for (i=1;i<max_lag;++i) {
			R[i] = R[i]/R[0];
		}
		R[0]=1;
    	}
    	return R;
}
/* Calculate the autocorrelation for hbond of each water molecule <h(0)H(t)> */
double_vec_t hb_continuous_ACF_func(double_vec_t x, int max_lag)
{
    	short size_x=x.size();

    	vector<double> R(max_lag);
    	double sum;
    	int i,j;

    	for (i=0;i<max_lag;++i) {
		sum=0;
		//#pragma omp parallel private(i, j)
		//{
		//#pragma omp for
        	for (j=0;j<size_x-i;++j) {
			double ttt = 1;
			for(int rrr=j; rrr<=j+i; ++rrr) 
			{ 
				if ( x[rrr] < 1.0 ) { ttt = 0; break; }
			}
			sum+=x[j] * ttt;
		}
		R[i]=sum/(size_x-i);
		//}
    	}
    	if (R[0] != 0){
		for (i=1;i<max_lag;++i) {
			R[i] = R[i]/R[0];
		}
		R[0]=1;
    	}

    	return R;
}
// exponential fitting ---out[0]*exp(out[1]*t)----------
vector<double> exp_fit(vector<double> x, vector<double> y)
{
	int n=x.size();
	int i;
	double a, b, c;
	vector<double> lny,exp_coeff_out;
	lny.resize(n);
	exp_coeff_out.resize(2);
	for (i=0;i<n;i++)                        //Calculate the values of ln(yi)
		lny[i]=log(y[i]);        
	double xsum=0,x2sum=0,ysum=0,xysum=0;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
	for (i=0;i<n;i++)
	{
		xsum=xsum+x[i];                        //calculate sigma(xi)
		ysum=ysum+lny[i];                        //calculate sigma(yi)
		x2sum=x2sum+pow(x[i],2);                //calculate sigma(x^2i)
		xysum=xysum+x[i]*lny[i];                    //calculate sigma(xi*yi)
	}
	a=(n*xysum-xsum*ysum)/(n*x2sum-xsum*xsum);            //calculate slope(or the the power of exp)
	b=(x2sum*ysum-xsum*xysum)/(x2sum*n-xsum*xsum);            //calculate intercept
	c=pow(2.71828,b);                        //since b=ln(c)
	exp_coeff_out[0]=c;
	exp_coeff_out[1]=a;
	return exp_coeff_out;
}

public:
int main(FILE*in,FILE*out){

	int first_frame;
	int last_frame;
	string xtcfile;
	string tprfile; 
	string ndxfile;
	string outfile; 
	int k_max;
	string direct_id;
	string weight_id; 
	string ring_id;
	string print_path_id;
	string shell_id;
	double pop_shell;
	string scramble_id;
	int last_pro_res_ndx; 
	float hb_dist_cut; 
	float hb_ang_cut;
	float hb_weight_dist_cut; 
	double dist_DA_D0;
	double dist_DA_R0;
	double dist_HA_D0;
	double dist_HA_R0;
	double ang_D0;
	double ang_R0;
	int max_lag;
	double r_shell_1_N_O;
	double r_shell_1_C;
	double r_shell_bulk;
	double scrambling_ti;
	double scrambling_tf;
	int WDIM;

	double **adj_matrix;
	int **connection_matrix;
	int *defect_in, *defect_out;
	int N_matrix;
	int n_NH;
	int n_CO;
	int n_OW;
	
	string n_hb_grp;
	string h_hb_grp;
	string o_hb_grp;
	string ow_grp;
        string c_grp;
	string active_grp;
	string print_adjmatrix_id;
	string print_defect_id;
	string sys_grp = "System";
	string fpt_interface_grp;

	string fpt_id;
	string fpt_jump_id;
	double fpt_box_epsilon;
	double fpt_box_width;
	float fpt_time_width;

	clock_t tStart = clock();
	
	int frame;
	int i, j, k, m, n, ow_i;
	bool earlycheck;

	read_input(stdin, first_frame, last_frame, xtcfile, tprfile, ndxfile, outfile, 
		k_max, direct_id, weight_id, ring_id, print_path_id, shell_id, pop_shell, scramble_id, last_pro_res_ndx,
		hb_dist_cut, hb_ang_cut, hb_weight_dist_cut, dist_DA_D0, dist_DA_R0, dist_HA_D0, dist_HA_R0, ang_D0, ang_R0, 
		max_lag, r_shell_1_N_O, r_shell_1_C, r_shell_bulk, scrambling_ti, scrambling_tf,
		WDIM, n_hb_grp, h_hb_grp, o_hb_grp, ow_grp, c_grp, active_grp, print_adjmatrix_id, print_defect_id, 
		fpt_id, fpt_interface_grp, fpt_jump_id, fpt_box_epsilon, fpt_box_width, fpt_time_width);

	if(pop_shell > 1.0 )
	{
		pop_shell = 1.0;
	}else if(pop_shell <= 0.0)
	{	
		pop_shell = 0.000000000001;
	}
	
	fprintf(stdout,"%s %d\n","First frame step           :",first_frame);
	fprintf(stdout,"%s %d\n","Last frame step            :",last_frame);
	fprintf(stdout,"%s %s\n","Trajectory file            :",xtcfile.c_str());
	fprintf(stdout,"%s %s\n","Topology file              :",tprfile.c_str());
	fprintf(stdout,"%s %s\n","Index file                 :",ndxfile.c_str());
	fprintf(stdout,"%s %s\n","Output file                :",outfile.c_str());
	fprintf(stdout,"%s %d\n","Maximum shortest path      :",k_max);
	fprintf(stdout,"%s %s\n","direct_id                  :",direct_id.c_str());
	fprintf(stdout,"%s %s\n","weight_id                  :",weight_id.c_str());
	fprintf(stdout,"%s %s\n","ring_id                    :",ring_id.c_str());
	fprintf(stdout,"%s %s\n","print_path_id              :",print_path_id.c_str());
	fprintf(stdout,"%s %s\n","print_adjmatrix            :",print_adjmatrix_id.c_str());
	fprintf(stdout,"%s %s\n","print_defect               :",print_defect_id.c_str());
	fprintf(stdout,"%s %s\n","shell_id                   :",shell_id.c_str());
	fprintf(stdout,"%s %f\n","pop_shell                  :",pop_shell);
	fprintf(stdout,"%s %s\n","scramble_id                :",scramble_id.c_str());
	fprintf(stdout,"%s %d\n","last_pro_res_ndx           :",last_pro_res_ndx);
	fprintf(stdout,"%s %f\n","hbond_dist                 :",hb_dist_cut);
	fprintf(stdout,"%s %f\n","hbond_ang                  :",hb_ang_cut);
	fprintf(stdout,"%s %f\n","hbond_weight_dist          :",hb_weight_dist_cut);
	fprintf(stdout,"%s %f\n","dist_DA_D0                 :",dist_DA_D0);
	fprintf(stdout,"%s %f\n","dist_DA_R0                 :",dist_DA_R0);
	fprintf(stdout,"%s %f\n","dist_HA_D0                 :",dist_HA_D0);
	fprintf(stdout,"%s %f\n","dist_HA_R0                 :",dist_HA_R0);
	fprintf(stdout,"%s %f\n","ang_D0                     :",ang_D0);
	fprintf(stdout,"%s %f\n","ang_R0                     :",ang_R0);
	fprintf(stdout,"%s %d\n","max_lag                    :",max_lag);
	fprintf(stdout,"%s %f\n","r_shell_1_N_O              :",r_shell_1_N_O);
	fprintf(stdout,"%s %f\n","r_shell_1_C                :",r_shell_1_C);
	fprintf(stdout,"%s %f\n","r_shell_bulk               :",r_shell_bulk);
	fprintf(stdout,"%s %f\n","scrambling_ti              :",scrambling_ti);
	fprintf(stdout,"%s %f\n","scrambling_tf              :",scrambling_tf);
	fprintf(stdout,"%s %d\n","WDIM                       :",WDIM);
	fprintf(stdout,"%s %s\n","Hbond Donor groups         :",n_hb_grp.c_str());
	fprintf(stdout,"%s %s\n","Hbond Hydrogen groups      :",h_hb_grp.c_str());
	fprintf(stdout,"%s %s\n","Hbond Acceptor groups      :",o_hb_grp.c_str());
	fprintf(stdout,"%s %s\n","Solvent (Water) groups     :",ow_grp.c_str());
	fprintf(stdout,"%s %s\n","Heavy atoms in Active group:",c_grp.c_str());
	fprintf(stdout,"%s %s\n","All atoms in Active group  :",active_grp.c_str());
	fprintf(stdout,"%s %s\n","FPT analysis               :",fpt_id.c_str());
	fprintf(stdout,"%s %s\n","Interface group of FPT     :",fpt_interface_grp.c_str());
	fprintf(stdout,"%s %s\n","Jump Direction(x,y,z) FPT  :",fpt_jump_id.c_str());
	fprintf(stdout,"%s %f\n","Epsilon of box slice       :",fpt_box_epsilon);
	fprintf(stdout,"%s %f\n","Width of box slice         :",fpt_box_width);
	fprintf(stdout,"%s %f\n","Interval of time slice     :",fpt_time_width);

	Index ndx(ndxfile);
	string fn_old_length, fn_old_path;
	ifstream length_old_file, path_old_file;
	
	string fn_scrambling_time, fn_scrambling_check_node, fn_scrambling_check_conn, fn_scrambling_check_num;
	string fn_scrambling_len_acf, fn_scrambling_len_change;
	ofstream scrambling_time_file, scrambling_check_node_file, scrambling_check_conn_file, scrambling_check_num_file;
	ofstream scrambling_len_acf_file, scrambling_len_change_file;
	if(scramble_id.compare("yes") == 0){
		std::cout << ">>>The scrambling time of the hbond network wires will be calculated" << std::endl;
		ring_id = "no";
		print_path_id = "no";
		fn_old_length = "scrambling_length_tot";
		length_old_file.open(fn_old_length.c_str()); 
		if (!length_old_file.is_open()) {
			std::cout << "ERROR: Cannot open \"length\" file to calculate scrambling time, named \"scrambling_length_tot\"." << std::endl;
			return 0;
		}
		fn_old_path = "scrambling_path_tot";
		path_old_file.open(fn_old_path.c_str()); 
		if (!path_old_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"path\" file to calculate scrambling time, named \"scrambling_path_tot\"." << std::endl;
			return 0;
		}
		fn_scrambling_time = outfile+".scramble_time";
		scrambling_time_file.open(fn_scrambling_time.c_str()); 
		if (!scrambling_time_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"scrambling time\" file." << std::endl;
			return 0;
		}
		fn_scrambling_check_node = outfile+".scramble_check_node";
		scrambling_check_node_file.open(fn_scrambling_check_node.c_str()); 
		if (!scrambling_check_node_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"scrambling check node\" file." << std::endl;
			return 0;
		}
		fn_scrambling_check_conn = outfile+".scramble_check_conn";
		scrambling_check_conn_file.open(fn_scrambling_check_conn.c_str()); 
		if (!scrambling_check_conn_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"scrambling check connection\" file." << std::endl;
			return 0;
		}
		fn_scrambling_check_num = outfile+".scramble_check_num";
		scrambling_check_num_file.open(fn_scrambling_check_num.c_str()); 
		if (!scrambling_check_num_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"scrambling check connection\" file." << std::endl;
			return 0;
		}
		fn_scrambling_len_acf = outfile+".scramble_len_acf";
		scrambling_len_acf_file.open(fn_scrambling_len_acf.c_str()); 
		if (!scrambling_len_acf_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"scrambling length acf\" file." << std::endl;
			return 0;
		}
		fn_scrambling_len_change = outfile+".scramble_len_change";
		scrambling_len_change_file.open(fn_scrambling_len_change.c_str()); 
		if (!scrambling_len_change_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"scrambling length acf\" file." << std::endl;
			return 0;
		}
	}else if(scramble_id.compare("no") == 0 ){
		std::cout << ">>>The scrambling time of the hbond network wires will not be calculated" << std::endl;
	}else{
		std::cout << "ERROR: Please input correct scrambling time option (yes, no)" << std::endl;
		return 0;
	}
	
	//Checking output of path nodes before to open the path files
	if ( print_path_id.compare("yes") == 0 ){
		std::cout << ">>>The path file will be printed into the file" << std::endl;
		
	}else if ( print_path_id.compare("no") == 0 ){
		std::cout << ">>>The path file will not be printed into the file" << std::endl;
	}else{
		std::cout << "ERROR: Please input correct path output parameters (yes, no)" << std::endl;
		return 0;
	}
	
	//Output filenames
	string fn_path,fn_defect, fn_length, fn_weight;
	ofstream path_file, length_file, weight_file, defect_file;
	if(( ring_id.compare("DA") == 0 ) || ( ring_id.compare("DADDAA") == 0 ) ){
		std::cout << ">>>PATH between NH and CO will be calculated" << std::endl;
		fn_length = outfile + ".length";
		length_file.open(fn_length.c_str()); 
		if (!length_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"length\" file." << std::endl;
			return 0;
		}
		if (print_path_id.compare("yes") == 0 ){
			fn_path = outfile + ".path";	
			path_file.open(fn_path.c_str()); 
			if (!path_file.is_open()) {
				std::cout << "ERROR: Cannot open output \"path\" file." << std::endl;
				return 0;
			}
		}
		if ( direct_id.compare("undirected") == 0 && print_defect_id.compare("yes") == 0 ){
			fn_defect = outfile + ".defect";
			defect_file.open(fn_defect.c_str()); 
			if (!defect_file.is_open()) {
				std::cout << "ERROR: Cannot open output \"defect\" file." << std::endl;
				return 0;
			}
		}
		if ( weight_id.compare("weight") == 0 ){
			fn_weight = outfile + ".weight";
			weight_file.open(fn_weight.c_str()); 
			if (!weight_file.is_open()) {
				std::cout << "ERROR: Cannot open output \"weight\" file." << std::endl;
				return 0;
			}
			weight_file << fixed << setprecision(3);
		}
	}else if((ring_id.compare("no") == 0) || (ring_id.compare("DD") == 0) || (ring_id.compare("AA") == 0)){
		std::cout << ">>>PATH calculation between NH and CO is canceled" << std::endl;
	}else{
		std::cout << "ERROR: Please input correct ring input parameters(DA, DD, AA, DADDAA, no)" << std::endl;
		return 0;
	}

	string fn_ringnhnh_path, fn_ringnhnh_length, fn_ringnhnh_defect, fn_ringnhnh_weight;
	ofstream ringnhnh_path_file, ringnhnh_length_file, ringnhnh_defect_file, ringnhnh_weight_file;
	if((ring_id.compare("DD") == 0) || (ring_id.compare("DADDAA") == 0) ){
		std::cout << ">>>Ring between NH and NH will be calculated" << std::endl;
		fn_ringnhnh_length = outfile + ".ringnhnh_length";
		ringnhnh_length_file.open(fn_ringnhnh_length.c_str()); 
		if (!ringnhnh_length_file.is_open()) {
		      std::cout << "ERROR: Cannot open output \"ringnhnh_length\" file." << std::endl;
		      return 0;
		}
		if ( print_defect_id.compare("yes") == 0 ){
			fn_ringnhnh_defect = outfile + ".ringnhnh_defect";
			ringnhnh_defect_file.open(fn_ringnhnh_defect.c_str()); 
			if (!ringnhnh_defect_file.is_open()) {
			      std::cout << "ERROR: Cannot open output \"ringnhnh_defect\" file." << std::endl;
			      return 0;
			}
		}
		if (print_path_id.compare("yes") == 0 ){
			fn_ringnhnh_path = outfile + ".ringnhnh_path";
			ringnhnh_path_file.open(fn_ringnhnh_path.c_str()); 
			if (!ringnhnh_path_file.is_open()) {
				std::cout << "ERROR: Cannot open output \"ringnhnh_path\" file." << std::endl;
				return 0;
			}
		}
		if ( weight_id.compare("weight") == 0 ){
			fn_ringnhnh_weight = outfile + ".ringnhnh_weight";
			ringnhnh_weight_file.open(fn_ringnhnh_weight.c_str()); 
			if (!ringnhnh_weight_file.is_open()) {
				std::cout << "ERROR: Cannot open output \"ringnhnh_weight\" file." << std::endl;
				return 0;
			}
			ringnhnh_weight_file << fixed << setprecision(3);
		}
	}else if((ring_id.compare("no") == 0) || (ring_id.compare("DA") == 0) || (ring_id.compare("AA") == 0)){
		std::cout << ">>>Ring calculation between NH and NH is canceled" << std::endl;
	}else{
		std::cout << "ERROR: Please input correct ring input parameters(DA, DD, AA, DADDAA, no)" << std::endl;
		return 0;
	}
	
	string fn_ringcoco_path, fn_ringcoco_length, fn_ringcoco_defect, fn_ringcoco_weight;
	ofstream ringcoco_path_file, ringcoco_length_file, ringcoco_defect_file, ringcoco_weight_file;
	if((ring_id.compare("AA") == 0) || (ring_id.compare("DADDAA") == 0) ){
		std::cout << ">>>Ring between CO and CO will be calculated" << std::endl;
		fn_ringcoco_length = outfile + ".ringcoco_length";
		ringcoco_length_file.open(fn_ringcoco_length.c_str()); 

		if (!ringcoco_length_file.is_open()) {
		      std::cout << "ERROR: Cannot open output \"ringcoco_length\" file." << std::endl;
		      return 0;
		}
		if ( print_defect_id.compare("yes") == 0 ){
			fn_ringcoco_defect = outfile + ".ringcoco_defect";
			ringcoco_defect_file.open(fn_ringcoco_defect.c_str()); 
			if (!ringcoco_defect_file.is_open()) {
			      std::cout << "ERROR: Cannot open output \"ringcoco_defect\" file." << std::endl;
			      return 0;
			}
		}
		if (print_path_id.compare("yes") == 0 ){
			fn_ringcoco_path = outfile + ".ringcoco_path";		
			ringcoco_path_file.open(fn_ringcoco_path.c_str()); 
			if (!ringcoco_path_file.is_open()) {
				std::cout << "ERROR: Cannot open output \"ringcoco_path\" file." << std::endl;
				return 0;
			}
		}
		if ( weight_id.compare("weight") == 0 ){
			fn_ringcoco_weight = outfile + ".ringcoco_weight";
			ringcoco_weight_file.open(fn_ringcoco_weight.c_str()); 
			if (!ringcoco_weight_file.is_open()) {
				std::cout << "ERROR: Cannot open output \"ringcoco_weight\" file." << std::endl;
				return 0;
			}
			ringcoco_weight_file << fixed << setprecision(3);
		}
	}else if((ring_id.compare("no") == 0) || (ring_id.compare("DA") == 0) || (ring_id.compare("DD") == 0)){
		std::cout << ">>>Ring calculation between CO and CO is canceled" << std::endl;
	}else{
		std::cout << "ERROR: Please input correct ring input parameters (DA, DD, AA, DADDAA, no)" << std::endl;
		return 0;
	}

	string fn_prohb, fn_spanning, fn_shell_num, fn_shell_ndx, fn_shell_defect, fn_shell_inout, fn_water_defect, fn_eig;
	string fn_dpl, fn_shell_hb_span, fn_shell_acf, fn_shell_q, fn_shell_lsi;
	ofstream shell_defect_file, shell_num_file, shell_inout_file, water_defect_file;
	ofstream dpl_file,  shell_ndx_file, hb_span_file, shell_acf_file, shell_q_file, shell_lsi_file;
	if(shell_id.compare("yes") ==0){
		std::cout << ">>>The properites of first shell water network will be calculated" << std::endl;
		fn_dpl		  = outfile + ".rad_dpl";
		fn_shell_hb_span  = outfile + ".shell_hb_span";;
		fn_shell_ndx   	  = outfile + ".shell_ndx";
		fn_shell_acf      = outfile + ".shell_acf";
		fn_shell_q	  = outfile + ".shell_q";
		fn_shell_lsi      = outfile + ".shell_lsi";
		dpl_file.open(fn_dpl.c_str()); 
		if (!dpl_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"Water radial dipole-dipole moment correlation function \" file." << std::endl;
			return 0;
		}
		dpl_file << fixed << setprecision(3);
		shell_ndx_file.open(fn_shell_ndx.c_str()); 
		if (!shell_ndx_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"Water index \" file." << std::endl;
			return 0;
		}
		hb_span_file.open(fn_shell_hb_span.c_str()); 
		if (!hb_span_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"Water H-bond around chain A\" file." << std::endl;
			return 0;
		}
		shell_acf_file.open(fn_shell_acf.c_str()); 
		if (!shell_acf_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"ACF\" file." << std::endl;
			return 0;
		}
		shell_q_file.open(fn_shell_q.c_str());
                if (!shell_q_file.is_open()) {
                        std::cout << "ERROR: Cannot open output \"q\" file." << std::endl;
                        return 0;
                }
                shell_lsi_file.open(fn_shell_lsi.c_str());
                if (!shell_lsi_file.is_open()) {
                        std::cout << "ERROR: Cannot open output \"LSI\" file." << std::endl;
                        return 0;
                }
	}else if(shell_id.compare("no") == 0 ){
		std::cout << ">>>The properites of water molecules and hbond network within 1st and 2nd shells and bulk will not be calculated" << std::endl;
	}else{
		std::cout << "ERROR: Please input correct shell input parameters (yes, no)" << std::endl;
		return 0;
	}
	
	string fn_defect_all; 
	ofstream defect_all_file;
	if ( print_defect_id.compare("yes") == 0 ){
                std::cout << ">>>The defect information of all donor and acceptor groups will be printed into the file" << std::endl;
		fn_defect_all = outfile + ".defect_all";
		defect_all_file.open(fn_defect_all.c_str());
		if (!defect_all_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"defect\" file." << std::endl;
			return 0;
		}
        }else if ( print_defect_id.compare("no") == 0 ){
                std::cout << ">>>The defct information of all donor and acceptor groups will not be printed into the file" << std::endl;
        }else{
                std::cout << "ERROR: Please input correct defect print parameters (T, F)" << std::endl;
                return 0;
        }
	
	string fn_adjmatrix;
	ofstream adjmatrix_file;
	if ( print_adjmatrix_id.compare("yes") == 0 ){
                std::cout << ">>>The adjcency matrix will be printed into the file" << std::endl;
		fn_adjmatrix = outfile + ".adjmatrix";	
		adjmatrix_file.open(fn_adjmatrix.c_str()); 
		if (!adjmatrix_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"adjcency matrix\" file." << std::endl;
			return 0;
		}
        }else if ( print_adjmatrix_id.compare("no") == 0 ){
                std::cout << ">>>The defct information of all donor and acceptor groups will not be printed into the file" << std::endl;
        }else{
                std::cout << "ERROR: Please input correct adjcency matrix print parameters (T, F)" << std::endl;
                return 0;
        }
	
	string fn_fpt_t, fn_fpt_sign, fn_fpt_info, fn_fpt_ndx;
	ofstream fpt_t_file, fpt_sign_file, fpt_info_file, fpt_ndx_file;
	if ( fpt_id.compare("yes") == 0 ){
                std::cout << ">>>First passage time along " << fpt_jump_id << "-direction" << " will be calculated" << std::endl;
		fn_fpt_t = outfile + ".fpt_" + fpt_jump_id + "_t";
		fpt_t_file.open(fn_fpt_t.c_str()); 
		if (!fpt_t_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"fpt_t\" file." << std::endl;
			return 0;
		}
		fn_fpt_sign = outfile + ".fpt_" + fpt_jump_id + "_sign";
		fpt_sign_file.open(fn_fpt_sign.c_str()); 
		if (!fpt_sign_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"fpt_sign\" file." << std::endl;
			return 0;
		}
		fn_fpt_info = outfile + ".fpt_" + fpt_jump_id + "_info";
		fpt_info_file.open(fn_fpt_info.c_str()); 
		if (!fpt_info_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"fpt_info\" file." << std::endl;
			return 0;
		}
		fpt_info_file << "# time_id, box_id, box_center, time_window, epsilon_box, width_box, initial_number, box_type" << std::endl;
		fn_fpt_ndx = outfile + ".fpt_" + fpt_jump_id + "_ndx";
		fpt_ndx_file.open(fn_fpt_ndx.c_str()); 
		if (!fpt_ndx_file.is_open()) {
			std::cout << "ERROR: Cannot open output \"fpt_num\" file." << std::endl;
			return 0;
		}
        }else if ( fpt_id.compare("no") == 0 ){
                std::cout << ">>>First passage time will not be calculated" << std::endl;
        }else{
                std::cout << "ERROR: Please fpt parameters (T, F)" << std::endl;
                return 0;
        }

	//Creating a Trajectory object as a pointer
	Trajectory *traj = new Trajectory(xtcfile, ndx);
	
	//Input values of Shortest Path
	n_NH = traj->GetNAtoms(n_hb_grp);
	unsigned int n_hNH = traj->GetNAtoms(h_hb_grp);
	if(n_NH != n_hNH){
		std::cout << "Please check number of N and H atoms of NH donor groups in index file \"which should be the same\"" << std::endl;
		return 0;
	}
	n_CO = traj->GetNAtoms(o_hb_grp);
	n_OW = (traj->GetNAtoms(ow_grp))/WDIM;
        int n_C = traj->GetNAtoms(c_grp);
	int n_SYS = traj->GetNAtoms(sys_grp);
			
	string fn_bulk_path1, fn_bulk_length1, fn_bulk_defect1, fn_bulk_weight1, fn_bulk_el_force1;
	string fn_bulk_path2, fn_bulk_length2, fn_bulk_defect2, fn_bulk_weight2, fn_bulk_el_force2;

	ofstream bulk_path_file1, bulk_length_file1, bulk_defect_file1, bulk_weight_file1, bulk_el_force_file1;
	ofstream bulk_path_file2, bulk_length_file2, bulk_defect_file2, bulk_weight_file2, bulk_el_force_file2;
		
	if( (n_SYS - n_OW*WDIM) == 0 )
	{
		if(n_NH == 0 && n_CO == 0)
		{
			ring_id = "bulk";
			std::cout << "Bulk water : Number of water molecules " << n_OW << std::endl;
			std::cout << "\"ring_id\" is set to \"bulk\""  << std::endl << std::endl;
			if (shell_id.compare("yes") == 0)
			{
				std::cout << " Shell analysis is cancelled out because the system is bulk water !!! " << std::endl;
			}	
			
			/*	
 			std::cout << ">>>PATH between each water molecule to its first and second shell water moelcues  will be calculated" << std::endl;
			fn_bulk_length1 = outfile + ".bulk_length_1st";
			bulk_length_file1.open(fn_bulk_length1.c_str());
			if(!bulk_length_file1.is_open()) {
				std::cout << "ERROR: Cannot open output \"first-shell bulk length\" file." << std::endl;
				return 0;
			}
			fn_bulk_length2 = outfile + ".bulk_length_2nd";
			bulk_length_file2.open(fn_bulk_length2.c_str());
			if(!bulk_length_file2.is_open()) {
				std::cout << "ERROR: Cannot open output \"second-shell bulk length\" file." << std::endl;
				return 0;
			}
			fn_bulk_el_force1 = outfile + ".bulk_el_force_1st";
			bulk_el_force_file1.open(fn_bulk_el_force1.c_str());
			if(!bulk_el_force_file1.is_open()) {
				std::cout << "ERROR: Cannot open output \"first-shell bulk electric field\" file." << std::endl;
				return 0;
			}
			fn_bulk_el_force2 = outfile + ".bulk_el_force_2nd";
			bulk_el_force_file2.open(fn_bulk_el_force2.c_str());
			if(!bulk_el_force_file2.is_open()) {
				std::cout << "ERROR: Cannot open output \"second-shell bulk electric field\" file." << std::endl;
				return 0;
			}
			if(print_path_id.compare("yes") == 0 ){
				fn_bulk_path1= outfile + ".bulk_path_1st";
				bulk_path_file1.open(fn_bulk_path1.c_str());
				if(!bulk_path_file1.is_open()) {
					std::cout << "ERROR: Cannot open output \"first-shell bulk path\" file." << std::endl;
					return 0;
				}
				fn_bulk_path2= outfile + ".bulk_path_2nd";
				bulk_path_file2.open(fn_bulk_path2.c_str());
				if(!bulk_path_file2.is_open()) {
					std::cout << "ERROR: Cannot open output \"second-shell bulk path\" file." << std::endl;
					return 0;
				}
			}
			if( direct_id.compare("undirected") == 0 ){
				fn_bulk_defect1 = outfile + ".bulk_defect_1st";
				bulk_defect_file1.open(fn_bulk_defect1.c_str());
				if(!bulk_defect_file1.is_open()) {
					std::cout << "ERROR: Cannot open output \"first-shell bulk defect\" file." << std::endl;
					return 0;
				}
				fn_bulk_defect2 = outfile + ".bulk_defect_2nd";
				bulk_defect_file2.open(fn_bulk_defect2.c_str());
				if(!bulk_defect_file2.is_open()) {
					std::cout << "ERROR: Cannot open output \"second-shell bulk defect\" file." << std::endl;
					return 0;
				}
			}
			if( weight_id.compare("weight") == 0 ){
				fn_bulk_weight1 = outfile + ".bulk_weight_1st";
				bulk_weight_file1.open(fn_bulk_weight1.c_str());
				if(!bulk_weight_file1.is_open()) {
					std::cout << "ERROR: Cannot open output \"first-shell bulk weight\" file." << std::endl;
					return 0;
				}
				bulk_weight_file1 << fixed << setprecision(3);
				fn_bulk_weight2 = outfile + ".bulk_weight_2nd";
				bulk_weight_file2.open(fn_bulk_weight2.c_str());
				if(!bulk_weight_file2.is_open()) {
					std::cout << "ERROR: Cannot open output \"second-shell bulk weight\" file." << std::endl;
					return 0;
				}
				bulk_weight_file2 << fixed << setprecision(3);
			}
			*/
		}
		else
		{
			std::cout << "Bulk water : Number of water molecules " << n_OW << std::endl;
			std::cout << "Please remove atoms in groups of " << n_hb_grp << " and " << o_hb_grp << " !!! "<< std::endl;
			return 0;
		}
	}
	else
	{
		std::cout << "Number of H-bond Donors       " << n_NH << std::endl;
		std::cout << "Number of H-bond Acceptors    " << n_CO << std::endl;
		std::cout << "Number of Water Molecules     " << n_OW << std::endl;
		std::cout << "Number of Hydrophobic atoms   " << n_C << std::endl;
		std::cout << std::endl;
	}
	std::cout << ">>>>>> Please pay attention for the starting index of water molecule <<<<" << std::endl;
	std::cout << "Residue Index of First Water  " << last_pro_res_ndx+1 << std::endl;
	std::cout << std::endl;
	std::cout << ">>>>>> To visualize the water nodes with residue mode in VMD <<<<" << std::endl;
	std::cout << std::endl;
	
	// Reading whole Trajectory
        traj->read(first_frame, 1, last_frame);
        double dt = traj->GetTime(1) - traj->GetTime(0);
        std::cout << " Time step [ps] : " << dt << std::endl;
        last_frame=first_frame+traj->GetNFrames(); //if last_frame is beyond real last frame, it's easily corrected by # of the loaded frames
        std::cout << " # of frames to be read : " << traj->GetNFrames() << " >>corrected last frame :" << last_frame << std::endl;
        int n_tot_frame = last_frame - first_frame;
        std::cout << " Reading FINISHED >>>"<< "TOTAL FRAMES = " << traj->GetNFrames() << " | TOTAL TIME = "<< (double)(clock() - tStart)/CLOCKS_PER_SEC << std::endl;

        //Initialization of Adjacent Matrix

	N_matrix = n_NH + n_CO + n_OW;
	adj_matrix = new double*[N_matrix];
	connection_matrix = new int*[N_matrix];
	for(i =0; i<N_matrix; ++i){
		adj_matrix[i] = new double[N_matrix];
		connection_matrix[i] = new int[N_matrix];
	}
	defect_in = new int[N_matrix];
	defect_out = new int[N_matrix];
	//int *spanning_node = new int[n_OW];
	
	vector <coordinates> xyz_hNH, xyz_nNH, xyz_CO, xyz_OW;
	vector <coordinates> xyz_heavy_C;
	triclinicbox box;
        
	int n_protein = traj->GetNAtoms(active_grp);
	Topology top(tprfile,ndx);	
	vector <double> mass_protein;
	vector <double> charge_protein;
	vector <coordinates> xyz_active;
	mass_protein = top.GetMass(active_grp);
	charge_protein = top.GetCharge(active_grp);
	coordinates com_protein;
	coordinates dpl_protein;
	
        
	
	//Initilization of Graph and Path
	Graph graph;
 	Graph::path_list paths; // Find the k shortest paths

 	int *path_len, *ring_node, *end_node_vec;
	int n_dd, n_aa;
	unsigned int node_ndx;

        int path_max = 0;
        int N_tot_path = n_NH*n_CO;

	//Initialization of shell analysis
        //int **first_shell, **second_shell, **bulk;
	vector<vector <int>> first_shell, first_shell_N, first_shell_O, first_shell_C;
	vector<vector <double>> ow_dist_matrix;
	if( shell_id.compare("yes") == 0 ){
            	//Initialization of shell coordinates
		first_shell.resize(n_tot_frame);
		first_shell_N.resize(n_tot_frame);
		first_shell_O.resize(n_tot_frame);
		first_shell_C.resize(n_tot_frame);
		//second_shell = new int*[n_tot_frame];
	        //bulk = new int*[n_tot_frame];
            	for(i =0; i<(n_tot_frame); ++i){
                    	first_shell[i].resize(n_OW);
                    	first_shell_N[i].resize(n_OW);
                    	first_shell_O[i].resize(n_OW);
                    	first_shell_C[i].resize(n_OW);
			//second_shell[i] = new int[n_OW];
                    	//bulk[i] = new int[n_OW];
            	}
		ow_dist_matrix.resize(n_OW);
	    	for(i =0; i<n_OW; ++i){
			ow_dist_matrix[i].resize(n_OW);
			for(j=0; j<n_OW;++j){
				ow_dist_matrix[i][j] = 0;
			}
	    	}
        }

        
	//Initialization of scrambling time
        int scrambling_tot_frame;
        scrambling_tot_frame = (scrambling_tf-scrambling_ti)/dt+1;
        int **old_length_data, ***old_path_data, ***old_path_data_origin;
        double **scramble_t, *scrambling_time_series;
        bool **scramble_t_bool;

	if(scramble_id.compare("yes") == 0){
                //Initialization of length and path coordinates for the scrambling time
                old_length_data = new int*[scrambling_tot_frame];
                scramble_t = new double*[scrambling_tot_frame];
                scramble_t_bool = new bool*[scrambling_tot_frame];
                scrambling_time_series = new double[scrambling_tot_frame];
                for(i=0; i<scrambling_tot_frame; ++i){
                        old_length_data[i] = new int[N_tot_path];
                        scramble_t[i] = new double[N_tot_path];
                        scramble_t_bool[i] = new bool[N_tot_path];
                        for(j=0; j<N_tot_path; ++j){ 
                                scramble_t[i][j] = 0; 
                                scramble_t_bool[i][j] = true; 
                        }
                }
            // 	//std::cout << " Scramble Length file is being allocated ... " << std::endl;
            // 	int ***scramble_len = new int**[n_tot_frame]; //this recorde all the scramble lengths of path on each frame
            // 											//len[ each frame ][ all previous frames ] [ path index 
            // 	for(i=0; i<n_tot_frame; ++i){
            // 		scramble_len[i] = new int*[scrambling_tot_frame];
            // 		for(j=0; j<scrambling_tot_frame; ++j){
            // 			scramble_len[i][j] = new int[N_tot_path];
            // 			for(k = 0; k<N_tot_path; ++k){
            // 				scramble_len[i][j][k] = 0;
            // 			}
            // 		}
            // 	}
            // 	
            // 	//std::cout << " Scramble Length file is allocated correltly !" << std::endl;
		std::cout << "******************* This scrambling time data should be checked carefully ****************" << std::endl;
		std::cout << ">> Initial Time (ps) of \"scrambling length and path data\" is     >> " << scrambling_ti << std::endl;
		std::cout << ">> Final Time (ps)  of \"scrambling length and path data\" is      >> " << scrambling_tf << std::endl;
		std::cout << ">> Time step  (ps)  of \"scrambling length and path data\" is      >> " << dt << std::endl;
		std::cout << ">> Total number of frames of scrambling length and path nodes file  : " << scrambling_tot_frame << std::endl << std::endl;
		path_max = 0;
		for(i=0; i < scrambling_tot_frame; ++i){ 
			double tt_t;
			length_old_file >> tt_t; length_old_file >> tt_t;
			for(j=0; j< N_tot_path; ++j){ 
				length_old_file >> old_length_data[i][j]; 
				if(path_max < old_length_data[i][j]){ path_max = old_length_data[i][j];}
				////std::cout << "frame" << frame << " path id " << i << " length " << old_length_data[frame][i] << std::endl;
			}
		}
		std::cout << "MAX PATH LENGTH " << path_max << std::endl;
                old_path_data = new int**[scrambling_tot_frame];
                old_path_data_origin = new int**[scrambling_tot_frame];
                for(i=0; i<(scrambling_tot_frame); ++i){
                        old_path_data[i] = new int*[N_tot_path];
                        old_path_data_origin[i] = new int*[N_tot_path];
                        for(j=0; j<N_tot_path; ++j){
                                old_path_data[i][j] = new int[path_max+1];
                                old_path_data_origin[i][j] = new int[path_max+1];
                                for(k = 0; k<(path_max+1); ++k){
                                        old_path_data[i][j][k] = 0;
                                        old_path_data_origin[i][j][k] = 0;
                                }
                        }
                }
                
		for(i=0; i<scrambling_tot_frame; ++i){
			int tt_t;
			string tt_str;
			path_old_file >> scrambling_time_series[i];
			////std::cout << " scrambling file index " <<  i << " " << scrambling_time_series[i] << " " << std::endl;
			for(j=0; j< N_tot_path; ++j){
				path_old_file >> tt_str; 
				////std::cout << tt_str << " ";
				path_old_file >> tt_str; 
				////std::cout << tt_str << " ";
				if( old_length_data[i][j] > 0) {
					path_old_file >> tt_t; path_old_file >> tt_str;
					////std::cout << tt_t << " " << tt_str << " ";
					path_old_file >> tt_t; 
					old_path_data[i][j][0] = (int)j/n_CO;
					old_path_data_origin[i][j][0] = tt_t;
					////std::cout << tt_t << " ";
					////std::cout << "time " << i << "path id " << j << "node " << old_path_data[i][j][0] << " original : " << tt_t << std::endl;
					for(k = 1; k < old_length_data[i][j]; ++k){
						path_old_file >> tt_t;
						old_path_data[i][j][k] = tt_t - (last_pro_res_ndx+1)+n_NH;
						old_path_data_origin[i][j][k] = tt_t;
						////std::cout << "time " << i << "path id " << j << "node " << old_path_data[i][j][k] << " original : " << tt_t << std::endl;
					}
					path_old_file >> tt_t; old_path_data[i][j][old_length_data[i][j]] = j % n_CO + n_NH + n_OW;
					old_path_data_origin[i][j][old_length_data[i][j]] = tt_t;
					////std::cout << "time " << i << "path id " << j << "node " << old_path_data[i][j][old_length_data[i][j]] << " original : " << tt_t << std::endl;
				}
			}
		}
	}
	
	//Real Calculation
	int scrambling_check_time = 0; 
	int scrambling_check_path = 1;
	double current_time;
	std::cout << "first frame " << first_frame << " last frame " << last_frame << std::endl;
	double max_box_size = 0.0;
	
	//HB-network, Hydration shell analysis, Scrambling analysis
	for (frame = 0; frame < n_tot_frame; ++frame){ 
		if((ring_id.compare( "no" ) == 0) && (shell_id.compare("no") == 0) && (scramble_id.compare("no") == 0)) break;
	    	//Getting information and trajectory files for each frame
		box = traj->GetBox(frame);
		double box_size = magnitude(coordinates(box(X,X),box(X,Y),box(X,Z))+coordinates(box(Y,X),box(Y,Y),box(Y,Z))+coordinates(box(Z,X),box(Z,Y),box(Z,Z)));
		if ( box_size > max_box_size ) max_box_size = box_size;
		current_time = traj->GetTime(frame);
		xyz_nNH = traj->GetXYZ(frame, n_hb_grp);
		xyz_hNH = traj->GetXYZ(frame, h_hb_grp);
	   	xyz_CO = traj->GetXYZ(frame, o_hb_grp);
	   	xyz_OW = traj->GetXYZ(frame, ow_grp);
		xyz_heavy_C = traj->GetXYZ(frame, c_grp);
		
		if(weight_id.compare( "binary") == 0){
			adjcent_matrix_calc(direct_id, box, xyz_nNH, xyz_hNH, xyz_CO, xyz_OW, adj_matrix, connection_matrix, defect_in, defect_out, hb_dist_cut, hb_ang_cut, WDIM);
			
			if(print_adjmatrix_id.compare("yes") == 0) {
                                for(j =0;j<n_CO;++j){   
                                        for(k =0;k<n_NH;k++){ 
                                                if( ndx.GetLocation(o_hb_grp, j) == ndx.GetLocation(n_hb_grp, k) )
                                                {       
                                                        std::cout << "the same atom of D and A : " << ndx.GetLocation(o_hb_grp, j) << " should be fixed as INF!!!" << std::endl;
                                                        adj_matrix[k][n_OW+n_NH+j] = INF; 
                                                        connection_matrix[k][n_OW+n_NH+j] = 0; 
                                                        if(direct_id.compare( "undirected") == 0){
                                                                adj_matrix[n_OW+n_NH+j][k] = INF;
                                                        }
                                                }
                                        }
                                }
                                for(j=0;j<N_matrix; ++j)
                                {       
                                        for(k=0;k<N_matrix; ++k)
                                        {       
                                                if(adj_matrix[j][k] >= INF )
                                                {       
                                                        //adjmatrix_file << j+1 << " " << k+1 << " 0 " << std::endl;
                                                }
                                                else    
                                                {       
                                                        adjmatrix_file << j+1 << " " << k+1 << " " << adj_matrix[j][k] << std::endl;
                                                }
                                        }
                                }
			}
			if(print_defect_id.compare("yes") == 0) {
                                for(j=0;j<N_matrix; ++j)
                                {       
					defect_all_file << j+1 << " " << defect_in[j] << " " << defect_out[j] << std::endl;
                                }
			}

		}else if(weight_id.compare( "weight") == 0){
	        	adjcent_matrix_weight_calc(direct_id, box, xyz_nNH, xyz_hNH, xyz_CO, xyz_OW, adj_matrix, connection_matrix, defect_in, defect_out, hb_weight_dist_cut, \
						  dist_DA_D0, dist_DA_R0, dist_HA_D0, dist_HA_R0, ang_D0, ang_R0, WDIM);
			if (print_adjmatrix_id.compare("yes") == 0) {
                                for(j =0;j<n_CO;++j){
                                        for(k =0;k<n_NH;k++){
                                                if( ndx.GetLocation(o_hb_grp, j) == ndx.GetLocation(n_hb_grp, k) )
                                                {
                                                        std::cout << "the same atom of D and A : " << ndx.GetLocation(o_hb_grp, j) << " should be fixed as INF!!!" << std::endl;
                                                        adj_matrix[k][n_OW+n_NH+j] = INF;
                                                        connection_matrix[k][n_OW+n_NH+j] = 0;
                                                        if(direct_id.compare( "undirected") == 0){
                                                                adj_matrix[n_OW+n_NH+j][k] = INF;
                                                        }
                                                }
                                        }
                                }
                                for(j=0;j<N_matrix; ++j)
                                {       
                                        for(k=0;k<N_matrix; ++k)
                                        {       
                                                if(adj_matrix[j][k] >= INF )
                                                {       
                                                        //adjmatrix_file << j+1 << " " << k+1 << " 0 " << std::endl;
                                                }
                                                else
                                                {
                                                        adjmatrix_file << j+1 << " " << k+1 << " " << double(1.0)/adj_matrix[j][k] << std::endl;
                                                }
                                        }
                                }
			}
			if(print_defect_id.compare("yes") == 0) {
                                for(j=0;j<N_matrix; ++j)
                                {       
					defect_all_file << j+1 << " " << defect_in[j] << " " << defect_out[j] << std::endl;
                                }
			}
		}
		
		//Calculating the adjacency matrix for D-A pairs
		if((ring_id.compare("DA") == 0) || (ring_id.compare("DADDAA") == 0) ){
			if( print_path_id.compare("yes") == 0 ) { path_file << current_time << std::endl; }
			length_file << current_time << " " ;
			if( print_defect_id.compare("yes") == 0 ) { defect_file << current_time << " "; }
			weight_file << current_time << " " ;
		}
		if((ring_id.compare("DD") == 0) || (ring_id.compare("DADDAA") == 0) ){
			if( print_path_id.compare("yes") == 0 ) { ringnhnh_path_file << current_time << std::endl;}
			ringnhnh_length_file << current_time << " " ;
			if( print_defect_id.compare("yes") == 0 ) {ringnhnh_defect_file << current_time << " "; }
			ringnhnh_weight_file << current_time << " " ;
		}
		if((ring_id.compare("AA") == 0) || (ring_id.compare("DADDAA") == 0) ){
			if( print_path_id.compare("yes") == 0 ) { ringcoco_path_file << current_time << std::endl; }
			ringcoco_length_file << current_time << " " ;
			if( print_defect_id.compare("yes") == 0 ) { ringcoco_defect_file << current_time << " "; }
			ringcoco_weight_file << current_time << " " ;
		}
		if( print_defect_id.compare("yes") == 0 ) { defect_all_file << "#t=" << current_time << std::endl; }
		if( print_adjmatrix_id.compare("yes") == 0 ) { adjmatrix_file << "#t=" << current_time << std::endl; }
		
		//Calculating the scrambling time
		if (scramble_id.compare("yes") == 0) {
			int current_step;
                        current_step = (int)((double)(current_time-scrambling_ti)/(double)dt);
                        //Saving the path data of the current step into the scramble_length data
//                         if(current_step < scrambling_tot_frame){
//                             for(j=0; j<N_tot_path; ++j) {
//                                 ////std::cout << "check " << frame << " >> " << current_step << "frame " << frame << " current time " << current_time << "ti " << scrambling_ti << " dt " << dt << std::endl;
//                                 //scramble_len[frame][current_step][j]=old_length_data[current_step][j];
//                             }
//                         }
			////std::cout << "current time " << current_time << "current step " << current_step << std::endl;  
			for(i=0; i< current_step; ++i) {    //to get the scrambling time of wire before that time step
				//scrambling_check_num_file << frame	<< " " << i << ": ";
				////std::cout << " i " << i << " scrambling_time step " << current_step-i << " time " << current_time-i*dt << std::endl;
				for(j=0; j<N_tot_path; ++j) {        //exploring whole water wires
					//if( (i == scrambling_check_time) && (j == scrambling_check_path) ){
						//scrambling_check_node_file << frame << " " << i << " " << j << ": ";	
						//scrambling_check_conn_file << frame	<< " " << i << " " << j << ": ";
					//}
					if(old_length_data[i][j] == 0) { 
						//if( (i == scrambling_check_time) && (j == scrambling_check_path) ){
							//scrambling_check_node_file <<  std::endl;
						//}
						//scrambling_check_num_file << "0 ";
					} else {					
						earlycheck = false;
						double scramble_sum = 0;
						for(k=0; k < old_length_data[i][j]; ++k){  //getting old length data and path nodes
							int node1 = old_path_data[i][j][k];
							int node2 = old_path_data[i][j][k+1];
							//if( (i == scrambling_check_time) && (j == scrambling_check_path) ){
								//scrambling_check_node_file << old_path_data_origin[i][j][k] << " " ;
								//scrambling_check_conn_file << connection_matrix[node1][node2] << " ";
							//}
							scramble_sum += connection_matrix[node1][node2];
							////std::cout << frame << " frame " << i << " path " << j << " node " << k << " len " << old_length_data[i][j] << " con " << connection_matrix[node1][node2] << " scrambling " << scramble_t_bool[i][j] << std::endl;
							if ((connection_matrix[node1][node2] == 1 ) && ( scramble_t_bool[i][j] == true)){ 
								scramble_t[i][j] = current_time - i*dt; 
								earlycheck = true; 
							}
						}
						//scrambling_check_num_file << scramble_sum << " ";
						if( !earlycheck ) {scramble_t_bool[i][j] = false; }
						//if( (i == scrambling_check_time) && (j == scrambling_check_path) ){
							//scrambling_check_node_file << old_path_data_origin[i][j][old_length_data[i][j]] << " bool : " << scramble_t_bool[i][j] << " >time: " << scramble_t[i][j] << std::endl;
						//}
// 						scramble_len[frame][i][j]=scramble_sum;
					}
					//if( (i == scrambling_check_time) && (j == scrambling_check_path) ){ 
						//scrambling_check_conn_file <<  std::endl; 
					//}
				}
				//scrambling_check_num_file << std::endl;
			}
		}
		if((ring_id.compare( "DA" ) == 0) || (ring_id.compare("DD") == 0) || (ring_id.compare("AA") == 0)  || (ring_id.compare("DADDAA") == 0)){
                	graph.Restart(adj_matrix, N_matrix);
		}
		//Calculating k-shortest path length for all possible D-A pairs
		if((ring_id.compare( "DA" ) == 0) || (ring_id.compare( "DADDAA") == 0) ){
			for (unsigned int start_node = 0; start_node < n_NH; ++start_node){
                		end_node_vec=new int[n_CO];
				for (unsigned int end_node = 0; end_node < n_CO; ++end_node){
                                    end_node_vec[end_node]=n_OW+n_NH+end_node;
                                }
                                path_len=new int[n_CO]; // path length
                                int tmp = graph.dijkstra_all_ends(start_node, end_node_vec, n_CO, paths, path_len);
                                for (unsigned int end_node = 0; end_node < n_CO; ++end_node){
                                        if (print_path_id.compare("yes") == 0) {path_file << start_node+1 << ", " << end_node+1 << ": "; }
					////std::cout << "end_node " << end_node << "path size " << paths.size() << std::endl;
					int pp = (path_len[end_node] <= 0) ? 0 : path_len[end_node]-1;
					////std::cout << " end_node " << pp << std::endl;
					length_file << pp << " ";
					ring_node = new int[path_len[end_node]]; 
					node_ndx = 0;
					if (print_path_id.compare("yes") == 0) { path_file << pp << " >> "; }
					if(pp<=0)
					{ 
						if (print_path_id.compare("yes") == 0) { path_file << std::endl; }
						continue;
					}
					while(paths[end_node].size())
					{
						ring_node[node_ndx] = paths[end_node].top(); 
						node_ndx += 1;
						if (print_path_id.compare("yes") == 0) {
							if (paths[end_node].top() >= n_OW+n_NH) { 	
								if ( paths[end_node].top()-(n_OW+n_NH) != end_node ) { path_file << "(s)" ;}
								path_file << ndx.GetLocation(o_hb_grp, paths[end_node].top()-(n_OW+n_NH))+1 << " ";
							} //interal index starts from zero
							else if( paths[end_node].top() < n_NH ) {
								if ( paths[end_node].top() != start_node ) {path_file << "(s)";}
								path_file << ndx.GetLocation(h_hb_grp,paths[end_node].top())+1 << " ";  
							} //interal index starts from zero
							else { 
								path_file << paths[end_node].top()-n_NH+(last_pro_res_ndx+1) << " ";  
							} //matching matrix index with original residue number of water molecule for ploting
						}
						paths[end_node].pop();
					}	
					if ( direct_id.compare("undirected") == 0 ){
						n_dd = 0;
						//n_aa = 0;
						for ( m =1; m<node_ndx-1; ++m){
							if((connection_matrix[ring_node[m]][ring_node[m-1]] + connection_matrix[ring_node[m]][ring_node[m+1]]) == 2) {n_dd += 1;}
							//if((connection_matrix[ring_node[m-1]][ring_node[m]] + connection_matrix[ring_node[m+1]][ring_node[m]]) == 2) {n_aa += 1;}
						}
						//defect_file << n_dd << " " << n_aa << " ";
						defect_file << n_dd << " " ;
					}
					if ( weight_id.compare("weight") == 0 ){
						double tot_weight = 0;
						for ( m =0; m<node_ndx-1; ++m){
							tot_weight += 1.0 / adj_matrix[ring_node[m]][ring_node[m+1]];
						}
						weight_file << tot_weight << " ";
					}
					delete[] ring_node;
					////std::cout << "start " << start_node << " end " << end_node << " k_max " << k_max << " paths.size " << paths.size() << std::endl;
					if (print_path_id.compare("yes") == 0) { path_file << std::endl; }
				}
				paths.clear();
				delete[] path_len;
				delete[] end_node_vec;
			}
			length_file << std::endl;
			if ( direct_id.compare("undirected") == 0 ){defect_file << std::endl;}
			if ( weight_id.compare("weight") == 0 ){ weight_file << std::endl;}
		}
                
		if((ring_id.compare("DD") == 0) || (ring_id.compare("AA") == 0)  || (ring_id.compare("DADDAA") == 0)){
			if(direct_id.compare( "undirected") != 0){
				#pragma omp parallel private(i, j)
                                {
                                        #pragma omp for
                                        for(i=0; i < N_matrix; ++i){
                                                for(j =i+1; j < N_matrix; ++j)
						{
                                                        if( adj_matrix[i][j] < adj_matrix[j][i] ) 
							{ 
								adj_matrix[j][i] = adj_matrix[i][j]; 
							}
                                                        else if( adj_matrix[i][j] > adj_matrix[j][i] ) 
							{ 
								adj_matrix[j][i] = adj_matrix[i][j]; 
							}
                                                }
                                        }
                                }
                		graph.Restart(adj_matrix, N_matrix);
			}
		}
		//Calculating k-shortest path length for the "RING" connecting  all possible D-D pairs
                if((ring_id.compare( "DD" ) == 0) || (ring_id.compare( "DADDAA") == 0) ){
                        for (unsigned int start_node = 0; start_node < n_NH-1; ++start_node){
                                end_node_vec=new int[n_NH-start_node-1];
                                for (unsigned int end_node = 0; end_node < n_NH-start_node-1; ++end_node){
                                    end_node_vec[end_node]=start_node+end_node+1;
                                }
                                path_len=new int[n_NH-start_node-1]; // path length
                                int tmp = graph.dijkstra_all_ends(start_node, end_node_vec, n_NH-start_node-1, paths, path_len);
                                for (unsigned int end_node = 0; end_node < n_NH-start_node-1; ++end_node){
                                        if (print_path_id.compare("yes") == 0) {ringnhnh_path_file << start_node+1 << ", " << start_node+end_node+2 << ": "; }
                                        ////std::cout << "end_node " << end_node << "path size " << paths.size() << std::endl;
                                        int pp = (path_len[end_node] <= 0) ? 0 : path_len[end_node]-1;
                                        ////std::cout << " end_node " << pp << std::endl;
                                        ringnhnh_length_file << pp << " ";
                                        ring_node = new int[path_len[end_node]];
                                        node_ndx = 0;
                                        if (print_path_id.compare("yes") == 0) { ringnhnh_path_file << pp << " >> "; }
                                        if(pp<=0)
                                        {
                                                if (print_path_id.compare("yes") == 0) { ringnhnh_path_file << std::endl; }
                                                continue;
                                        }
                                        while(paths[end_node].size())
                                        {
                                                ring_node[node_ndx] = paths[end_node].top();
                                                node_ndx += 1;
                                                if (print_path_id.compare("yes") == 0) {
                                                        if (paths[end_node].top() >= n_OW+n_NH) {
                                                                ringnhnh_path_file << "(s)" ;
                                                                ringnhnh_path_file << ndx.GetLocation(o_hb_grp, paths[end_node].top()-(n_OW+n_NH))+1 << " ";
                                                        } //interal index starts from zero
                                                        else if( paths[end_node].top() < n_NH ) {
                                                                if ((paths[end_node].top() != start_node) && (paths[end_node].top() != start_node+end_node+1)) {ringnhnh_path_file << "(s)";}
                                                                ringnhnh_path_file << ndx.GetLocation(h_hb_grp,paths[end_node].top())+1 << " ";
                                                        } //interal index starts from zero
                                                        else {
                                                                ringnhnh_path_file << paths[end_node].top()-n_NH+(last_pro_res_ndx+1) << " ";
                                                        } //matching matrix index with original residue number of water molecule for ploting
                                                }
						paths[end_node].pop();
                                        }
					if ( direct_id.compare("undirected") == 0 ){
                                                n_dd = 0;
                                                //n_aa = 0;
                                                for ( m =1; m<node_ndx-1; ++m){
                                                        if((connection_matrix[ring_node[m]][ring_node[m-1]] + connection_matrix[ring_node[m]][ring_node[m+1]]) == 2) {n_dd += 1;}
                                                        //if((connection_matrix[ring_node[m-1]][ring_node[m]] + connection_matrix[ring_node[m+1]][ring_node[m]]) == 2) {n_aa += 1;}
                                                }
                                                //ringnhnh_defect_file << n_dd << " " << n_aa << " ";
                                                ringnhnh_defect_file << n_dd << " " ;
                                        }
                                        if ( weight_id.compare("weight") == 0 ){
                                                double tot_weight = 0;
                                                for ( m =0; m<node_ndx-1; ++m){
                                                        tot_weight += 1.0 / adj_matrix[ring_node[m]][ring_node[m+1]];
                                                }
                                                ringnhnh_weight_file << tot_weight << " ";
                                        }
                                        delete[] ring_node;
                                        ////std::cout << "start " << start_node << " end " << end_node << " k_max " << k_max << " paths.size " << paths.size() << std::endl;
                                        if (print_path_id.compare("yes") == 0) { ringnhnh_path_file << std::endl; }
                                }
                                paths.clear();
                                delete[] path_len;
                                delete[] end_node_vec;
                        }
                        ringnhnh_length_file << std::endl;
                        if ( direct_id.compare("undirected") == 0 ){ringnhnh_defect_file << std::endl;}
                        if ( weight_id.compare("weight") == 0 ){ ringnhnh_weight_file << std::endl;}
                }

		//Calculating k-shortest path length for the "RING" connecting all possible A-A pairs
                if((ring_id.compare( "AA" ) == 0) || (ring_id.compare( "DADDAA") == 0) ){
                        for (unsigned int start_node = 0; start_node < n_CO-1; ++start_node){
                                end_node_vec=new int[n_CO-start_node-1];
                                for (unsigned int end_node = 0; end_node < n_CO-start_node-1; ++end_node){
                                    end_node_vec[end_node]=n_OW+n_NH+start_node+end_node+1;
                                }
                                path_len=new int[n_CO-start_node-1]; // path length
                                int tmp = graph.dijkstra_all_ends(n_OW+n_NH+start_node, end_node_vec, n_CO-start_node-1, paths, path_len);
                                for (unsigned int end_node = 0; end_node < n_CO-start_node-1; ++end_node){
                                        if (print_path_id.compare("yes") == 0) {ringcoco_path_file << start_node+1 << ", " << start_node+end_node+2 << ": "; }
                                        ////std::cout << "end_node " << end_node << "path size " << paths.size() << std::endl;
                                        int pp = (path_len[end_node] <= 0) ? 0 : path_len[end_node]-1;
                                        ////std::cout << " end_node " << pp << std::endl;
                                        ringcoco_length_file << pp << " ";
                                        ring_node = new int[path_len[end_node]];
                                        node_ndx = 0;
                                        if (print_path_id.compare("yes") == 0) { ringcoco_path_file << pp << " >> "; }
                                        if(pp<=0)
                                        {
                                                if (print_path_id.compare("yes") == 0) { ringcoco_path_file << std::endl; }
                                                continue;
                                        }
                                        while(paths[end_node].size())
                                        {
                                                ring_node[node_ndx] = paths[end_node].top();
                                                node_ndx += 1;
                                                if (print_path_id.compare("yes") == 0) {
                                                        if (paths[end_node].top() >= n_OW+n_NH) {
                                                                if ((paths[end_node].top()-(n_OW+n_NH) != start_node) && (paths[end_node].top()-(n_OW+n_NH) != start_node+end_node+1) ) { ringcoco_path_file << "(s)" ;}
                                                                ringcoco_path_file << ndx.GetLocation(o_hb_grp, paths[end_node].top()-(n_OW+n_NH))+1 << " ";
                                                        } //interal index starts from zero
                                                        else if( paths[end_node].top() < n_NH ) {
                                                                ringcoco_path_file << "(s)";
                                                                ringcoco_path_file << ndx.GetLocation(h_hb_grp,paths[end_node].top())+1 << " ";
                                                        } //interal index starts from zero
                                                        else {
                                                                ringcoco_path_file << paths[end_node].top()-n_NH+(last_pro_res_ndx+1) << " ";
                                                        } //matching matrix index with original residue number of water molecule for ploting
                                                }
						paths[end_node].pop();
					}
                                        if ( direct_id.compare("undirected") == 0 ){
                                                n_dd = 0;
                                                //n_aa = 0;
                                                for ( m =1; m<node_ndx-1; ++m){
                                                        if((connection_matrix[ring_node[m]][ring_node[m-1]] + connection_matrix[ring_node[m]][ring_node[m+1]]) == 2) {n_dd += 1;}
                                                        //if((connection_matrix[ring_node[m-1]][ring_node[m]] + connection_matrix[ring_node[m+1]][ring_node[m]]) == 2) {n_aa += 1;}
                                                }
                                                //ringcoco_defect_file << n_dd << " " << n_aa << " ";
                                                ringcoco_defect_file << n_dd << " " ;
                                        }
                                        if ( weight_id.compare("weight") == 0 ){
                                                double tot_weight = 0;
                                                for ( m =0; m<node_ndx-1; ++m){
                                                        tot_weight += 1.0 / adj_matrix[ring_node[m]][ring_node[m+1]];
                                                }
                                                ringcoco_weight_file << tot_weight << " ";
                                        }
                                        delete[] ring_node;
                                        ////std::cout << "start " << start_node << " end " << end_node << " k_max " << k_max << " paths.size " << paths.size() << std::endl;
                                        if (print_path_id.compare("yes") == 0) { ringcoco_path_file << std::endl; }
                                }
                                paths.clear();
                                delete[] path_len;
                                delete[] end_node_vec;
                        }
                        ringcoco_length_file << std::endl;
                        if ( direct_id.compare("undirected") == 0 ){ringcoco_defect_file << std::endl;}
                        if ( weight_id.compare("weight") == 0 ){ ringcoco_weight_file << std::endl;}
                }
		////////////////////////////////////////////////////////////////////////////////////////////
		//-----------Definition of hydration shell------------------------------------------------//
		////////////////////////////////////////////////////////////////////////////////////////////
		if(shell_id.compare("yes") == 0){
			for(i=0; i<n_OW; ++i){
				double hb_dist, hb_ang;
				earlycheck=false;
				first_shell[frame][i] = 0;
				first_shell_N[frame][i] = 0;
				first_shell_O[frame][i] = 0;
				first_shell_C[frame][i] = 0;
				for (j=0; j<n_CO; ++j){
					if(earlycheck) break;
					hb_dist = distance(xyz_CO.at(j),xyz_OW.at(i*WDIM),box);
					if(hb_dist < hb_dist_cut){
						for(m=1; m<3; ++m){
							hb_ang = bond_angle(xyz_OW.at(i*WDIM+m), xyz_OW.at(i*WDIM), xyz_CO.at(j), box);
							if(hb_ang < hb_ang_cut) { 
								first_shell[frame][i] = 1;
								first_shell_O[frame][i] = 1;
								earlycheck = true;
								break;
							}
						}
					}
				}
				if(!earlycheck){
					for (j=0; j<n_NH; ++j){
						hb_dist = distance(xyz_nNH.at(j), xyz_OW.at(i*WDIM), box);
						if(hb_dist < hb_dist_cut){
							hb_ang = bond_angle(xyz_hNH.at(j), xyz_nNH.at(j), xyz_OW.at(i*WDIM), box);
							if(hb_ang < hb_ang_cut) 
							{ 
								first_shell[frame][i] = 1;
								first_shell_N[frame][i] = 1;
								earlycheck = true;
								break;
							}
						}
					}
				}
				if(!earlycheck){
					for (j=0; j<n_C; ++j){
						if(distance(xyz_OW.at(i*WDIM), xyz_heavy_C.at(j), box) < r_shell_1_C ){ 
							first_shell[frame][i] = 1;
							first_shell_C[frame][i] = 1;
							break;
						}
					}
				}
			}
		}
                
	}
	

	if(shell_id.compare("yes") ==0){       
		//Initialization of radial dipole-dipole spatial correlation function
		double N_radial_dpl_acf = 60.0;
		double q_OW=-0.834;
		double q_HW=0.417;
		vector < double > radial_dpl_acf, radial_long_dpl_acf, r_dpl_acf;
		radial_dpl_acf.resize(N_radial_dpl_acf);    //F(r) : radial dipole-dipole spatial correlation function
		radial_long_dpl_acf.resize(N_radial_dpl_acf); //L(r) : longitudinal radial dipole-dipole spatial correlation function
		r_dpl_acf.resize(N_radial_dpl_acf);	      //vector_r
		vector<double> tmp_num_ow_hist, tmp_dot_ow_hist, tmp_long_dot_ow_hist;
		tmp_num_ow_hist.resize(N_radial_dpl_acf);	
		tmp_dot_ow_hist.resize(N_radial_dpl_acf);	
		tmp_long_dot_ow_hist.resize(N_radial_dpl_acf);	

		double delta_r_dpl_acf = max_box_size/(double)(2.0*N_radial_dpl_acf);
		dpl_file << "# r of radial dipole-dipole spatial correlation function " << std::endl;
		for(j=0; j<N_radial_dpl_acf; ++j)
		{
			r_dpl_acf[j] = delta_r_dpl_acf*(j+0.5);	
			radial_dpl_acf[j] = 0.0;
			radial_long_dpl_acf[j] = 0.0;
			dpl_file << r_dpl_acf[j] << " ";
		}
		dpl_file << std::endl;
		//Solvation Shell water dynamics
	        vector<vector <coordinates>> dpl_first_shell;
	        vector<vector <coordinates>> dpl_first_shell_N;
	        vector<vector <coordinates>> dpl_first_shell_O;
	        vector<vector <coordinates>> dpl_first_shell_C;
        	vector<vector <double>> hb_first_shell_p_p;
	        vector<vector <double>> hb_first_shell_p_w;
        	vector<vector <double>> hb_first_shell_w_w;

        	dpl_first_shell.resize(n_tot_frame);
        	dpl_first_shell_N.resize(n_tot_frame);
        	dpl_first_shell_O.resize(n_tot_frame);
        	dpl_first_shell_C.resize(n_tot_frame);
	        hb_first_shell_p_p.resize(n_tot_frame);
	        hb_first_shell_p_w.resize(n_tot_frame);
        	hb_first_shell_w_w.resize(n_tot_frame);
	        int ndx_first_shell[n_OW];
        	int ndx_first_shell_N[n_OW];
	        int ndx_first_shell_O[n_OW];
	        int ndx_first_shell_C[n_OW];
        	int num_first_shell = 0;
	        int num_first_shell_N = 0;
        	int num_first_shell_O = 0;
	        int num_first_shell_C = 0;
		vector<int> num_first_shell_count_OW;
		num_first_shell_count_OW.resize(n_OW);
        	std::cout << " >> Creating the water molecules indices in first hydration shell ... " <<std::endl;
		ofstream fn_ow_freq_shell;
		fn_ow_freq_shell.open ("fn_ow_freq_shell.txt");
		fn_ow_freq_shell << "#OW-resid freq_shell freq_shell_O freq_shell_N freq_C" << std::endl; 
	    	for(i=0; i<n_OW; ++i){
			int sum_first_shell = 0;
			int sum_first_shell_N = 0;
			int sum_first_shell_O = 0;
			int sum_first_shell_C = 0;
			for(j=0; j<n_tot_frame;++j){
				sum_first_shell += first_shell[j][i];
				sum_first_shell_N += first_shell_N[j][i];
				sum_first_shell_O += first_shell_O[j][i];
				sum_first_shell_C += first_shell_C[j][i];
			}
			num_first_shell_count_OW[i] = sum_first_shell;
			fn_ow_freq_shell << i+(last_pro_res_ndx+1) << " " << sum_first_shell << " " << sum_first_shell_O << " " << sum_first_shell_N << " " << sum_first_shell_C << std::endl;
			if( sum_first_shell >= pop_shell*n_tot_frame ) { ndx_first_shell[num_first_shell] = i ; num_first_shell += 1; }
			if( sum_first_shell_N >= pop_shell*n_tot_frame ) { ndx_first_shell_N[num_first_shell_N] = i ; num_first_shell_N += 1; }
			if( sum_first_shell_O >= pop_shell*n_tot_frame ) { ndx_first_shell_O[num_first_shell_O] = i ; num_first_shell_O += 1; }
			if( sum_first_shell_C >= pop_shell*n_tot_frame ) { ndx_first_shell_C[num_first_shell_C] = i ; num_first_shell_C += 1; }
		}
  		fn_ow_freq_shell.close();
		shell_ndx_file << "# survival percentage to be considered in the shell " << pop_shell*100 << " %" << std::endl;
		shell_ndx_file << "# first hydration shell : " << "N, O, or S >> " << r_shell_1_N_O << " : C >> " << r_shell_1_C << std::endl;
		shell_ndx_file << "# r_shell_bulk : beyond " << r_shell_bulk << std::endl;
		shell_ndx_file << "# The water in the hydration shell is selected in the original time. " << std::endl;
		shell_ndx_file << "# Number of water molecules in first shell:All,O,N,C. " << std::endl;
		shell_ndx_file << num_first_shell << " " << num_first_shell_O << " " << num_first_shell_N << " " << num_first_shell_C << std::endl << std::endl;
		shell_ndx_file << "# VMD : first shell water indices (residue index by \"resid\"):" << std::endl;
	        shell_ndx_file << "# >>>> ALL : first solvation shell indices : " << std::endl;
		for(i=0; i < num_first_shell; ++i) { shell_ndx_file << ndx_first_shell[i] + (last_pro_res_ndx+1) << " "; }
		shell_ndx_file << std::endl;
	        shell_ndx_file << "# >>>> O : first solvation shell indices : " << std::endl;
		for(i=0; i < num_first_shell_O; ++i) { shell_ndx_file << ndx_first_shell_O[i] + (last_pro_res_ndx+1) << " "; }
		shell_ndx_file << std::endl;
        	shell_ndx_file << "# >>>> N : first solvation shell indices : " << std::endl;
		for(i=0; i < num_first_shell_N; ++i) { shell_ndx_file << ndx_first_shell_N[i] + (last_pro_res_ndx+1) << " "; }
		shell_ndx_file << std::endl;
	        shell_ndx_file << "# >>>> C : first solvation shell indices : " << std::endl;
		for(i=0; i < num_first_shell_C; ++i) { shell_ndx_file << ndx_first_shell_C[i] + (last_pro_res_ndx+1) << " "; }
		shell_ndx_file << std::endl;
	
		shell_ndx_file << "#GROMACS : first hydration shell OW : " << std::endl;
	        shell_ndx_file << "[ first_shell ] " << std::endl;//+1 because internal index starts from 0
		for(i=0; i < num_first_shell; ++i) { shell_ndx_file << ndx.GetLocation(ow_grp, ndx_first_shell[i]*WDIM)+1 << " "; } 
		shell_ndx_file << std::endl;
        	shell_ndx_file << "[ first_shell_O ] " << std::endl;//+1 because internal index starts from 0
		for(i=0; i < num_first_shell_O; ++i) { shell_ndx_file << ndx.GetLocation(ow_grp, ndx_first_shell_O[i]*WDIM)+1 << " "; } 
		shell_ndx_file << std::endl;
	        shell_ndx_file << "[ first_shell_N ] " << std::endl;//+1 because internal index starts from 0
		for(i=0; i < num_first_shell_N; ++i) { shell_ndx_file << ndx.GetLocation(ow_grp, ndx_first_shell_N[i]*WDIM)+1 << " "; } 
		shell_ndx_file << std::endl;
        	shell_ndx_file << "[ first_shell_C ] " << std::endl;//+1 because internal index starts from 0
		for(i=0; i < num_first_shell_C; ++i) { shell_ndx_file << ndx.GetLocation(ow_grp, ndx_first_shell_C[i]*WDIM)+1 << " "; } 
		shell_ndx_file << std::endl;
        	std::cout << ">> Creating is finished ... " <<std::endl;
        
		std::cout << ">> Extracting dipole moments and protein and/or water H-bond for the water molecules of the first solvation shell ... " << std::endl;
		double hb_dist, hb_ang;
		//Dipole and Hbond Calculation
		hb_span_file << "#time, # of hydration shell, #max size, #min1, #min2,..." << std::endl;
		for(frame=0; frame<n_tot_frame; ++frame){ 
			//Getting information and trajectory files for each frame
			box = traj->GetBox(frame);
			current_time = traj->GetTime(frame);
			xyz_nNH = traj->GetXYZ(frame, n_hb_grp);
			xyz_hNH = traj->GetXYZ(frame, h_hb_grp);
			xyz_CO = traj->GetXYZ(frame, o_hb_grp);
			xyz_OW = traj->GetXYZ(frame, ow_grp);
			xyz_heavy_C = traj->GetXYZ(frame, c_grp);
			xyz_active = traj->GetXYZ(frame, active_grp);
			com_protein = center_of_mass(xyz_active, mass_protein);
			for(k=0; k<3; ++k)
			{
				dpl_protein[k] = 0.0;
				for(j=0; j<n_protein; ++j)
				{
					dpl_protein[k] += charge_protein[j] * xyz_active[j][k];
				}
			}
			dpl_protein /= magnitude(dpl_protein);
			for(j=0; j<N_radial_dpl_acf; ++j)
			{
				tmp_num_ow_hist[j] = 0.;
				tmp_dot_ow_hist[j] = 0.;
				tmp_long_dot_ow_hist[j] = 0.;
			}	
			for(j=0; j<n_OW; ++j)
			{
				double tmp_ow_dist = distance(xyz_OW.at(j*WDIM), com_protein, box);
				int n_hist_dpl_acf = floor( tmp_ow_dist/delta_r_dpl_acf);
				if(n_hist_dpl_acf >= N_radial_dpl_acf ) n_hist_dpl_acf = N_radial_dpl_acf - 1;
				tmp_num_ow_hist[n_hist_dpl_acf] += 1;
				//q_HW = 1, q_OW = -2 * q_HW
				//coordinates vec_ow_dpl = xyz_OW.at(j*WDIM+1) + xyz_OW.at(j*WDIM+2) - 2.0*xyz_OW.at(j*WDIM); 
				coordinates vec_ow_dpl = bond_vector(xyz_OW.at(j*WDIM+1),xyz_OW.at(j*WDIM),box)+bond_vector(xyz_OW.at(j*WDIM+2),xyz_OW.at(j*WDIM),box); 
				//coordinates vec_ow_dpl = bond_vector(xyz_OW.at(j*WDIM+1), xyz_OW.at(j*WDIM), box); 
				vec_ow_dpl /= magnitude(vec_ow_dpl);
				coordinates vec_rij = bond_vector(xyz_OW.at(j*WDIM), com_protein, box);
				vec_rij /= magnitude(vec_rij);
			 	tmp_dot_ow_hist[n_hist_dpl_acf] += dot(dpl_protein, vec_ow_dpl);
				tmp_long_dot_ow_hist[n_hist_dpl_acf] += dot(dpl_protein, vec_rij)*dot(vec_ow_dpl,vec_rij);
			}
			for(j=0; j<N_radial_dpl_acf; ++j)
			{
				if (tmp_num_ow_hist[j] != 0 )
				{ 
					radial_dpl_acf[j] += tmp_dot_ow_hist[j]/tmp_num_ow_hist[j];
					radial_long_dpl_acf[j] += tmp_long_dot_ow_hist[j]/tmp_num_ow_hist[j];
				}
			}	
			//Hydration Shell Analysis	
			dpl_first_shell[frame].resize(num_first_shell);
			dpl_first_shell_O[frame].resize(num_first_shell_O);
			dpl_first_shell_N[frame].resize(num_first_shell_N);
			dpl_first_shell_C[frame].resize(num_first_shell_C);
			hb_first_shell_p_p[frame].resize(n_NH*n_CO);
			hb_first_shell_p_w[frame].resize(num_first_shell*(n_NH+n_CO));
			hb_first_shell_w_w[frame].resize((num_first_shell-1)*num_first_shell);
                        //Check the dipole moments of first hydration shell 
			for(j=0; j<num_first_shell; ++j){
				dpl_first_shell[frame][j] = bond_vector(xyz_OW.at(ndx_first_shell[j]*WDIM+1), xyz_OW.at(ndx_first_shell[j]*WDIM), box); 
			}
			for(j=0; j<num_first_shell_O; ++j){
				dpl_first_shell_O[frame][j] = bond_vector(xyz_OW.at(ndx_first_shell_O[j]*WDIM+1), xyz_OW.at(ndx_first_shell_O[j]*WDIM), box); 
			}
			for(j=0; j<num_first_shell_N; ++j){
				dpl_first_shell_N[frame][j] = bond_vector(xyz_OW.at(ndx_first_shell_N[j]*WDIM+1), xyz_OW.at(ndx_first_shell_N[j]*WDIM), box);
			}
			for(j=0; j<num_first_shell_C; ++j){
				dpl_first_shell_C[frame][j] = bond_vector(xyz_OW.at(ndx_first_shell_C[j]*WDIM+1), xyz_OW.at(ndx_first_shell_C[j]*WDIM), box); 
			}
			//Hbond connectivity between protein NH and CO 
                        int count_hb_first_shell_p_p =0;
			for(j=0; j<n_NH; ++j){
				for(k=0; k<n_CO; ++k){
					hb_first_shell_p_p[frame][count_hb_first_shell_p_p] = 0; 
					hb_dist = distance(xyz_nNH.at(j),xyz_CO.at(k),box);
					if(hb_dist < hb_dist_cut){
						hb_ang = bond_angle(xyz_hNH.at(j), xyz_nNH.at(j), xyz_CO.at(k), box);
						if(hb_ang < hb_ang_cut) { hb_first_shell_p_p[frame][count_hb_first_shell_p_p] = 1; }
					}
					count_hb_first_shell_p_p +=1;
				}
			}
			//Hbond connectivity between protein and first Solvation shell
                        int count_hb_first_shell_p_w =0;
			//Hbonds between NH group and OW in first shell
			for(j=0; j<n_NH; ++j){
				for(k=0; k<num_first_shell; ++k){
					hb_first_shell_p_w[frame][count_hb_first_shell_p_w] = 0; 
					hb_dist = distance(xyz_nNH.at(j), xyz_OW.at(ndx_first_shell[k]*WDIM), box);
					if(hb_dist < hb_dist_cut){
						hb_ang = bond_angle(xyz_hNH.at(j), xyz_nNH.at(j), xyz_OW.at(ndx_first_shell[k]*WDIM), box);
						if(hb_ang < hb_ang_cut) { hb_first_shell_p_w[frame][count_hb_first_shell_p_w] = 1; }
					}
					count_hb_first_shell_p_w +=1;
				}
			}
			//Hbonds between CO group and OW in first shell
			for(j=0; j<n_CO; ++j){
				for(k=0; k<num_first_shell; ++k){
					hb_first_shell_p_w[frame][count_hb_first_shell_p_w] = 0; 
					hb_dist = distance(xyz_CO.at(j),xyz_OW.at(ndx_first_shell[k]*WDIM),box);
					if(hb_dist < hb_dist_cut) {
						for(m=1; m<3; ++m){
							hb_ang = bond_angle(xyz_OW.at(ndx_first_shell[k]*WDIM+m), xyz_OW.at(ndx_first_shell[k]*WDIM), xyz_CO.at(j), box);
							if(hb_ang < hb_ang_cut) { hb_first_shell_p_w[frame][count_hb_first_shell_p_w] = 1; break;}
						}
					}
					count_hb_first_shell_p_w +=1;
				}
			}
			//Hbonds between OWs in first shell
                        int count_hb_first_shell_w_w =0;
			for(j=0; j<num_first_shell; ++j){
				for(k=0; k<num_first_shell; ++k){
					if(k == j) continue;
					hb_first_shell_w_w[frame][count_hb_first_shell_w_w]=0;
					hb_dist = distance(xyz_OW.at(ndx_first_shell[j]*WDIM), xyz_OW.at(ndx_first_shell[k]*WDIM), box);
					if(hb_dist < hb_dist_cut) {
						for(m =1; m<3; ++m){
							hb_ang = bond_angle(xyz_OW.at(ndx_first_shell[j]*WDIM+m), xyz_OW.at(ndx_first_shell[j]*WDIM), xyz_OW.at(ndx_first_shell[k]*WDIM), box);
							if(hb_ang < hb_ang_cut) { 
								hb_first_shell_w_w[frame][count_hb_first_shell_w_w] = 1; 
								break;
							}
						}
					}
					count_hb_first_shell_w_w += 1;
				}
			}
			///////////////////////////////////////////////////////////
			//--------------------Spanning Network Analysis-----------
			///////////////////////////////////////////////////////////
			//Hbonds between OWs in the instantaneous first hydration shell
			double **ow_hb_first_array = new double*[num_first_shell]; //HB network adjmatrix between first shell water molecules
			for(j=0; j<num_first_shell; ++j){
			 	ow_hb_first_array[j] = new double[num_first_shell];
				for(k=0; k<num_first_shell; ++k){
					ow_hb_first_array[j][k] = INF;
				}
			}
			for(j=0; j<num_first_shell; ++j){
				for(k=0; k<num_first_shell; ++k){
					if(k == j) continue;
					hb_dist = distance(xyz_OW.at(ndx_first_shell[j]*WDIM), xyz_OW.at(ndx_first_shell[k]*WDIM), box);
					if(hb_dist < hb_dist_cut) {
						for(m =1; m<3; ++m){
							hb_ang = bond_angle(xyz_OW.at(ndx_first_shell[j]*WDIM+m), xyz_OW.at(ndx_first_shell[j]*WDIM), xyz_OW.at(ndx_first_shell[k]*WDIM), box);
							if(hb_ang < hb_ang_cut) {
								ow_hb_first_array[j][k] = 1; 
								ow_hb_first_array[k][j] = 1; 
								break;
							}
						}
					}
				}
			}
			graph.Restart(ow_hb_first_array, num_first_shell);
			vector<int> ndx_span_tree=graph.prims_MST();
			vector<int> hist_span_tree(num_first_shell, 0);
			for(j=0; j<num_first_shell; ++j){
				hist_span_tree[ndx_span_tree[j]] += 1;
			}
			std::stable_sort(hist_span_tree.begin(), hist_span_tree.end());
			hb_span_file << current_time << " " << num_first_shell << " " ;
			//printing out from the last of the ascending order vector
			for(j=num_first_shell-1; j >= 0; j-=1){
				if(hist_span_tree[j] <=  0) break; 
				hb_span_file << hist_span_tree[j] << " ";
			}
			hb_span_file << endl;
			//----------------------------------------------------------
			//Orientational Tetrahedral order paprameter and LSI
			//-------------------------------------------------------
			//Initialization of distance matrix between the water molecules
			for(j=0; j<n_OW; ++j){
				ow_dist_matrix[j][j] = 0;
				for(k=j+1; k<n_OW; ++k){
					ow_dist_matrix[j][k] = distance(xyz_OW.at(j*WDIM), xyz_OW.at(k*WDIM), box);
					ow_dist_matrix[k][j] = ow_dist_matrix[j][k];
				}
			}
			//std::cout << " << Starting Orientational Tetrahedral Oder q and LSI " << std::endl;
			shell_q_file << current_time << " ";
			shell_lsi_file << current_time << " ";
			double r_nn_first;
			vector<int> ndx_nn_first;
			double_vec_t r_sort_delta;
			ndx_nn_first.resize(4);
			r_sort_delta.resize(n_OW);
			for(j=0; j<num_first_shell; ++j){
				vector<double> r_sort_all;
				r_sort_all.resize(n_OW);
				for(k=0; k<n_OW; ++k){
					r_sort_all[k] = ow_dist_matrix[ndx_first_shell[j]][k];
				}
				std::stable_sort(r_sort_all.begin(), r_sort_all.end());
				r_nn_first = r_sort_all[4]; //fourth largest distance
				//index of 4 nn
				int size_nn_tmp = 0;
				for(k =0; k<n_OW; ++k){
					if(k == ndx_first_shell[j]) continue;
					if(ow_dist_matrix[ndx_first_shell[j]][k] <= r_nn_first+0.000001) {
						if(size_nn_tmp > 3) break;
						ndx_nn_first[size_nn_tmp]=k;
						size_nn_tmp += 1;
					}
				}
				//q calculation
				double q_tmp = 0;
				for(k=0; k<3; ++k){
					for(m=k+1; m<4; ++m){
						double q_ang = bond_angle(xyz_OW.at(ndx_nn_first[k]*WDIM), xyz_OW.at(ndx_first_shell[j]*WDIM), xyz_OW.at(ndx_nn_first[m]*WDIM), box);
                                        	q_tmp += pow(cos(q_ang)+1.0/(double)3.0, 2);
					}
                                }
				q_tmp *= 3.0/(double)8.0;
				shell_q_file << 1.0 - q_tmp << " ";
				//LSI calculation
				double r_sort_delta_ave = 0;
				double r_sort_delta_ave2 = 0;
				int r_sort_num = 0;
				for (k=0; k<n_OW; ++k){
					if(r_sort_all[k] > 0.37 ) break; 
					if(k>0) 
					{
						r_sort_num = k; 
						r_sort_delta[k-1] = r_sort_all[k] - r_sort_all[k-1]; 
						r_sort_delta_ave += r_sort_delta[k-1]; 
					}
				}
				if(r_sort_num > 0) 
				{
					r_sort_delta_ave /= (double) r_sort_num;
					for (k=0; k<r_sort_num; ++k){
						r_sort_delta_ave2 += pow(r_sort_delta[k]-r_sort_delta_ave, 2);
					}
					shell_lsi_file << r_sort_delta_ave2/(double) r_sort_num << " ";
				}
				else {
					shell_lsi_file << 0.0 << " ";
				}
			}	
			shell_q_file << std::endl;
			shell_lsi_file << std::endl;
			//std::cout << " << Ending Orientational Tetrahedral Oder q and LSI " << std::endl;
		}
                
		dpl_file << " # Radial diplole-dipole spatial correlation function " << std::endl;
		for(j=0; j<N_radial_dpl_acf; ++j)
		{
			radial_dpl_acf[j] /= (double) n_tot_frame;
		 	dpl_file << radial_dpl_acf[j] << " " ;
		}
		dpl_file << std::endl;
		dpl_file << " # Longitudinal diplole-dipole spatial correlation function " << std::endl;
		for(j=0; j<N_radial_dpl_acf; ++j)
		{
			radial_long_dpl_acf[j] /= (double) n_tot_frame;
		 	dpl_file << radial_long_dpl_acf[j] << " " ;
		}
		dpl_file << std::endl;
		if(max_lag > n_tot_frame) {max_lag = n_tot_frame;}
                std::cout << ">> : Maximum lag (unit:number of frames) for correlation function: " << max_lag << std::endl;
		
		double_vec_t temp_hb;
		temp_hb.resize(n_tot_frame);
		vector<coordinates> temp_dpl;
		temp_dpl.resize(n_tot_frame);
                        
                double_vec_t temp_ACF, mean_ACF;
		temp_ACF.resize(max_lag);
                mean_ACF.resize(max_lag);
		//Time of ACF
		vector<double> t_ACF, exp_fit_coeff;
		t_ACF.resize(max_lag);
		exp_fit_coeff.resize(2);
                for(i=0; i< max_lag; ++i){ t_ACF[i] = dt*i;}
                //-----------------------------------------------------------
                //first and second order dipole auto correlation function
                //-----------------------------------------------------------        
                std::cout << ">> Starting 1st rank dipole-dipole correlation function ... " << std::endl;
		shell_acf_file << "#####-------1st dipole correlation function------------" << std::endl;
		shell_acf_file << "#1st dipole tau all" << std::endl;
                //First rank orientational correlation function for all water molecules in the first shell
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
                for(j=0; j<num_first_shell; ++j){
                    	for(i=0; i< n_tot_frame; ++i){ temp_dpl[i] = dpl_first_shell[i][j];}
             		temp_ACF = rot_first_ACF_func(temp_dpl, max_lag);
			if(temp_ACF[0] != 0 )
			{
				/*vector<double> tmp_non_negative_ACF;
				vector<double> tmp_non_negative_t;
	                        for(i=0; i< max_lag; ++i){ 
					if(temp_ACF[i] < 0) break;
					tmp_non_negative_ACF.push_back(temp_ACF[i]);
					tmp_non_negative_t.push_back(i*dt);
				}
				exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
				shell_acf_file << ndx_first_shell[j]+(last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell[j]] << " " << -1.0/exp_fit_coeff[1] << " ";*/
				shell_acf_file << ndx_first_shell[j]+(last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell[j]] << " " ;
	                        for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " ";}
				shell_acf_file << std::endl;
			}
                }
		shell_acf_file << std::endl;
                shell_acf_file << "#1st dipole averaged all " << std::endl;
                if(mean_ACF[0] !=0){
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
		shell_acf_file << std::endl;
		shell_acf_file << "#1st dipole tau O" << std::endl;
                //First rank orientational correlation function for all water molecules in the first shell near O
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
                for(j=0; j<num_first_shell_O; ++j){
                    	for(i=0; i< n_tot_frame; ++i){ temp_dpl[i] = dpl_first_shell_O[i][j];}
             		temp_ACF = rot_first_ACF_func(temp_dpl, max_lag);
			if(temp_ACF[0] != 0 )
			{
				/*vector<double> tmp_non_negative_ACF;
				vector<double> tmp_non_negative_t;
	                        for(i=0; i< max_lag; ++i){ 
					if(temp_ACF[i] < 0) break;
					tmp_non_negative_ACF.push_back(temp_ACF[i]);
					tmp_non_negative_t.push_back(i*dt);
				}
				exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
				shell_acf_file << ndx_first_shell_O[j]+(last_pro_res_ndx+1) << " " <<num_first_shell_count_OW[ndx_first_shell_O[j]] << " " << -1.0/exp_fit_coeff[1] << " ";*/
				shell_acf_file << ndx_first_shell_O[j]+(last_pro_res_ndx+1) << " " <<num_first_shell_count_OW[ndx_first_shell_O[j]] << " " ;
	                        for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " "; }
				shell_acf_file << std::endl;
			}
                }
                shell_acf_file << "#1st dipole averaged O " << std::endl;
                if(mean_ACF[0] !=0){
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
		shell_acf_file << std::endl;
		shell_acf_file << "#1st dipole tau N" << std::endl;
                //First rank orientational correlation function for all water molecules in the first shell near N
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
                for(j=0; j<num_first_shell_N; ++j){
                    	for(i=0; i< n_tot_frame; ++i){ temp_dpl[i] = dpl_first_shell_N[i][j];}
             		temp_ACF = rot_first_ACF_func(temp_dpl, max_lag);
			if(temp_ACF[0] != 0 )
			{
				/*vector<double> tmp_non_negative_ACF;
				vector<double> tmp_non_negative_t;
	                        for(i=0; i< max_lag; ++i){ 
					if(temp_ACF[i] < 0) break;
					tmp_non_negative_ACF.push_back(temp_ACF[i]);
					tmp_non_negative_t.push_back(i*dt);
				}
				exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
				shell_acf_file << ndx_first_shell_N[j] + (last_pro_res_ndx+1) << " " <<num_first_shell_count_OW[ndx_first_shell_N[j]] << " " << -1.0/exp_fit_coeff[1] << " ";*/
				shell_acf_file << ndx_first_shell_N[j] + (last_pro_res_ndx+1) << " " <<num_first_shell_count_OW[ndx_first_shell_N[j]] << " " ;
	                        for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " "; }
				shell_acf_file << std::endl;
			}
                }
                shell_acf_file << "#1st dipole averaged N " << std::endl;
                if(mean_ACF[0] !=0){
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
		shell_acf_file << std::endl;
		shell_acf_file << "#1st dipole tau C" << std::endl;
                //First rank orientational correlation function for all water molecules in the first shell near C
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
                for(j=0; j<num_first_shell_C; ++j){
                    	for(i=0; i< n_tot_frame; ++i){ temp_dpl[i] = dpl_first_shell_C[i][j];}
             		temp_ACF = rot_first_ACF_func(temp_dpl, max_lag);
			if(temp_ACF[0] != 0 )
			{
				/*vector<double> tmp_non_negative_ACF;
				vector<double> tmp_non_negative_t;
	                        for(i=0; i< max_lag; ++i){ 
					if(temp_ACF[i] < 0) break;
					tmp_non_negative_ACF.push_back(temp_ACF[i]);
					tmp_non_negative_t.push_back(i*dt);
				}
				exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
				shell_acf_file << ndx_first_shell_C[j] + (last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell_C[j]] << " " <<-1.0/exp_fit_coeff[1] << " ";*/
				shell_acf_file << ndx_first_shell_C[j] + (last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell_C[j]] << " ";
	                        for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " "; }
				shell_acf_file << std::endl;
			}
                }
                shell_acf_file << "#1st dipole averaged C " << std::endl;
                if(mean_ACF[0] !=0){
                    for(i=0; i<max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                    for(i=0; i<max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
                std::cout << "<< The calculation of the first rank dipole-dipole correlation function is finished ... " << std::endl;
                    
                std::cout << " >> Starting 2nd rank dipole-dipole correlation function ... " << std::endl;
		shell_acf_file << "#####-------2nd dipole correlation function------------" << std::endl;
		shell_acf_file << "#2nd dipole tau all" << std::endl;
                //Second rank orientational correlation function for all water molecules in the first shell 
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
                for(j=0; j<num_first_shell; ++j){
                    	for(i=0; i< n_tot_frame; ++i){ temp_dpl[i] = dpl_first_shell[i][j];}
             		temp_ACF = rot_second_ACF_func(temp_dpl, max_lag);
			if(temp_ACF[0] != 0 )
			{
				/*vector<double> tmp_non_negative_ACF;
				vector<double> tmp_non_negative_t;
	                        for(i=0; i< max_lag; ++i){ 
					if(temp_ACF[i] < 0) break;
					tmp_non_negative_ACF.push_back(temp_ACF[i]);
					tmp_non_negative_t.push_back(i*dt);
				}
				exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
				shell_acf_file << ndx_first_shell[j]+(last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell[j]] << " " << -1.0/exp_fit_coeff[1] << " ";*/
				shell_acf_file << ndx_first_shell[j]+(last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell[j]] << " ";
	                        for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " "; }
				shell_acf_file << std::endl;
			}
                }
		shell_acf_file << std::endl;
                shell_acf_file << "#2nd dipole averaged all " << std::endl;
                if(mean_ACF[0] !=0){
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
		shell_acf_file << std::endl;
                //Second rank orientational correlation function for all water molecules in the first shell near O
		shell_acf_file << "#2nd dipole tau O" << std::endl;
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
                for(j=0; j<num_first_shell_O; ++j){
                    	for(i=0; i< n_tot_frame; ++i){ temp_dpl[i] = dpl_first_shell_O[i][j];}
             		temp_ACF = rot_second_ACF_func(temp_dpl, max_lag);
			if(temp_ACF[0] != 0 )
			{
				/*vector<double> tmp_non_negative_ACF;
				vector<double> tmp_non_negative_t;
	                        for(i=0; i< max_lag; ++i){ 
					if(temp_ACF[i] < 0) break;
					tmp_non_negative_ACF.push_back(temp_ACF[i]);
					tmp_non_negative_t.push_back(i*dt);
				}
				exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
				shell_acf_file << ndx_first_shell_O[j]+(last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell_O[j]] << " " << -1.0/exp_fit_coeff[1] << " ";*/
				shell_acf_file << ndx_first_shell_O[j]+(last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell_O[j]] << " " ;
	                        for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " "; }
				shell_acf_file << std::endl;
			}
                }
                shell_acf_file << "#2nd dipole averaged O" << std::endl;
                if(mean_ACF[0] !=0){
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
		shell_acf_file << std::endl;
                //Second rank orientational correlation function for all water molecules in the first shell near N
		shell_acf_file << "#2nd dipole tau N" << std::endl;
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
                for(j=0; j<num_first_shell_N; ++j){
                    	for(i=0; i< n_tot_frame; ++i){ temp_dpl[i] = dpl_first_shell_N[i][j];}
             		temp_ACF = rot_second_ACF_func(temp_dpl, max_lag);
			if(temp_ACF[0] != 0 )
			{
				/*vector<double> tmp_non_negative_ACF;
				vector<double> tmp_non_negative_t;
	                        for(i=0; i< max_lag; ++i){ 
					if(temp_ACF[i] < 0) break;
					tmp_non_negative_ACF.push_back(temp_ACF[i]);
					tmp_non_negative_t.push_back(i*dt);
				}
				exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
				shell_acf_file << ndx_first_shell_N[j]+(last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell_N[j]] << " " << -1.0/exp_fit_coeff[1] << " ";*/
				shell_acf_file << ndx_first_shell_N[j]+(last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell_N[j]] << " " ;
	                        for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " "; }
				shell_acf_file << std::endl;
			}
                }
                shell_acf_file << "#2nd dipole averaged N " << std::endl;
                if(mean_ACF[0] !=0){
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
		shell_acf_file << std::endl;
                //Second rank orientational correlation function for all water molecules in the first shell near C
		shell_acf_file << "#2nd dipole tau C" << std::endl;
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
                for(j=0; j<num_first_shell_C; ++j){
                    	for(i=0; i< n_tot_frame; ++i){ temp_dpl[i] = dpl_first_shell_C[i][j];}
             		temp_ACF = rot_second_ACF_func(temp_dpl, max_lag);
			if(temp_ACF[0] != 0 )
			{
				/*vector<double> tmp_non_negative_ACF;
				vector<double> tmp_non_negative_t;
	                        for(i=0; i< max_lag; ++i){ 
					if(temp_ACF[i] < 0) break;
					tmp_non_negative_ACF.push_back(temp_ACF[i]);
					tmp_non_negative_t.push_back(i*dt);
				}
				exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
				shell_acf_file << ndx_first_shell_C[j]+(last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell_C[j]] << " " << -1.0/exp_fit_coeff[1] << " ";*/
				shell_acf_file << ndx_first_shell_C[j]+(last_pro_res_ndx+1) << " " << num_first_shell_count_OW[ndx_first_shell_C[j]] << " " ;
	                        for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " "; }
                		shell_acf_file << std::endl;
			}
                }
                shell_acf_file << "#2nd dipole averaged C" << std::endl;
                if(mean_ACF[0] !=0){
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
                std::cout << "<< The calculation of the second rank dipole-dipole correlation function is finished ... " << std::endl;

                //------------------------------------------------------------------------
                //Hbond correlation function between the protein NH and CO 
                //------------------------------------------------------------------------
                std::cout << ">> Starting the intermittent H-bond correlation function between protein NH and CO... " << std::endl;
                //H-bond intermittent correlation function P-P : <h(0)h(t)>
                shell_acf_file << "#tau : <h(0)h(t)> protein---Protein" << std::endl;
                for(i=0; i<max_lag; ++i){ mean_ACF[i] = 0;}
		int count_hb_first_shell_p_p=0;
                for(j=0; j<n_NH; ++j){
                	for(k=0; k<n_CO; ++k){
                        	for(i=0; i<n_tot_frame; ++i){ temp_hb[i] = hb_first_shell_p_p[i][count_hb_first_shell_p_p]; }
	         		temp_ACF = gen_ACF_func(temp_hb,max_lag);
				if(temp_ACF[0] != 0 )
				{
					/*vector<double> tmp_non_negative_ACF;
					vector<double> tmp_non_negative_t;
		                        for(i=0; i< max_lag; ++i){ 
						if(temp_ACF[i] < 0) break;
						tmp_non_negative_ACF.push_back(temp_ACF[i]);
						tmp_non_negative_t.push_back(i*dt);
					}
					exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
					shell_acf_file << "N " << ndx.GetLocation(n_hb_grp,j)+1 << " O " << ndx.GetLocation(o_hb_grp,k)+1 << " " << -1.0/exp_fit_coeff[1] << " ";*/
					shell_acf_file << "N " << ndx.GetLocation(n_hb_grp,j)+1 << " O " << ndx.GetLocation(o_hb_grp,k)+1 << " " ;
		                        for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " "; }
					shell_acf_file << std::endl;
				}
				count_hb_first_shell_p_p +=1;
			}
                }
                shell_acf_file << "#averaged : <h(0)h(t)> protein---Protein" << std::endl;
                if(mean_ACF[0] !=0){
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
                std::cout << "<< The calculation of the intermittent H-bond correlation function between protein NH and CO is finished ... " << std::endl;
                std::cout << ">> Starting the continuous H-bond correlation function between protein NH and CO ... " << std::endl;
                //H-bond continuous correlation function P-P : <h(0)H(t)>
                shell_acf_file << "#tau : <h(0)H(t)> protein---Protein" << std::endl;
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
		count_hb_first_shell_p_p=0;
                for(j=0; j<n_NH; ++j){
                	for(k=0; k<n_CO; ++k){
                        	for(i=0; i<n_tot_frame; ++i){ temp_hb[i] = hb_first_shell_p_p[i][count_hb_first_shell_p_p]; }
         			temp_ACF = hb_continuous_ACF_func(temp_hb,max_lag);
				if(temp_ACF[0] != 0 )
				{
					/*vector<double> tmp_non_negative_ACF;
					vector<double> tmp_non_negative_t;
	                        	for(i=0; i<max_lag; ++i){ 
						if(temp_ACF[i] < 0) break;
						tmp_non_negative_ACF.push_back(temp_ACF[i]);
						tmp_non_negative_t.push_back(i*dt);
					}
					exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
					shell_acf_file << "N " << ndx.GetLocation(n_hb_grp,j)+1 << " O " << ndx.GetLocation(o_hb_grp,k)+1 << " " << -1.0/exp_fit_coeff[1] << " ";*/
					shell_acf_file << "N " << ndx.GetLocation(n_hb_grp,j)+1 << " O " << ndx.GetLocation(o_hb_grp,k)+1 << " " ;
		                        for(i=0; i<max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " "; }
					shell_acf_file << std::endl;
				}
				count_hb_first_shell_p_p +=1;
			}
                }
                shell_acf_file << "#averaged :  <h(0)H(t)> protein---Protein" << std::endl;
                if(mean_ACF[0] !=0){
                      for(i=0; i<max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                      for(i=0; i<max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
                std::cout << "<< The calculation of the continuous H-bond correlation function between protein NH and CO is finished ... " << std::endl;

                //------------------------------------------------------------------------
                //Hbond correlation function between the protein and the first shell water
                //------------------------------------------------------------------------
                std::cout << ">> Starting the intermittent H-bond correlation function between protein and water of first shell... " << std::endl;
                //H-bond intermittent correlation function of protein and water in the first shell <h(0)h(t)>
                shell_acf_file << "#tau NH : <h(0)h(t)> protein---first shell water" << std::endl;
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
		int count_hb_first_shell_p_w =0;
		//Hbonds between NH group and OW in first shell
		for(j=0; j<n_NH; ++j){
			for(k=0; k<num_first_shell; ++k){
                        	for(i=0; i< n_tot_frame; ++i){ temp_hb[i] = hb_first_shell_p_w[i][count_hb_first_shell_p_w]; }
	         		temp_ACF = gen_ACF_func(temp_hb,max_lag);
				if(temp_ACF[0] != 0 )
				{
					/*vector<double> tmp_non_negative_ACF;
					vector<double> tmp_non_negative_t;
		                        for(i=0; i< max_lag; ++i){ 
						if(temp_ACF[i] < 0) break;
						tmp_non_negative_ACF.push_back(temp_ACF[i]);
						tmp_non_negative_t.push_back(i*dt);
					}
					exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
					shell_acf_file << "N " << ndx.GetLocation(n_hb_grp,j)+1 << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " " << -1.0/exp_fit_coeff[1] << " ";*/
					shell_acf_file << "N " << ndx.GetLocation(n_hb_grp,j)+1 << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " " ;
	        	                for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " ";}
					shell_acf_file << std::endl;
				}
				count_hb_first_shell_p_w +=1;
			}
                }
                shell_acf_file << "#averaged NH : <h(0)h(t)> protein---first shell water" << std::endl;
                if(mean_ACF[0] !=0){
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
		vector<double> mean_co_hh_ACF;
		mean_co_hh_ACF.resize(max_lag);
                shell_acf_file << "#tau CO : <h(0)h(t)> protein---first shell water" << std::endl;
                for(i=0; i< max_lag; ++i){ mean_co_hh_ACF[i] = 0;}
		//Hbonds between CO group and OW in first shell
		for(j=0; j<n_CO; ++j){
			for(k=0; k<num_first_shell; ++k){
                        	for(i=0; i< n_tot_frame; ++i){ temp_hb[i] = hb_first_shell_p_w[i][count_hb_first_shell_p_w]; }
	         		temp_ACF = gen_ACF_func(temp_hb,max_lag);
				if(temp_ACF[0] != 0 )
				{
					/*vector<double> tmp_non_negative_ACF;
					vector<double> tmp_non_negative_t;
		                        for(i=0; i< max_lag; ++i){ 
						if(temp_ACF[i] < 0) break;
						tmp_non_negative_ACF.push_back(temp_ACF[i]);
						tmp_non_negative_t.push_back(i*dt);
					}
					exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
					shell_acf_file << "O " << ndx.GetLocation(o_hb_grp,j)+1 << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " " << -1.0/exp_fit_coeff[1] << " ";*/
					shell_acf_file << "O " << ndx.GetLocation(o_hb_grp,j)+1 << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " " ;
	        	                for(i=0; i< max_lag; ++i)
					{ 
						mean_ACF[i] += temp_ACF[i]; 
						mean_co_hh_ACF[i] += temp_ACF[i];
						shell_acf_file << temp_ACF[i] << " ";
					}
					shell_acf_file << std::endl;
				}
				count_hb_first_shell_p_w +=1;
			}
                }
                shell_acf_file << "#averaged CO : <h(0)h(t)> protein---first shell water" << std::endl;
                if(mean_co_hh_ACF[0] !=0){
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_co_hh_ACF[i]/mean_co_hh_ACF[0] << " ";}
                }else{
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_co_hh_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
                shell_acf_file << "#averaged NH and CO : <h(0)h(t)> protein---first shell water" << std::endl;
                if(mean_ACF[0] !=0){
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
                std::cout << "<< The calculation of the intermittent H-bond correlation function between protein and water of first shell is finished ... " << std::endl;
                        
                std::cout << ">> Starting the continuous H-bond correlation function between protein and water of first shell ... " << std::endl;
                //H-bond continuous correlation function of protein and water in the first shell <h(0)H(t)>
                shell_acf_file << "#tau NH : <h(0)H(t)> protein---first shell water" << std::endl;
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
		count_hb_first_shell_p_w =0;
		//Hbonds between NH group and OW in first shell
		for(j=0; j<n_NH; ++j){
			for(k=0; k<num_first_shell; ++k){
                        	for(i=0; i< n_tot_frame; ++i){ temp_hb[i] = hb_first_shell_p_w[i][count_hb_first_shell_p_w]; }
         			temp_ACF = hb_continuous_ACF_func(temp_hb,max_lag);
				if(temp_ACF[0] != 0 )
				{
					/*vector<double> tmp_non_negative_ACF;
					vector<double> tmp_non_negative_t;
		                        for(i=0; i< max_lag; ++i){ 
						if(temp_ACF[i] < 0) break;
						tmp_non_negative_ACF.push_back(temp_ACF[i]);
						tmp_non_negative_t.push_back(i*dt);
					}
					exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
					shell_acf_file << "N " << ndx.GetLocation(n_hb_grp,j)+1 << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " " << -1.0/exp_fit_coeff[1] << " ";*/
					shell_acf_file << "N " << ndx.GetLocation(n_hb_grp,j)+1 << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " " ;
	        	                for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " ";}
					shell_acf_file << std::endl;
				}
				count_hb_first_shell_p_w +=1;
			}
                }
                shell_acf_file << "#averaged NH : <h(0)H(t)> protein---first shell water" << std::endl;
                if(mean_ACF[0] !=0){
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
		vector<double> mean_co_hH_ACF;
		mean_co_hH_ACF.resize(max_lag);
                shell_acf_file << "#tau CO : <h(0)H(t)> protein---first shell water" << std::endl;
                for(i=0; i< max_lag; ++i){ mean_co_hH_ACF[i] = 0;}
		//Hbonds between CO group and OW in first shell
		for(j=0; j<n_CO; ++j){
			for(k=0; k<num_first_shell; ++k){
                        	for(i=0; i< n_tot_frame; ++i){ temp_hb[i] = hb_first_shell_p_w[i][count_hb_first_shell_p_w]; }
         			temp_ACF = hb_continuous_ACF_func(temp_hb,max_lag);
				if(temp_ACF[0] != 0 )
				{
					/*vector<double> tmp_non_negative_ACF;
					vector<double> tmp_non_negative_t;
		                        for(i=0; i< max_lag; ++i){ 
						if(temp_ACF[i] < 0) break;
						tmp_non_negative_ACF.push_back(temp_ACF[i]);
						tmp_non_negative_t.push_back(i*dt);
					}
					exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
					shell_acf_file << "O " << ndx.GetLocation(o_hb_grp,j)+1 << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " " << -1.0/exp_fit_coeff[1] << " ";*/
					shell_acf_file << "O " << ndx.GetLocation(o_hb_grp,j)+1 << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " " ;
	        	                for(i=0; i< max_lag; ++i)
					{ 
						mean_ACF[i] += temp_ACF[i]; 
						mean_co_hH_ACF[i] += temp_ACF[i];
						shell_acf_file << temp_ACF[i] << " ";
					}
					shell_acf_file << std::endl;
				}
				count_hb_first_shell_p_w +=1;
			}
                }
                shell_acf_file << "#averaged CO : <h(0)H(t)> protein---first shell water" << std::endl;
                if(mean_co_hH_ACF[0] !=0){
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_co_hH_ACF[i]/mean_co_hH_ACF[0] << " ";}
                }else{
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_co_hH_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
                shell_acf_file << "#averaged NH and CO : <h(0)H(t)> protein---first shell water" << std::endl;
                if(mean_ACF[0] !=0){
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
                std::cout << "<< The calculation of the continuous H-bond correlation function between protein and water of first shell is finished ... " << std::endl;
                
		//-----------------------------------------------------------
		//Hbond correlation function of first shell water
		//---------------------------------------------------------
                std::cout << ">> Starting the intermittent H-bond correlation function between first shell water... " << std::endl;
                //H-bond intermittent correlation function of water in the first shell <h(0)h(t)>
                shell_acf_file << "#tau : <h(0)h(t)> first shell water---water" << std::endl;
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
		int count_hb_first_shell_w_w =0;
                for(j=0; j<num_first_shell; ++j){
	                for(k=0; k<num_first_shell; ++k){
				if(k == j) continue;
                        	for(i=0; i< n_tot_frame; ++i){ temp_hb[i] = hb_first_shell_w_w[i][count_hb_first_shell_w_w]; }
	         		temp_ACF = gen_ACF_func(temp_hb,max_lag);
				if(temp_ACF[0] != 0 )
				{
					/*vector<double> tmp_non_negative_ACF;
					vector<double> tmp_non_negative_t;
		                        for(i=0; i< max_lag; ++i){ 
						if(temp_ACF[i] < 0) break;
						tmp_non_negative_ACF.push_back(temp_ACF[i]);
						tmp_non_negative_t.push_back(i*dt);
					}
					exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
					shell_acf_file << "OW " << ndx_first_shell[j]+(last_pro_res_ndx+1) << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " " << -1.0/exp_fit_coeff[1] << " ";*/
					shell_acf_file << "OW " << ndx_first_shell[j]+(last_pro_res_ndx+1) << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " " ;
	        	                for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " "; }
					shell_acf_file << std::endl;
				}
				count_hb_first_shell_w_w += 1;
			}
                }
                shell_acf_file << "#averaged : <h(0)h(t)> first shell water---water" << std::endl;
                if(mean_ACF[0] !=0){
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                    for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
                std::cout << "<< The calculation of the intermittent H-bond correlation function between first shell water is finished ... " << std::endl;
                        
                std::cout << ">> Starting the continuous H-bond correlation function between first shell water ... " << std::endl;
                //H-bond continuous correlation function of water in the first shell <h(0)H(t)>
                shell_acf_file << "#tau : <h(0)H(t)> first shell water---water" << std::endl;
                for(i=0; i< max_lag; ++i){ mean_ACF[i] = 0;}
		count_hb_first_shell_w_w =0;
                for(j=0; j<num_first_shell; ++j){
	                for(k=0; k<num_first_shell; ++k){
				if(k == j) continue;
                        	for(i=0; i< n_tot_frame; ++i){ temp_hb[i] = hb_first_shell_w_w[i][count_hb_first_shell_w_w]; }
         			temp_ACF = hb_continuous_ACF_func(temp_hb,max_lag);
				if(temp_ACF[0] != 0 )
				{
					/*vector<double> tmp_non_negative_ACF;
					vector<double> tmp_non_negative_t;
	                        	for(i=0; i< max_lag; ++i){ 
						if(temp_ACF[i] < 0) break;
						tmp_non_negative_ACF.push_back(temp_ACF[i]);
						tmp_non_negative_t.push_back(i*dt);
					}
					exp_fit_coeff = exp_fit(tmp_non_negative_t, tmp_non_negative_ACF);
					shell_acf_file << "OW " << ndx_first_shell[j]+(last_pro_res_ndx+1) << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " " << -1.0/exp_fit_coeff[1] << " ";*/
					shell_acf_file << "OW " << ndx_first_shell[j]+(last_pro_res_ndx+1) << " OW " << ndx_first_shell[k]+(last_pro_res_ndx+1) << " ";
		                        for(i=0; i< max_lag; ++i){ mean_ACF[i] += temp_ACF[i]; shell_acf_file << temp_ACF[i] << " "; }
					shell_acf_file << std::endl;
				}
				count_hb_first_shell_w_w += 1;
			}
                }
                shell_acf_file << "#averaged : <h(0)H(t)> first shell water---water" << std::endl;
                if(mean_ACF[0] !=0){
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i]/mean_ACF[0] << " ";}
                }else{
                      for(i=0; i< max_lag; ++i){ shell_acf_file << mean_ACF[i] << " ";}
                }
                shell_acf_file << std::endl;
                std::cout << "<< The calculation of the continuous H-bond correlation function between first shell water is finished ... " << std::endl;
		
		//------------------------------------------------
		//survival probability in the first shell
		//---------------------------------------------
                std::cout << ">> Starting H-bond survival probability function ... " << std::endl;
                //Solvation shell survival probability
                shell_acf_file << "# <N_w(t)> of first hydration shell " << std::endl;
                for (i=0;i<max_lag;i++) { //i means lag
                        double sum=0;
                        for (j=0;j<n_tot_frame-i;j++) {
                                for(ow_i=0; ow_i<n_OW; ++ow_i){
                                        if( first_shell[j][ow_i] != 0){
                                                double ttt = 1;
                                                for(int rrr=j; rrr<=j+i; ++rrr) { 
							if(first_shell[rrr][ow_i] == 0) 
							{
								ttt=0; 
								break;
							}
						}
                                  		sum += ttt;
                                        }
                                }
                        }
			temp_ACF[i]=sum/double(n_tot_frame-i);
                        shell_acf_file << sum/double(n_tot_frame-i) << " ";
                }
                shell_acf_file << std::endl ;
		exp_fit_coeff = exp_fit(t_ACF, temp_ACF);
                shell_acf_file << "# tau : <N_w(t)> " << std::endl ;
		shell_acf_file << -1.0/exp_fit_coeff[1] << " ";
                shell_acf_file << std::endl ;
                std::cout << "<< The calculation of H-bond survival probability function is finished ... " << std::endl;
		dpl_first_shell.clear();
	        dpl_first_shell_O.clear();
        	dpl_first_shell_N.clear();
		dpl_first_shell_C.clear();
        	hb_first_shell_p_p.clear();
		hb_first_shell_p_w.clear();
	        hb_first_shell_w_w.clear();
        }
        //	int max_scramble_len_change=10;

//	if(max_lag > n_tot_frame) {max_lag = n_tot_frame;}
//	std::cout << " Maximum lag (unit: number of frames) for correlation function for scrambling path length : " << max_lag << std::endl;
	//Printing out scrambling time 
	if(scramble_id.compare("yes")==0){
//                 int_vec_t scramble_len_change;
//                 double_vec_t temp_sc, temp_sc_ACF;
//                 temp_sc.resize(n_tot_frame);
//                 temp_sc_ACF.resize(max_lag);
//                 scramble_len_change.resize(max_scramble_len_change);
//                 for(j=0; j < scrambling_tot_frame; ++j){
// 			for(k=0; k<N_tot_path; ++k){
// 				scrambling_len_change_file << j*dt << " " << k << " ";
// 				scrambling_len_acf_file << j*dt << " " << k << " ";
// 				for(i=0; i< max_scramble_len_change; ++i){ scramble_len_change[i] = 0; }
// 				for(i=0; i< max_lag; ++i) { temp_sc_ACF[i] = 0;}
// 				if( old_length_data[j][k] != 0)
// 				{
// 					for(i=j; i< (n_tot_frame-1); ++i){
// 						int ttt = (int)abs(scramble_len[i+1][j][k]-scramble_len[i][j][k]);
//                                                 
//                                                 if(ttt >= max_scramble_len_change ) 
//                                                 { 
//                                                    // //std::cout << "frame " << i << "scramble " << j << "path " << k << ": " << scramble_len[i][j][k] << " : " << scramble_len[i+1][j][k] <<">len-change" << ttt << std::endl;
//                                                     //std::cout << "Maximum value of scrambling path length change is set as : " << max_scramble_len_change-1 <<std::endl;
//                                                     //std::cout << "Please increase this value in the code " << std::endl;
//                                                     ttt = max_scramble_len_change -1;
//                                                 }
// 						scramble_len_change[ttt] += 1;
// 					}
// 					for(i=0; i<n_tot_frame; ++i){
// 						temp_sc[i] = scramble_len[i][j][k];
// 					}
//  					temp_sc_ACF = gen_ACF_func(temp_sc, max_lag);
// 				}
// 				for(i=0; i< max_scramble_len_change; ++i) { scrambling_len_change_file << scramble_len_change[i] << " ";}
// 				for(i=0; i< max_lag; ++i) { scrambling_len_acf_file << temp_sc_ACF[i] << " ";}
// 				scrambling_len_change_file << std::endl;
//                                 scrambling_len_acf_file << std::endl;
// 			}
// 		}
// 		
		////std::cout << " dt " << dt;
		for(i=0; i < scrambling_tot_frame; ++i){ 
			scrambling_time_file << i*dt << " ";
			for(j=0; j<N_tot_path; ++j){
				scrambling_time_file << scramble_t[i][j] << " ";
			}
			scrambling_time_file << std::endl;
		}
	}
	
	
	// Initialization of FPT analysis
	coordinates norm_vec;
	int n_fpt_time_slice;
	int n_fpt_left_box_slice;
	int n_fpt_right_box_slice;
	int n_fpt_box_slice;
	coordinates box_for_fpt;
	double set_box_boundary = 0.15; //accounting for PBC which will make empty space
	int n_fpt_interface_grp =0;
	vector<vector <coordinates>> xyz_OW_all;
	vector<vector <coordinates>> xyz_fpt_interface_grp;
	if( fpt_id.compare("yes") == 0 )
	{
                if(fpt_time_width > n_tot_frame)
                {
                        std::cout << "Time width ( " << fpt_time_width << " ) is set to be larger than " << n_tot_frame << std::endl;
                        n_fpt_time_slice = 1;
                        fpt_time_width = n_tot_frame;
                        std::cout << "# of time windows for FPT " << n_fpt_time_slice << std::endl;
                }else
                {
                        n_fpt_time_slice = floor( n_tot_frame/fpt_time_width );
                        std::cout << "# of time windows for FPT " << n_fpt_time_slice << std::endl;
                }

                box = traj->GetBox(0);
                box_for_fpt[X] = box(X,X);
		box_for_fpt[Y] = box(Y,Y);
		box_for_fpt[Z] = box(Z,Z);
                if(fpt_jump_id.compare("x") == 0 || fpt_jump_id.compare("X") == 0 )
                {
                        norm_vec[X]=1;
                        norm_vec[Y]=0;
                        norm_vec[Z]=0;
                }
                if(fpt_jump_id.compare("y") == 0 || fpt_jump_id.compare("Y") == 0 )
                {
                        norm_vec[X]=0;
                        norm_vec[Y]=1;
                        norm_vec[Z]=0;
                }
                if(fpt_jump_id.compare("z") == 0 || fpt_jump_id.compare("Z") == 0 )
                {
                        norm_vec[X]=0;
                        norm_vec[Y]=0;
                        norm_vec[Z]=1;
                }
		std::cout << "FPT : norm vector " << norm_vec << std::endl;
		xyz_OW_all.resize(n_tot_frame);
		if(fpt_interface_grp.compare("no") == 0)
		{
			if(ring_id.compare("bulk") != 0)
			{
				std::cout << "Please input correct \"fpt_interface_grp\" in input file which is now set to \"no\" " << std::endl;
				std::cout << "!!! system is not bulk!!!" << std::endl;
				return 0;
			}
			else
			{
				n_fpt_interface_grp =0;
			}
		}
		else
		{
			n_fpt_interface_grp = traj->GetNAtoms(fpt_interface_grp);
			xyz_fpt_interface_grp.resize(n_tot_frame);
		}
	
		for(frame=0; frame<n_tot_frame; ++frame)
                {       
                        xyz_OW_all[frame] = traj->GetXYZ(frame, ow_grp);
                        if(fpt_id.compare("yes") == 0 && n_fpt_interface_grp != 0 )
                        {       
                                xyz_fpt_interface_grp[frame] = traj->GetXYZ(frame, fpt_interface_grp);
                        }
                }
		//std::cout << "FPT : n_fpt_interface_grp : " << n_fpt_interface_grp << " fpt_id " << fpt_id << std::endl;
		double fpt_interface_bound_min, fpt_interface_bound_max;
		double real_box_edge = dot(box_for_fpt, norm_vec);
		//std::cout << "FPT : real_box_edge " << real_box_edge << std::endl;
                for(int id_fpt_time_width=0; id_fpt_time_width < n_fpt_time_slice; ++id_fpt_time_width)
                {
                        if (n_fpt_interface_grp == 0)
                        {
                                fpt_interface_bound_min = 0.0;
                                fpt_interface_bound_max = real_box_edge;
                        }
                        else
                        {
                                int start_time = id_fpt_time_width * fpt_time_width;
                                fpt_interface_bound_min = 0;
                                fpt_interface_bound_max = 0;
				for(k=start_time; k<start_time+fpt_time_width; ++k)
                                {
                                        double xxx_min=INF;
                                        double xxx_max=0;
                                        for (j=0; j<n_fpt_interface_grp; ++j)
                                        {
                                                double xxx_xxx = dot(xyz_fpt_interface_grp[k].at(j),norm_vec);
                                                if(xxx_xxx < xxx_min) { xxx_min = xxx_xxx; }
                                                if(xxx_xxx > xxx_max) { xxx_max = xxx_xxx; }
                                        }
					//std::cout << "id_fpt_time_width " << id_fpt_time_width << " xxx_min " << xxx_min << " xxx_max " << xxx_max << std::endl;
					fpt_interface_bound_min += xxx_min;
					fpt_interface_bound_max += xxx_max;
                                }
                                fpt_interface_bound_min /= fpt_time_width;
                                fpt_interface_bound_max /= fpt_time_width;
                        }
			//std::cout << "FPT : id_fpt_time_width " << id_fpt_time_width << " fpt_interface_bound_min " << fpt_interface_bound_min << 
			//" fpt_interface_bound_max " << fpt_interface_bound_max << std::endl;
                        double fpt_left_size = fpt_interface_bound_min;
                        double fpt_right_size = real_box_edge - fpt_interface_bound_max;
			//std::cout << "FPT : fpt_left_size " << fpt_left_size << " fpt_right_size " << fpt_right_size << std::endl;

			//Left slice
			if ( fpt_left_size > set_box_boundary )
			{
                        	n_fpt_left_box_slice = floor((fpt_left_size - set_box_boundary)/fpt_box_width);
			}
			else
			{	
				n_fpt_left_box_slice = 0;
			}
			//std::cout << "FPT : n_fpt_left_box_slice " << n_fpt_left_box_slice << std::endl;

                        for(int id_fpt_box_width=0; id_fpt_box_width < n_fpt_left_box_slice; ++id_fpt_box_width)
                        {
                                int start_time = id_fpt_time_width * fpt_time_width;
                                vector<int> ndx_all;
                                for(j=0; j<n_OW; ++j)
                                {
                                        if(dot(xyz_OW_all[start_time].at(j*WDIM),norm_vec) >= fpt_interface_bound_min - (id_fpt_box_width + 0.5) * fpt_box_width - fpt_box_epsilon && \
                                        dot(xyz_OW_all[start_time].at(j*WDIM),norm_vec) <= fpt_interface_bound_min - (id_fpt_box_width + 0.5) * fpt_box_width + fpt_box_epsilon )
                                        {
                                                ndx_all.push_back(j);
                                        }
                                }
				fpt_info_file << id_fpt_time_width+1 << " " << id_fpt_box_width+1 << " " << fpt_interface_bound_min - (id_fpt_box_width + 0.5) * fpt_box_width << " "
				<< fpt_time_width << " " << fpt_box_epsilon << " " << fpt_box_width << " " << ndx_all.size() << " L " << std::endl;
				fpt_ndx_file << "[ " << id_fpt_time_width+1 << " L " << id_fpt_box_width+1 << " ] [ ";
                                for(j=0; j<ndx_all.size(); ++j)
                                {
                                        fpt_ndx_file << ndx.GetLocation(ow_grp, ndx_all[j]*WDIM)+1 << " ";
                                }
                                fpt_ndx_file << "]" << std::endl;
                                vector<double> fpt_time_all;
                                fpt_time_all.resize(ndx_all.size());
                                vector<int> fpt_time_all_sign;
                                fpt_time_all_sign.resize(ndx_all.size());
				
                                for(j=0; j<ndx_all.size(); ++j)
                                {       
                                        fpt_time_all_sign[j]=0;
                                        for(k=start_time; k<start_time + fpt_time_width; ++k)
                                        {
                                                if(dot(xyz_OW_all[k].at(ndx_all[j]*WDIM),norm_vec) > fpt_interface_bound_min - id_fpt_box_width * fpt_box_width)
                                                {
                                                        fpt_time_all[j] = (k - start_time) * dt;
                                                        fpt_time_all_sign[j] = +1;
                                                        break;
                                                }
                                                if(dot(xyz_OW_all[k].at(ndx_all[j]*WDIM),norm_vec) < fpt_interface_bound_min - (id_fpt_box_width+1) * fpt_box_width)
                                                {
                                                        fpt_time_all[j] = (k - start_time) * dt;
                                                        fpt_time_all_sign[j] = -1;
                                                        break;
                                                }
                                        }       
                                        if(fpt_time_all_sign[j] == 0) { fpt_time_all[j] = fpt_time_width * dt; }
                                }               
                                //fpt_t_file << "[ " << id_fpt_time_width+1 << " L " << id_fpt_box_width+1 << " ] [ ";
                                fpt_t_file << "[ ";
                                for(j=0; j<ndx_all.size(); ++j)
                                {
                                        fpt_t_file << fpt_time_all[j] << " ";
                                }               
                                fpt_t_file << "]" << std::endl;
                                //fpt_sign_file << "[ " << id_fpt_time_width+1 << " L " << id_fpt_box_width+1 << " ] [ ";     
                                fpt_sign_file << "[ ";
                                for(j=0; j<ndx_all.size(); ++j)
                                {               
                                        fpt_sign_file << fpt_time_all_sign[j] << " ";
                                }               
                                fpt_sign_file << "]" << std::endl;
                                                        
                                ndx_all.clear();        
                                fpt_time_all.clear();
                                fpt_time_all_sign.clear();
			}
		
			//Right side	
			if ( fpt_right_size > set_box_boundary )
			{
                        	n_fpt_right_box_slice = floor((fpt_right_size - set_box_boundary)/fpt_box_width);
			}
			else
			{	
				n_fpt_right_box_slice = 0;
			}
			//std::cout << "FPT : n_fpt_right_box_slice " << n_fpt_right_box_slice << std::endl;

                        for(int id_fpt_box_width=0; id_fpt_box_width < n_fpt_right_box_slice; ++id_fpt_box_width)
                        {
                                int start_time = id_fpt_time_width * fpt_time_width;
                                vector<int> ndx_all;
                                for(j=0; j<n_OW; ++j)
                                {
                                        if(dot(xyz_OW_all[start_time].at(j*WDIM),norm_vec) >= fpt_interface_bound_max + (id_fpt_box_width + 0.5) * fpt_box_width - fpt_box_epsilon && 
                                        dot(xyz_OW_all[start_time].at(j*WDIM),norm_vec) <= fpt_interface_bound_max + (id_fpt_box_width + 0.5) * fpt_box_width + fpt_box_epsilon )
                                        {
                                                ndx_all.push_back(j);
                                        }
                                }
				fpt_info_file << id_fpt_time_width+1 << " " << id_fpt_box_width+1 << " " << fpt_interface_bound_max + (id_fpt_box_width + 0.5) * fpt_box_width << " "
				<< fpt_time_width << " " << fpt_box_epsilon << " " << fpt_box_width << " " << ndx_all.size() << " R " << std::endl;
				fpt_ndx_file << "[ " << id_fpt_time_width+1 << " R " << id_fpt_box_width+1 << " ] [ ";
                                for(j=0; j<ndx_all.size(); ++j)
                                {
                                        fpt_ndx_file << ndx.GetLocation(ow_grp, ndx_all[j]*WDIM)+1 << " ";
                                }
                                fpt_ndx_file << "]" << std::endl;
                                vector<double> fpt_time_all;
                                fpt_time_all.resize(ndx_all.size());
                                vector<int> fpt_time_all_sign;
                                fpt_time_all_sign.resize(ndx_all.size());
				
                                for(j=0; j<ndx_all.size(); ++j)
                                {       
                                        fpt_time_all_sign[j]=0;
                                        for(k=start_time; k<start_time + fpt_time_width; ++k)
                                        {
                                                if(dot(xyz_OW_all[k].at(ndx_all[j]*WDIM),norm_vec) > fpt_interface_bound_max + (id_fpt_box_width+1) * fpt_box_width)
                                                {
                                                        fpt_time_all[j] = (k - start_time) * dt;
                                                        fpt_time_all_sign[j] = +1;
                                                        break;
                                                }
                                                if(dot(xyz_OW_all[k].at(ndx_all[j]*WDIM),norm_vec) < fpt_interface_bound_max + id_fpt_box_width * fpt_box_width)
                                                {
                                                        fpt_time_all[j] = (k - start_time) * dt;
                                                        fpt_time_all_sign[j] = -1;
                                                        break;
                                                }
                                        }       
                                        if(fpt_time_all_sign[j] == 0) { fpt_time_all[j] = fpt_time_width * dt; }
                                }               
                                //fpt_t_file << "[ " << id_fpt_time_width+1 << " R " << id_fpt_box_width+1 << " ] [ ";
                                fpt_t_file << "[ ";
                                for(j=0; j<ndx_all.size(); ++j)
                                {
                                        fpt_t_file << fpt_time_all[j] << " ";
                                }               
                                fpt_t_file << "]" << std::endl;
                                //fpt_sign_file << "[ " << id_fpt_time_width+1 << " R " << id_fpt_box_width+1 << " ] [ ";
                                fpt_sign_file << "[ ";
                                for(j=0; j<ndx_all.size(); ++j)
                                {               
                                        fpt_sign_file << fpt_time_all_sign[j] << " ";
                                }               
                                fpt_sign_file << "]" << std::endl;
                                                        
                                ndx_all.clear();        
                                fpt_time_all.clear();
                                fpt_time_all_sign.clear();
			}
			//If the system is bulk or no interface along given direction
			if ( n_fpt_left_box_slice == 0 && n_fpt_right_box_slice == 0 )
			{
				if ( real_box_edge > 2.0 * set_box_boundary )
				{
					n_fpt_box_slice = floor((real_box_edge - 2.0 * set_box_boundary)/fpt_box_width);
				}
				else
				{
					n_fpt_box_slice = 0;
					std::cout << "Error! The box size is smaller than 2.0 * set_box_boundary : Too small " << std::endl;
					return 0;
				}
				for(int id_fpt_box_width=0; id_fpt_box_width < n_fpt_box_slice; ++id_fpt_box_width)
                		{
					int start_time = id_fpt_time_width * fpt_time_width;
					vector<int> ndx_all;
					for(j=0; j<n_OW; ++j)
					{
						if(dot(xyz_OW_all[start_time].at(j*WDIM),norm_vec) >= set_box_boundary + (id_fpt_box_width + 0.5) * fpt_box_width - fpt_box_epsilon && 
						dot(xyz_OW_all[start_time].at(j*WDIM),norm_vec) <= set_box_boundary + (id_fpt_box_width + 0.5) * fpt_box_width + fpt_box_epsilon )
						{
							ndx_all.push_back(j);
						}
					}
					fpt_info_file << id_fpt_time_width+1 << " " << id_fpt_box_width+1 << " " << set_box_boundary + (id_fpt_box_width + 0.5) * fpt_box_width << " "
					<< fpt_time_width << " " << fpt_box_epsilon << " " << fpt_box_width << " " << ndx_all.size() << " B " << std::endl;
					fpt_ndx_file << "[ " << id_fpt_time_width+1 << " B " << id_fpt_box_width+1 << " ] [ ";
					for(j=0; j<ndx_all.size(); ++j)
					{
						fpt_ndx_file << ndx.GetLocation(ow_grp, ndx_all[j]*WDIM)+1 << " ";
					}
					fpt_ndx_file << "]" << std::endl;
					vector<double> fpt_time_all;
					fpt_time_all.resize(ndx_all.size());
					vector<int> fpt_time_all_sign;
					fpt_time_all_sign.resize(ndx_all.size());
					
					for(j=0; j<ndx_all.size(); ++j)
					{       
						fpt_time_all_sign[j]=0;
						for(k=start_time; k<start_time + fpt_time_width; ++k)
						{
							if(dot(xyz_OW_all[k].at(ndx_all[j]*WDIM),norm_vec) > set_box_boundary + (id_fpt_box_width+1) * fpt_box_width)
							{
								fpt_time_all[j] = (k - start_time) * dt;
								fpt_time_all_sign[j] = +1;
								break;
							}
							if(dot(xyz_OW_all[k].at(ndx_all[j]*WDIM),norm_vec) < set_box_boundary + id_fpt_box_width * fpt_box_width)
							{
								fpt_time_all[j] = (k - start_time) * dt;
								fpt_time_all_sign[j] = -1;
								break;
							}
						}       
						if(fpt_time_all_sign[j] == 0) { fpt_time_all[j] = fpt_time_width * dt; }
					}               
					fpt_t_file << "[ ";
					for(j=0; j<ndx_all.size(); ++j)
					{
						fpt_t_file << fpt_time_all[j] << " ";
					}               
					fpt_t_file << "]" << std::endl;
					fpt_sign_file << "[ ";
					for(j=0; j<ndx_all.size(); ++j)
					{               
						fpt_sign_file << fpt_time_all_sign[j] << " ";
					}               
					fpt_sign_file << "]" << std::endl;
								
					ndx_all.clear();        
					fpt_time_all.clear();
					fpt_time_all_sign.clear();
				}
			}
			
		}
		xyz_OW_all.clear();
	}

	//Deleting all variables
	for(i =0; i<N_matrix; ++i) {
        	delete[] adj_matrix[i]; 
			delete[] connection_matrix[i]; 
        }
	delete[] adj_matrix;
	delete[] connection_matrix;
        delete[] defect_in;
	delete[] defect_out;
        
	//if(shell_id.compare("yes") == 0){
        //    for(i =0; i< n_tot_frame; ++i){
                    //delete[] bulk[i];
                    // 		for(j=0; j< scrambling_tot_frame; ++j){
                    // 			delete[] scramble_len[i][j];
                    // 		}
                    // 		delete[] scramble_len[i];
        //    }
            //delete[] bulk;
        //}
        if(scramble_id.compare("yes") == 0){
            for(i=0; i< scrambling_tot_frame; ++i){
                    delete[] scramble_t[i];
                    delete[] scramble_t_bool[i];
                    delete[] old_length_data[i];
                    for(j=0; j< N_tot_path; ++j){
                            delete[] old_path_data[i][j];
                            delete[] old_path_data_origin[i][j];
                    }
                    delete[] old_path_data[i];
                    delete[] old_path_data_origin[i];
            }
            
            delete[] old_length_data;
            delete[] scramble_t;
            delete[] scramble_t_bool;
            delete[] old_path_data;
            delete[] old_path_data_origin;
            // 	delete[] scramble_len;
 
        }
	//Closing all files
	path_file.close();
	length_file.close();
	defect_file.close();
	weight_file.close();
	ringnhnh_path_file.close();
	ringnhnh_length_file.close();
	ringnhnh_defect_file.close();
	ringnhnh_weight_file.close();
	ringcoco_path_file.close();
	ringcoco_length_file.close();
	ringcoco_defect_file.close();
	ringcoco_weight_file.close();

	bulk_path_file1.close();
        bulk_el_force_file1.close();
        bulk_length_file1.close();
        bulk_defect_file1.close();
        bulk_weight_file1.close();
        bulk_path_file2.close();
        bulk_el_force_file2.close();
        bulk_length_file2.close();
        bulk_defect_file2.close();
        bulk_weight_file2.close();
        
	defect_all_file.close();
        adjmatrix_file.close();
	
	shell_num_file.close();
	shell_ndx_file.close();
	shell_inout_file.close();
	water_defect_file.close();
	shell_acf_file.close();
        shell_q_file.close();
	shell_lsi_file.close(); 
	dpl_file.close();
	hb_span_file.close();	
	
	length_old_file.close();
	path_old_file.close();
	
	scrambling_check_conn_file.close();
	scrambling_check_node_file.close();
	scrambling_check_num_file.close();
	scrambling_len_acf_file.close();
	scrambling_len_change_file.close();

	fpt_t_file.close();
	fpt_sign_file.close();
	fpt_info_file.close();
	fpt_ndx_file.close();
	
	std::cout << " Calculation FINISHED !"<<  " TOTAL TIME: "<< (double)(clock() - tStart)/CLOCKS_PER_SEC << std::endl;
	return 0;
}

};

int main(){
	HBNetwork test;
  	return test.main(stdin,stdout);
}
