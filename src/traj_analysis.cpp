//traj_an.cpp

/*********************************************************************/
/* This source code file is part of LQQCMMtools software project created 
 * by Igor Barden Grillo at Federal University of Paraíba. 
 * barden.igor@gmail.com ( Personal e-mail ) 
 * igor.grillo@acad.pucrs.br ( Academic e-mail )
 * quantum-chem.pro.br ( group site )
 * IgorChem ( Git Hub account )
 */ 

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
/*********************************************************************/

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <cmath>

#include "../include/traj_analysis.h"
#include "../include/global.h"
#include "../include/Line.h"
#include "../include/pdbAtom.h"
#include "../include/pdbModel.h"
#include "../include/PDB.h"
#include <experimental/filesystem>


using std::string;
using std::vector;
using std::cout;
using std::endl;
namespace fs = std::experimental::filesystem;

/***********************************************************************/
traj_an::traj_an(){
	
}
/***********************************************************************/
traj_an::traj_an(string file_name):
	traj_file(file_name)		{
		
	py_name = remove_extension( file_name.c_str() );
	py_name += "_mdtraj.py";
			
	python_script.open( py_name.c_str() );
	
	R_name = remove_extension( file_name.c_str() );
	R_name += "_rscript.R";
	
	R_script.open( R_name.c_str() );
}
/***********************************************************************/
traj_an::traj_an(string file_name, vector<int> atoms):
	traj_file(file_name)					   ,
	atoms_pairs(atoms)						   {
		
	py_name = remove_extension( file_name.c_str() );
	py_name += "_mdtraj.py";
	python_script.open( py_name.c_str() );
			
	R_name = remove_extension( file_name.c_str() );
	R_name += "_rscript.R";
	
	R_script.open( R_name.c_str() );
}
/***********************************************************************/
traj_an::~traj_an(){
}
/***********************************************************************/
void traj_an::mdtraj_geo(){
	
	string topname = change_extension( traj_file.c_str(),".pdb" );
	python_script << " #/usr/bin/env python \n"
				  << " # -*- coding: utf-8 -*- \n\n"
				  << "import mdtraj as md \n"
				  << "import os\n"
				  << "trj = md.load('" << traj_file << "',top='" << topname << "')\n"
				  << "rg = md.rmsd(trj,trj)\n"
				  << "rmsd = md.compute_rg(trj)\n"
				  << "file_txt = 'Time RMSD RG'+os.linesep\n"
				  << "file_rmsd = open('md_geoD.txt','w')\n"
				  << "for i in range(len(rmsd)):\n"
				  << "	file_txt += str(trj.time[i]) +' '+ str(rmsd[i]) +' '+str(rg[i]) +os.linesep\n"
				  << "file_rmsd.write(file_txt)\n"
				  << "file_rmsd.close()\n"
				  << "print('md_traj finished')";
				  
	
	python_script.close();
	
	string comand = "python3 " + py_name; 
	system( comand.c_str() );	

	vector<double> time;
	vector<double> rmsd;
	vector<double> rg;
	
	int line = 0;
	char tmp_line[50];
	double tmp_doub;
	int col = 0;
	string dat_name = "md_geoD.txt";
	std::ifstream geo_data( dat_name.c_str() );
	while( !geo_data.eof() ){
		geo_data.getline(tmp_line,50);
		if ( line > 0 ){
			std::stringstream stream(tmp_line);
			while( stream >> tmp_doub ){
				if ( col == 0 ){
					time.push_back(tmp_doub);
					col++;
				}else if ( col == 1 ){
					rmsd.push_back(tmp_doub);
					col++;
				}else if ( col == 2 ){
					rg.push_back(tmp_doub);
					col = 0;
				}
			}
		}
		line++;
	}
	
	int frame_mp = this->bi_most_probable_point(rmsd,rg);
	cout << "Frame number with most probable values of RMSD and RG: " << frame_mp << endl;
	cout << "RMSD = " << rmsd[frame_mp] << endl;
	cout << "RG = " << rg[frame_mp] << endl;

	R_script << "library(ggplot2)\n"
			 << "a <-read.table('md_geoD.txt',header=T)\n"
			 << "gp1 <-ggplot(a,aes(x=Time,y=RG))+\n"
			 << "theme_minimal()+ \n"
			 << "geom_line()+ \n"
			 << "ylab('RG')+\n"
			 << "xlab('Time (ps)')\n"
			 << "png('" << traj_file << "_RG.png',units='in',res=600,width=6,height=4) \n"
			 << "gp1\n"
			 << "dev.off()\n" 
			 << "a <-read.table('md_geoD.txt',header=T)\n"
			 << "gp2 <-ggplot(a,aes(x=Time,y=RMSD))+\n"
			 << "theme_minimal()+ \n"
			 << "geom_line()+ \n"
			 << "ylab('RMSD')+\n"
			 << "xlab('Time (ps)')\n"
			 << "png('" << traj_file << "_RMSD.png',units='in',res=600,width=6,height=4) \n"
			 << "gp2\n"
			 << "dev.off()\n"
			 << "aa <-data.frame(scale(a))\n"
			 << "gp3 <- ggplot(aa, aes(x=RG, y=RMSD)) +\n"
			 << "geom_point()+ \n"
			 << "geom_density_2d() + \n"
			 << "stat_density_2d(aes(fill = ..level..), geom='polygon')+\n"	
			 << "scale_fill_gradient(low='blue', high='red')+ \n"
			 << "ylab('RMSD')+\n"
			 << "xlab('RG') \n"
			 << "png('RMSD_rg_biplot.png',units='in',res=600,width=6,height=4)\n"
			 << "gp3 \n"
			 << "dev.off()";  
	
	R_script.close();
	
	comand = "Rscript " + R_name; 
	system( comand.c_str() );
}
/***********************************************************************/
void traj_an::calc_distances(const char* pdb_file){
	
	char tmp_line[80];
	double tmp_doub;
	int line = 0;
	
	
	PDB traj_pdb(pdb_file);
		
	vector<vector <double> > dist;
	dist.resize( atoms_pairs.size()/2 );
	
	for(int i=0;i<traj_pdb.nModels;i++){
		for ( int j=0;j<dist.size();j++){
			dist[j].emplace_back( traj_pdb.models[i].atom_distance( atoms_pairs[j], atoms_pairs[j+1] ) );
		}
	}

	string frame_name = remove_extension(pdb_file);

	frame_name += "_dists";
	std::ofstream dist_file;
	dist_file.open( frame_name.c_str() );
	
	for( int i=0;i<dist.size();i++){
		dist_file << "pair_" << std::to_string(i) << " ";
	}
	dist_file << endl;
	
	for( int j=0;j<dist[0].size();j++){
		for( int i=0;i<dist.size();i++){
			dist_file << dist[i][j] << " ";
		}
		dist_file << endl;
	}
	
	int frame = 0;
	if ( dist.size() == 2 ){
		frame = this->bi_most_probable_point(dist[0],dist[1]);	
		cout << "Frame number with most probable values of pair distances: " << frame << endl;
		cout << "D pair 1 = " << dist[0][frame] << endl;
		cout << "D pair 2 = " << dist[1][frame] << endl;
		cout << "Extracting frame and writing it in a PDB file!" << endl;
		this->extract_frame(pdb_file,frame);
	}	
	
	R_script << "b <-read.table('" <<frame_name<<"',header=T)\n"
			 << "gp4 <- ggplot(b, aes(x=pair_0, y=pair_1)) + \n"
			 << "geom_point()+ \n"
			 << "geom_density_2d() + \n"
			 << "stat_density_2d(aes(fill = ..level..), geom='polygon')+ \n"
			 << "scale_fill_gradient(low='blue', high='red')+ \n"
			 << "ylab('Pair 1')+ \n"
			 << "xlab('Pair 2') \n"
			 << "png('"<<frame_name<<"_biplotDist.png',units='in',res=600,width=6,height=4)\n"
			 << "gp4\n"
			 << "dev.off()";

	dist_file.close();
	
	python_script.close();
	R_script.close();
	
}
/***********************************************************************/
void traj_an::extract_frame(const char* pdb_file, int frames){
	
	PDB traj_pdb(pdb_file);
	string frame = std::to_string(frames);
	string frame_name = remove_extension(pdb_file);
	frame_name += "_#";
	frame_name += frame;
	frame_name += ".pdb";
	
	traj_pdb.models[frames].write_model(frame_name);
}
/***********************************************************************/
void traj_an::extract_frames(const char* pdb_file, int interval,int fr_sz){
	for( int i=1;i<=fr_sz;i++){
		if ( i%interval == 0 ){
			this->extract_frame(pdb_file,i);
		}
	}
}
/***********************************************************************/
int traj_an::bi_most_probable_point( vector<double> v1, vector<double> v2 ){
	int result = 0;
	vector<double> sc_v1 = scale_dvec(v1);
	vector<double> sc_v2 = scale_dvec(v2);
	vector<double> diff_tot( sc_v1.size() );
	for( int i=0;i<sc_v1.size();i++){
		diff_tot[i] = std::abs(sc_v1[i]) + std::abs(sc_v2[i]);
	}
	double hold_value = diff_tot[0];
	for( int i=0;i<sc_v1.size();i++){
		if ( diff_tot[i] < hold_value ){
			result = i;
			hold_value = diff_tot[i];
		}
	}
	return result;
}
/***********************************************************************/
std::ostream& operator<<(std::ostream& out, const traj_an& obj){
	
}
/***********************************************************************/
void traj_an::print(){
	
}
/***********************************************************************/
void UnitTest_traj_an(){}

/***********************************************************************/
