//interface.cpp

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
#include <string>
#include <vector>

#include "../include/global.h"
#include "../include/interface.h"
#include "../include/traj_analysis.h"
#include "../include/QCPinput.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::to_string;
using std::stoi;

/**********************************************************************/
interface::interface(){}
/**********************************************************************/
interface::~interface(){}
/**********************************************************************/
interface::interface(int argc, char* argv[]):
	m_argc(argc)							{
	
	for(int i =0;i<m_argc;i++){
		m_argv.emplace_back(argv[i]);
	}	
	cout << "Starting interface++ ! "	<< endl;
	
}
/**********************************************************************/
void interface::run(){
	
	/********************************/
	//Trajectory analysis with mdtraj
	if ( m_argv[1] == "-traj_geo"){
		traj_an Trajectory( m_argv[2] );
		Trajectory.mdtraj_geo();
	}
	/********************************/
	//Distance analysis from PDB traj
	else if ( m_argv[1] == "-dist_an" ){
		vector<int> ats;
		for(int i=3;i<m_argc;i++){
			ats.push_back( std::stoi(m_argv[i]) );
		}
		traj_an Trajectory(m_argv[2],ats);
		Trajectory.calc_distances( m_argv[2].c_str() );
	}
	/********************************/
	//Extract frame from a trajectory
	else if ( m_argv[1] == "-extract_frame" ){
		int arg3 = stoi(m_argv[3]);
		traj_an Trajectory;
		Trajectory.extract_frame(m_argv[2].c_str(),arg3);
	}
	/********************************/
	//Extract frames from a trajectory
	else if ( m_argv[1] == "-extract_frames" ){
		int arg3 = stoi(m_argv[3]);
		int arg4 = stoi(m_argv[4]);
		traj_an Trajectory;
		Trajectory.extract_frames(m_argv[2].c_str(),arg3,arg4);
	}
	/********************************/
	//Check if is ok the output from QCP 
	else if ( m_argv[1] == "-check_QCP_outs" ){
	}
	/********************************/
	//creat input from QCP calculations
	else if ( m_argv[1] == "-QCP_inp" ){
		if ( m_argv[2] == "gamess"){
			this->input_gamess();
		}
	}
	/********************************/
	else if ( m_argv[1] == "-MDprep" ){
		
	}
}
/**********************************************************************/
void interface::input_gamess(){
	// Always required options
	// - geometry extension
	// - QM method
	// - Basis Set
	//Default options 
	package prog 	= package::GAMESS;
	string geo_ext	= ".xyz";
	string QM_		= "am1";
	string basis_	= "am1";
	string mode		= "Normal"; // or FD
	int bcharge		= 0;
	int bmulti		= 1;
	string rtype	= "Energy";
	int charge_diff	= 1;
	
	for(int i =0;i<m_argc;i++){
		if ( m_argv[i] == "-gExt" ){
			geo_ext = m_argv[i+1];
		}
		else if ( m_argv[i] == "-QMm" ){
			QM_ = m_argv[i+1];
		}
		else if ( m_argv[i] == "-basis" ){
			basis_ = m_argv[i+1];
		}
		else if ( m_argv[i] == "-FDcalc" ){
			mode = "FD";
		}
		else if ( m_argv[i] == "-bchg" ){
			bcharge = stoi(m_argv[i+1]);
		}
		else if ( m_argv[i] == "-bmulti" ){
			bmulti = stoi(m_argv[i+1] );
		}
		else if ( m_argv[i] == "-runOP" ){
			rtype = m_argv[i+1];
		}
		else if ( m_argv[i] == "-chgDiff"){
			charge_diff = stoi(m_argv[i+1]);
		}
	}
	//this->print_options();
	QCPinput run_gms(geo_ext,rtype,basis_,QM_);
	if ( mode == "FD" ){
		run_gms.make_input_from_folder_FD(prog,bmulti,bcharge,charge_diff);
	}else{
		run_gms.make_input_from_folder(prog,bmulti,bcharge);
	}
	
}
/**********************************************************************/
void interface::print_options(){
	for(unsigned int i=0;i<m_argc;i++){
		cout << "#" << to_string(i+1) << " " << m_argv[i] << endl;
	}
}
/**********************************************************************/
void interface::help(){
	
}
/**********************************************************************/
void interface::test(){
	QCPinput inpsFromFolder(".xyz","optimize","TZV","MP2");
	inpsFromFolder.make_input_from_folder_FD(package::GAMESS,1,0,1);
}
////////////////////////////////////////////////////////////////////////