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

using std::cout;
using std::endl;
using std::string;
using std::vector;

/**********************************************************************/
interface::interface(){}
/**********************************************************************/
interface::~interface(){}
/**********************************************************************/
interface::interface(int argc, char* argv[]):
	m_argc(argc)						{
	
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
		traj Trajectory( m_argv[2] );
		Trajectory.mdtraj_geo();
	}
	/********************************/
	//Distance analysis from PDB traj
	else if ( m_argv[1] == "-dist_an" ){
		vector<int> ats;
		for(int i=3;i<m_argc;i++){
			ats.push_back( std::stoi(m_argv[i]) );
		}
		traj Trajectory(m_argv[2],ats);
		Trajectory.calc_distances( m_argv[2].c_str() );
	}
	/********************************/
	//Extract frame from a trajectory
	else if ( m_argv[1] == "-extract_frame" ){
		int arg3 = stoi(m_argv[3]);
		traj Trajectory;
		Trajectory.extract_frame(m_argv[2].c_str(),arg3);
	}
	/********************************/
	//Extract frames from a trajectory
	else if ( m_argv[1] == "-extract_frames" ){
		int arg3 = stoi(m_argv[3]);
		int arg4 = stoi(m_argv[4]);
		traj Trajectory;
		Trajectory.extract_frames(m_argv[2].c_str(),arg3,arg4);
	}
	/********************************/
	//Check if is ok the output from QCP 
	else if ( m_argv[1] == "-check_QCP_outs" ){
		
	}
	/********************************/
	//creat input from QCP calculations
	else if ( m_argv[1] == "-built_QCP_inp" ){
		
	}
	/********************************/
	else if ( m_argv[1] == "-MDprep" ){
		
	}
}
/**********************************************************************/
void interface::print_options(){
	
}
/**********************************************************************/
void interface::help(){
	
}
////////////////////////////////////////////////////////////////////////