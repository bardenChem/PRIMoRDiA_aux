//ooccupy.cpp

/*********************************************************************/
/* This source code file is part of OOCCuPy++ software project created 
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
#include <ctime>
#include <string>
#include <vector>

#include "../include/global.h"
#include "../include/ooccupy.h"
#include "../include/traj.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

/**********************************************************************/
ooccupy::ooccupy(){
	
}
/**********************************************************************/
ooccupy::~ooccupy(){
	
}
/**********************************************************************/
ooccupy::ooccupy(int argc, char* argv[]):
	m_argc(argc)						{
	
	for(int i =0;i<m_argc;i++){
		m_argv.emplace_back(argv[i]);
	}
	
	time_t now	= time(0);
	char* dt	= ctime(&now);
	cout << "Starting OOCCuPy++ ! "	<< endl;
	
}
/**********************************************************************/
void ooccupy::run(){
	if ( m_argv[1] == "-traj_geo"){
		traj Trajectory( m_argv[2] );
		Trajectory.mdtraj_geo();
	}
	else if ( m_argv[1] == "-dist_an" ){
		vector<int> ats;
		for(int i=3;i<m_argc;i++){
			ats.push_back( std::stoi(m_argv[i]) );
		}
		traj Trajectory(m_argv[2],ats);
		Trajectory.calc_distances( m_argv[2].c_str() );
	}
	else if ( m_argv[1] == "-extract_frame"){
		int arg3 = stoi(m_argv[3]);
		traj Trajectory;
		Trajectory.extract_frame(m_argv[2].c_str(),arg3);
	}
	else if ( m_argv[1] == "-extract_frames"){
		int arg3 = stoi(m_argv[3]);
		int arg4 = stoi(m_argv[4]);
		traj Trajectory;
		Trajectory.extract_frames(m_argv[2].c_str(),arg3,arg4);
	}
}
/**********************************************************************/