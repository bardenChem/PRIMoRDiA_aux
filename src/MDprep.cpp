//MDprep.cpp

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

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <"../include/global.h"
#include <"../include/system.h"
#include <"../include/geometry.h"
#include <"../include/MDprep.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::move;

//=================================
MDprep::MDprep()		:
	mdPKG(undef)		,
	mdTask(undefined)	,
	net_chg(0)			,
	ref_temp(310.15)	,
	prod_time(10)		,
	top_file("")		,
	crd_file("")		,
	sys_basename("")	{	
}
/******************************************************/
MDprep::MDprep(std::string bs_name	,
				mdPKG pkg			,
				mdTask Job			):
	mdPKG(pkg)						,
	mdTask(Job)						,
	net_chg(0)						,
	ref_temp(310.15)				,
	prod_time(10)					,
	sys_basename(bs_name)			{	
	
	switch ( pkg ){
		case AMBER:
			top_file = change_extension( sys_basename.c_str(), ".prmtop" );
			crd_file = change_extension( sys_basename.c_str(), ".crd" );
		break;
		case GROMACS:
			top_file = change_extension( sys_basename.c_str(), ".top" );
			crd_file = change_extension( sys_basename.c_str(), ".gro" );
		break;
		case default:
			top_file = "_";
			crd_file = "_";
	}
}
/******************************************************/
MDprep::~MDprep(){}
/******************************************************/
void MDprep::prepare_ligand(int nligand		, 
							bool ambTools)	{	
}
/******************************************************/
void MDprep::prepare_complex(int nligand){
	
}
/******************************************************/
void MDprep::built_topology(){
	
}
/******************************************************/
void MDprep::prepare_minimization(){
	
}
/******************************************************/
void MDprep::prepare_equilibration(){
	
}
/******************************************************/
void MDprep::prepare_production(){
	
}
/******************************************************/
void MDprep::organize_dir_files(){
	
}
////////////////////////////////////////////////////////