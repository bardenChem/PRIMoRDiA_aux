//read_traj.cpp

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
#include <iostream>
#include <fstream>

#include "../include/global.h"
#include "../include/pdbAtom.h"
#include "../include/residue.h"
#include "../include/pdbModel.h"
#include "../include/PDB.h"
#include "../include/read_traj.h"

#include <chemfiles.hpp>

using std::vector;
using std::string;
using std::move;
using std::cout;
using std::endl;


/*********************************************************/
ReadTraj::ReadTraj()		:
	natoms(0)				,
	nframes(0)				,
	traj_file("nonamed")	,
	Type( traj_type::NONE ) {
}
/*********************************************************/
ReadTraj::ReadTraj( const char* file_name )	:
	natoms(0)								,
	nframes(0)								,
	traj_file(file_name)					{
		
	if ( check_file_ext(".dcd",file_name) ){
		cout << "No topology file provided! " << endl;
	}else if( check_file_ext(".xtc",file_name) ){
		cout << "No topology file provided! " << endl;
	}else if( check_file_ext(".pdb",file_name) ){
		Type = traj_type::pdb;
		Positions = PDB(file_name);
		nframes = Positions.nModels;
	}
}
/*********************************************************/
ReadTraj::ReadTraj( const char* file_name, const char* topol_file ){
	traj_file = file_name;
	if( check_file_ext( ".pdb", topol_file) ){
		Type = traj_type::pdb;
		Positions = PDB(topol_file);
		nframes = Positions.nModels;
		Positions.MULTI = true;
	}	
	if ( check_file_ext(".dcd",file_name) ){
		Type = traj_type::dcd;
	}else if ( check_file_ext(".xtc",file_name) ){
		Type = traj_type::xtc;
	}
}
/*********************************************************/
ReadTraj::ReadTraj( const ReadTraj& rhs ):
	traj_file( rhs.traj_file )			,
	natoms( rhs.natoms )				,
	nframes( rhs.nframes )				,
	Positions( rhs.Positions )			,
	Type( rhs.Type )					{
}
/*********************************************************/
ReadTraj& ReadTraj::operator=( const ReadTraj& rhs ){
	if ( this != &rhs ){
		traj_file	= rhs.traj_file;
		natoms		= rhs.natoms;
		nframes		= rhs.nframes;
		Positions	= rhs.Positions;
		Type		= rhs.Type;
	}
	return *this;
}
/*********************************************************/
ReadTraj::ReadTraj( ReadTraj&& rhs ) noexcept:
	traj_file( move(rhs.traj_file) )		,
	natoms( move(rhs.natoms) )				,
	nframes( move(rhs.nframes ) )			,
	Positions( move(rhs.Positions ) )		,
	Type( move(rhs.Type) )					{
}
/*********************************************************/
ReadTraj& ReadTraj::operator=( ReadTraj&& rhs ) noexcept{
	if ( this != &rhs ){
		traj_file	= move(rhs.traj_file);
		natoms		= move(rhs.natoms);
		nframes		= move(rhs.nframes);
		Positions	= move(rhs.Positions);
		Type		= move(rhs.Type);
	}
	return *this;
}
/*********************************************************/
void ReadTraj::parse(){
	chemfiles::Trajectory traj( traj_file.c_str(), 'r' );
	nframes = traj.nsteps();
	Positions.init_models(nframes);
	chemfiles::Frame frame = traj.read();
	auto positions = frame.positions();
	unsigned poss = 0;
	for( unsigned i=0; i<nframes; i++){
		frame = traj.read_step(i);
		positions = frame.positions();
		for( unsigned m =0; m<Positions.models[i].monomers.size(); m++){
			for( unsigned n =0; n<Positions.models[i].monomers[m].r_atoms.size(); n++){
				Positions.models[i].monomers[m].r_atoms[n].xc = positions[poss][0];
				Positions.models[i].monomers[m].r_atoms[n].yc = positions[poss][1];
				Positions.models[i].monomers[m].r_atoms[n].zc = positions[poss++][2];
			}
		}
		poss = 0;
	}
}
/*********************************************************/
PDB ReadTraj::sample(unsigned interval){
	PDB sampled;
	sampled.basename = Positions.basename;
	sampled.MULTI = true;
	for(unsigned int i=0; i<nframes; i++ ){
		if ( i%interval == 0 ){
			sampled.add_model( Positions.models[i] );
			sampled.nModels++;
		}
	}
	return sampled;
}
/*********************************************************/
std::ostream& operator<<(std::ostream& out, const ReadTraj& obj){
	out << "Outputting information about 'ReadTraj' instanced object!\n"
		<< "File trajectory name: "	<< obj.traj_file
		<< "\nFile topology name: "	<< obj.Positions.basename
		<< "\nNumber of frames: "	<< obj.nframes;
	
	return out;
}
/*********************************************************/
void ReadTraj::print(){
	std::cout << *this << std::endl;
}
/*********************************************************/
void UnitTest_ReadTraj(){
	const char* dcd_file = "/home/igorchem/primordia-code/PRIMoRDiA_aux/test_data/structure/scan1d.dcd";
	const char* top_file = "/home/igorchem/primordia-code/PRIMoRDiA_aux/test_data/structure/frame0.pdb";
	const char* pdb_traj = "/home/igorchem/primordia-code/PRIMoRDiA_aux/test_data/structure/1l2y.pdb";
	
	ut_log.input_line("=======================================");
	ut_log.input_line("Starting unit test for 'pdbModel class!");
	ut_log.input_line("Default constructor: ");
	ReadTraj _traj_a;
	ut_log.data << _traj_a << std::endl;
	
	ut_log.input_line("PDB file constructor: ");
	ReadTraj _traj_b(pdb_traj);
	ut_log.data << _traj_b << std::endl;
	
	ut_log.input_line("Binary files and topology constructor constructor: ");
	ReadTraj _traj_c(dcd_file,top_file);
	_traj_c.parse();
	PDB _samp = _traj_c.sample(3);
	_samp.write_pdb("sampled_test.pdb");
	system("pymol sampled_test.pdb");
	ut_log.data << _traj_c << std::endl;
	
}
/////////////////////////////////////////////////////////////