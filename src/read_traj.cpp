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
#include "../include/geometry.h"
#include "../include/read_traj.h"

#include <chemfiles.hpp>

using std::vector;
using std::string;
using std::move;
using std::cout;
using std::endl;

/*********************************************************/
positions3D::positions3D()	:
	natoms(0)				,
	nframes(0)				{
}
/*********************************************************/
positions3D::~positions3D(){}
/*********************************************************/
positions3D::positions3D( const positions3D& rhs )	:
	natoms(rhs.natoms)								,
	nframes(0)										,
	xc(rhs.xc)										,
	yc(rhs.yc)										,
	zc(rhs.zc)										{
}
/*********************************************************/
positions3D& positions3D::operator=( const positions3D& rhs){
	if ( this != &rhs ){
		natoms	= rhs.natoms;
		nframes	= rhs.nframes;
		xc		= rhs.xc;
		yc		= rhs.yc;
		zc		= rhs.zc;
	}
	return *this;
}
/*********************************************************/
positions3D::positions3D( positions3D&& rhs ) noexcept	:
	natoms(rhs.natoms)									,
	nframes(rhs.nframes)								,
	xc( move(rhs.xc) )									,
	yc( move(rhs.yc) )									,
	zc( move(rhs.zc) )									{
}
/*********************************************************/
positions3D& positions3D::operator=( positions3D&& rhs ) noexcept{
	if ( this != &rhs ){
		natoms	= rhs.natoms;
		nframes	= rhs.nframes;
		xc		= move(rhs.xc);
		yc		= move(rhs.yc);
		zc		= move(rhs.zc);
	}
	return *this;
}
/*********************************************************/
double positions3D::get_distance(unsigned atom1, unsigned atom2, unsigned frame){
	double dx = xc[frame][atom1] - xc[frame][atom2]; 
	double dy = yc[frame][atom1] - yc[frame][atom2]; 
	double dz = zc[frame][atom1] - zc[frame][atom2];
	return sqrt( dx*dx + dy*dy + dz*dz );
}
/*********************************************************/
std::vector<double> positions3D::get_distances(unsigned atom1, unsigned atom2){
	std::vector<double> distances;
	double temp = 0.0;
	for( unsigned i=0; i<nframes; i++){
		temp = this->get_distance(atom1, atom2, i);
		distances.push_back(temp);
	}
	return distances;
}
/*********************************************************/
double positions3D::avg_distance(unsigned atom1, unsigned atom2){
	std::vector<double> distances = this->get_distances(atom1,atom2);
	double avg = mean_dvec(distances);
	return avg;
}
/*********************************************************/
void positions3D::resize(unsigned nframes){
	xc.resize(nframes);
	yc.resize(nframes);
	zc.resize(nframes);
}
/*********************************************************/
std::ostream& operator<<(std::ostream& out, const positions3D& obj){
	
}
/*********************************************************/
void positions3D::print(){
	
}
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
	}else if( check_file_ext(".gro",file_name) ){
		Type = traj_type::gro;
		geometry top( file_name );
		topology = top;
	}else if( check_file_ext(".pdb",file_name) ){
		Type = traj_type::pdb;
		PDB top(file_name);
		PDB frame0;
		frame0.add_model( top.models[0] );
		topology.pdb = frame0;		
		coordinates.resize(topology.pdb.nModels);
		unsigned c = 0;
		for(unsigned i=0; i<topology.pdb.nModels; i++ ){
			for(unsigned j=0; j<topology.pdb.models[i].nResidues; j++ ){
				for(unsigned k=0; k<topology.pdb.models[i].monomers[j].nAtoms; k++ ){
					coordinates.xc[i][c]	= topology.pdb.models[i].monomers[j].r_atoms[k].xc;
					coordinates.yc[i][c]	= topology.pdb.models[i].monomers[j].r_atoms[k].yc;
					coordinates.zc[i][c++]	= topology.pdb.models[i].monomers[j].r_atoms[k].zc;
				}
			}
		}
	}
}
/*********************************************************/
ReadTraj::ReadTraj( const char* file_name, const char* topol_file ){
	traj_file = file_name;
	geometry top( topol_file );
	topology = top;
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
	coordinates( rhs.coordinates)		,
	topology( rhs.topology )			,
	Type( rhs.Type )					{
}
/*********************************************************/
ReadTraj& ReadTraj::operator=( const ReadTraj& rhs ){
	if ( this != &rhs ){
		traj_file	= rhs.traj_file;
		natoms		= rhs.natoms;
		nframes		= rhs.nframes;
		coordinates = rhs.coordinates;
		topology	= rhs.topology;
		Type		= rhs.Type;
	}
	return *this;
}
/*********************************************************/
ReadTraj::ReadTraj( ReadTraj&& rhs ) noexcept:
	traj_file( move(rhs.traj_file) )		,
	natoms( move(rhs.natoms) )				,
	nframes( move(rhs.nframes ) )			,
	coordinates( move(rhs.coordinates) )	,
	topology( move(rhs.topology ) )			,
	Type( move(rhs.Type) )					{
}
/*********************************************************/
ReadTraj& ReadTraj::operator=( ReadTraj&& rhs ) noexcept{
	if ( this != &rhs ){
		traj_file	= move(rhs.traj_file);
		natoms		= move(rhs.natoms);
		nframes		= move(rhs.nframes);
		coordinates = move(rhs.coordinates);
		topology	= move(rhs.topology);
		Type		= move(rhs.Type);
	}
	return *this;
}
/*********************************************************/
void ReadTraj::parse(){
	chemfiles::Trajectory file( traj_file.c_str() );
	chemfiles::Frame frame = file.read();
	auto positions = frame.positions();
}
/*********************************************************/
PDB ReadTraj::sample(unsigned interval){
	
}
/*********************************************************/
std::ostream& operator<<(std::ostream& out, const ReadTraj& obj){
	
}
/*********************************************************/
void ReadTraj::print(){
	
}
/*********************************************************/
void UnitTest_positions3D(){
	
}
void UnitTest_ReadTraj(){
	
}