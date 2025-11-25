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
#include <omp.h>


#include "../include/global.h"
#include "../include/pdbAtom.h"
#include "../include/residue.h"
#include "../include/pdbModel.h"
#include "../include/PDB.h"
#include "../include/read_traj.h"
#include "../include/GRO.h"

#include <chemfiles.hpp>

using std::vector;
using std::string;
using std::move;
using std::cout;
using std::endl;

/********************************************************/
TrjCRD::TrjCRD(){}
/********************************************************/
TrjCRD::TrjCRD( pdbModel _topology):
	topology(_topology){
}
/********************************************************/
TrjCRD::TrjCRD( const TrjCRD& rhs ):
	xc(rhs.xc)						,
	yc(rhs.yc)						,
	zc(rhs.zc)						,
	topology(rhs.topology)			{	
}
/********************************************************/
TrjCRD& TrjCRD::operator=( const TrjCRD& rhs ){
	if ( this != &rhs ){
		xc = rhs.xc;
		yc = rhs.yc;
		zc = rhs.zc;
		topology = rhs.topology;
	}
	return *this;
}
/********************************************************/
TrjCRD::TrjCRD(  TrjCRD&& rhs ) noexcept:
	xc( move (rhs.xc) ),
	yc( move (rhs.yc) ),
	zc( move (rhs.zc) ),
	topology( move (rhs.topology) ){		
}
/********************************************************/
TrjCRD& TrjCRD::operator=( TrjCRD&& rhs ) noexcept{
	if ( this != &rhs ){
		xc = move(rhs.xc);
		yc = move(rhs.yc);
		zc = move(rhs.zc);
		topology = move(rhs.topology);
	}
	return *this;
}
/********************************************************/
void TrjCRD::init_frames(unsigned _nFrames){
	xc.resize(_nFrames);
	yc.resize(_nFrames);
	zc.resize(_nFrames);
	
	if (topology.nAtoms > 0){
		for(unsigned i=0;i<_nFrames;i++){
			xc[i].resize(topology.nAtoms);
			yc[i].resize(topology.nAtoms);
			zc[i].resize(topology.nAtoms);
		}
	}
}
/********************************************************/
TrjCRD::~TrjCRD(){}
/********************************************************/
pdbModel TrjCRD::create_pdb( unsigned _frame){
	pdbModel new_pdb(topology);
	unsigned pos = 0;
	for(unsigned i=0; i<new_pdb.monomers.size();i++){
		for(unsigned j=0; j<new_pdb.monomers[i].r_atoms.size(); j++){
			new_pdb.monomers[i].r_atoms[j].xc = xc[_frame][pos];
			new_pdb.monomers[i].r_atoms[j].yc = yc[_frame][pos];
			new_pdb.monomers[i].r_atoms[j].zc = zc[_frame][pos++];
		}	
	}
	return new_pdb;	
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
	}else if( check_file_ext(".pdb",file_name) ){
		Type = traj_type::pdb;
		PDB _topology(file_name);
		Positions = TrjCRD(_topology.models[0]);
	}
}
/*********************************************************/
ReadTraj::ReadTraj( const char* file_name, const char* topol_file ){
	traj_file = file_name;
	if( check_file_ext( ".pdb", topol_file) ){
		PDB _topology(topol_file);
		Positions = TrjCRD(_topology.models[0]);
		nframes   = 0;
	}else if ( check_file_ext( ".gro", topol_file) ){
		GRO _topology(topol_file);
		Positions = TrjCRD(_topology.get_pdb_from_gro()); 
	}
	natoms = Positions.topology.nAtoms;
	if ( check_file_ext(".dcd",file_name) ){
		Type = traj_type::dcd;
	}else if ( check_file_ext(".xtc",file_name) ){
		Type = traj_type::xtc;
	}else if ( check_file_ext(".pdb",file_name) ){
		Type = traj_type::pdb;
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
	Positions.init_frames(nframes);
	chemfiles::Frame frame = traj.read();
	auto positions = frame.positions();
	unsigned poss = 0;
	for( unsigned i=0; i<nframes; i++){
		frame = traj.read_step(i);
		positions = frame.positions();
		for( unsigned m =0; m<Positions.topology.monomers.size(); m++){
			for( unsigned n =0; n<Positions.topology.monomers[m].r_atoms.size(); n++){
				Positions.xc[i][poss] = positions[poss][0];
				Positions.yc[i][poss] = positions[poss][1];
				Positions.zc[i][poss] = positions[poss][2];
				poss++;
			}
		}
		poss = 0;
	}
}
/*********************************************************/
PDB ReadTraj::sample(unsigned interval){
	PDB sampled;
	sampled.MULTI = true;
	for(unsigned int i=0; i<nframes; i++ ){
		if ( i%interval == 0 ){
			sampled.add_model( Positions.create_pdb(i) );
			sampled.nModels++;
		}
	}
	return sampled;
}
/*********************************************************/
PDB ReadTraj::sample_chunk(unsigned _init, unsigned _final){
	PDB sampled;
	sampled.MULTI = true;
	
	if (_final > nframes){
		cout << "Final number of the chunk greater than the number of frames!!" << endl;
		return sampled;
	}		
	for(unsigned int i=_init; i<_final; i++ ){
		sampled.add_model( Positions.create_pdb(i) );
		sampled.nModels++;
	}
	return sampled;
}
/*********************************************************/
void ReadTraj::analysis_ac_from_molecules(unsigned _res_indx, std::string _molecule_name, double _radius, double _prune){
	
	std::vector<unsigned> _list_molecules = Positions.topology.get_res_list(_molecule_name);
	std::vector< std::vector<double> > distances;
	distances.resize(Positions.xc.size());
	
	unsigned i = 0;
	#pragma omp parallel
	{
	#pragma omp for private(i)
		for(i=0;i<Positions.xc.size();i++){
		
			pdbModel model_frame = Positions.create_pdb(i);
			for (unsigned j=0; j<_list_molecules.size();j++){			
				double dist = model_frame.monomers[_res_indx-1].smallest_distance( model_frame.monomers[ _list_molecules[j] ] ); 
				distances[i].push_back(dist);
			}
			std::sort(distances[i].begin(), distances[i].end());
		}
	}
	
	std::ofstream outfile("sorted_distances.dat");
    if(!outfile.is_open()) {
        std::cerr << "Error: Could not open file for writing!" << std::endl;
        return;
    }
    
    // Write header (optional)
    outfile << "#Frame ";
	
	for(unsigned i=0; i< _list_molecules.size(); i++){
		outfile << Positions.topology.monomers[_list_molecules[i]].name << _list_molecules[i] << " ";
	}
	outfile << "\n";    
    // Write data - one frame per line
    for(unsigned i=0; i < distances.size(); i++) {
        outfile << i << " ";  // Frame number
        for(double dist : distances[i]) {
            outfile << dist << " ";
        }
        outfile << "\n";
    }
    outfile.close();
	
	PDB sampled;
	sampled.MULTI = true;
	
	for(unsigned i =0;i<distances.size();i++){
		if ( distances[i][0]  < _radius ) {
			sampled.add_model( Positions.create_pdb(i) );
		}
	}
	string res_name = Positions.topology.monomers[_res_indx-1].name;
	if ( res_name[2] ==' '){
		res_name = res_name.substr(0,2);
	}
	unsigned center = Positions.topology.monomers[_res_indx-1].r_atoms[0].indx -1;
	sampled.basename = "analysis"+res_name+"_"+_molecule_name;
	if (_prune > 0.0 ) {
		string res_str = std::to_string(center);
		string radius  = std::to_string(_prune);
		std::vector<string> _parameters = { res_str, radius, "false","false"};
		sampled.iterate_models("prune_atoms",_parameters);
	}
	
	sampled.split_models_in_files();
		
	//std::ofstream R_script;
	
}
/*********************************************************/
std::ostream& operator<<(std::ostream& out, const ReadTraj& obj){
	out << "Outputting information about 'ReadTraj' instanced object!\n"
		<< "File trajectory name: "	<< obj.traj_file
		<< "\nNumber of frames: "	<< obj.nframes;
	
	return out;
}
/*********************************************************/
void ReadTraj::print(){
	std::cout << *this << std::endl;
}
/*********************************************************/
void UnitTest_ReadTraj(){
	
	const char* xtc_file = "/home/igorchem/CCDIR/PETROBRAS-F3/AC_DM/AC/traj.xtc";
	const char* gro_file = "/home/igorchem/CCDIR/PETROBRAS-F3/AC_DM/AC/traj.gro";
	
	ReadTraj _xtc_test(xtc_file,gro_file);
	_xtc_test.parse();
	_xtc_test.analysis_ac_from_molecules(466,"CO2",4.2,30.0);
	
	/*
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
	*/
	
}
/////////////////////////////////////////////////////////////