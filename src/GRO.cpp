//GRO.cpp

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
#include <iomanip>
#include <iostream>

#include "../include/global.h"
#include "../include/atom.h"
#include "../include/molecule.h"
#include "../include/pdbAtom.h"
#include "../include/residue.h"
#include "../include/pdbModel.h"
#include "../include/GRO.h"
#include "../include/Line.h"
using std::vector;
using std::string;
using std::move;
using std::endl; 

////////////////////////////////////////////////////////////////////////
GRO::GRO():
	nAtoms(0),
	sides(3),
	basename("noname"){	    
}
/**********************************************************************/
GRO::~GRO(){}
/**********************************************************************/
GRO::GRO(const char* gro_file):
	nAtoms(0)			      ,
	basename(gro_file)        {
		
	sides.resize(3);
	
	if ( !check_file_ext(".gro",gro_file) )	{
		std::cout << "Warning! The file has wrong extension name!" << std::endl;
	}		
	
	vector<pdbAtom> tmp_atoms;
	
	char gro_line[100];

	if ( IF_file(gro_file) ){
		
		std::ifstream buf(gro_file);
		
		buf.getline(gro_line,100);
		buf.getline(gro_line,100);
		string _grostring(gro_line);
		nAtoms = std::stoi(_grostring);
		
		
		for(unsigned i=0; i<nAtoms; ++i){
			buf.getline(gro_line,100);
			string _grostring(gro_line);	
			unsigned _resi =  std::stoi( _grostring.substr(0,5) );
			string 	 _resn =  _grostring.substr(5,3);
			string   _ele = _grostring.substr(11,4);
			double   _tmpXC = std::stod(_grostring.substr(20,8) );
			double   _tmpYC = std::stod(_grostring.substr(28,8) );
			double   _tmpZC = std::stod(_grostring.substr(36,8) );
			pdbAtom  _atom(_resn,_resi, _ele, _tmpXC,_tmpYC,_tmpZC);
			_atom.indx = i+1;
			atoms.emplace_back( std::move(_atom) ); 
			
		}
		buf.getline(gro_line,100);
		std::istringstream iss(gro_line);
		iss >> sides[0];
		iss >> sides[1];
		iss >> sides[2];		
	}else{
		std::cout << "GRO file not open!" << std::endl;
	}
}
/**********************************************************************/
GRO::GRO(const GRO& rhs):
	nAtoms(rhs.nAtoms),
	sides(rhs.sides)  , 
	atoms(rhs.atoms)  ,
	basename(rhs.basename){	
}
/**********************************************************************/
GRO& GRO::operator=(const GRO& rhs){
	if (this != &rhs){
		nAtoms 		= rhs.nAtoms;
		sides  		= rhs.sides;
		atoms  		= rhs.atoms;
		basename 	= rhs.basename;		
	}
	return *this;
	
}
/**********************************************************************/
GRO::GRO(GRO&& rhs) noexcept:
	nAtoms(rhs.nAtoms),
	sides( std::move(rhs.sides) ), 
	atoms( std::move(rhs.atoms) ),
	basename( std::move(rhs.basename))
	{	
}
/**********************************************************************/
GRO& GRO::operator=(GRO&& rhs) noexcept{
	if (this != &rhs){
		nAtoms 		= rhs.nAtoms;
		sides  		= std::move(rhs.sides);
		atoms  		= std::move(rhs.atoms);
		basename 	= std::move(rhs.basename);
	}
	return *this;	
}
/**********************************************************************/
void GRO::write_gro(std::string out_name){
	std::ofstream gro_file;
	gro_file.open( out_name.c_str() );
	gro_file << std::fixed;
	gro_file << "GRO file written by LQQCMMtools created by barden.igor@gmail.com" << endl;
	
	gro_file << atoms.size() << endl;
	
	for(unsigned i=0;i<atoms.size();i++){
		
		gro_file << std::setw(5) << atoms[i].res_indx;
        gro_file << std::left << std::setw(5) << atoms[i].res_name << std::right;
        gro_file << std::setw(5) << atoms[i].atom_name;
        gro_file << std::setw(5) << i;
        gro_file << std::fixed << std::setprecision(3);
        gro_file << std::setw(8) << atoms[i].xc;
        gro_file << std::setw(8) << atoms[i].yc;
        gro_file << std::setw(8) << atoms[i].zc;
		gro_file << endl;
	}
	gro_file << sides[0] <<  "  " << sides[1] << "  " << sides[2] << endl;
	gro_file.close();

}
/**********************************************************************/
void GRO::init_from_system(const molecule& mol){
	
}
/**********************************************************************/
molecule GRO::get_system_from_gro(unsigned int model){
	
}
/**********************************************************************/
std::vector<molecule> GRO::get_systems(){
	
}
/**********************************************************************/
pdbModel GRO::get_pdb_from_gro(){	
	for(unsigned i=0;i<nAtoms;i++){
			atoms[i].xc *=10;
			atoms[i].yc *=10;
			atoms[i].zc *=10;		
	}
	pdbModel _pdb(this->atoms);	
	return _pdb;
}
/**********************************************************************/
std::ostream& operator<<(std::ostream& out, const GRO& obj){
	
}
/**********************************************************************/
void GRO::print(){

}
/******************************************************************************/
void UnitTest_GRO(){
	const char* gro_md1  = "/home/igorchem/CCDIR/PETROBRAS-F3/AC_DM/AC/md_1.gro";
	GRO _gro_test(gro_md1);
	_gro_test.write_gro("gro_test_write.gro");
	pdbModel _pdbFile = _gro_test.get_pdb_from_gro();
	_pdbFile.write_model("pdb_from_gro.pdb");
	
}