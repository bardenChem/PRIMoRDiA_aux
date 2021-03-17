//geometry.cpp

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
#include <iostream>

#include "../include/global.h"
#include "../include/atom.h"
#include "../include/molecule.h"
#include "../include/XYZ.h"
#include "../include/PDB.h"
#include "../include/geometry.h"


using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::move;


///////////////////////////////////////////////////////////////////////
MOL2::MOL2(){}
/*********************************************************************/
MOL2::~MOL2(){}
/*********************************************************************/
MOL2::MOL2(const MOL2& rhs){}
/*********************************************************************/
MOL2& MOL2::operator=(const MOL2& rhs){}
/*********************************************************************/
MOL2::MOL2(MOL2&& rhs) noexcept{}
/*********************************************************************/
MOL2& MOL2::operator=(MOL2&& rhs) noexcept{}
///////////////////////////////////////////////////////////////////////
geometry::geometry():
	type(INVALID)	,
	cUnit(Ang)		{	
}
/*********************************************************************/
geometry::geometry(	const char* file_name)	:
	type(INVALID)							,
	cUnit(Ang)								{	
	
	if ( check_file_ext(".xyz",file_name ) ){
		type = xyz_;
	}else if ( check_file_ext(".pdb",file_name ) ){
		type = pdb_;
	}
	switch( type ){
		case xyz_:
			xyz = XYZ(file_name);
			Molecule = xyz.get_molecule();
			break;
		case pdb_:
			pdb = PDB(file_name);
			Molecule = pdb.get_system_from_model(0);
			break;
	}
		
}
/*********************************************************************/
geometry::~geometry(){}
/*********************************************************************/
geometry::geometry(const geometry& rhs) :
	type(rhs.type)						,
	xyz(rhs.xyz)						,
	pdb(rhs.pdb)						,
	mol2(rhs.mol2)						,
	Molecule(rhs.Molecule)				,
	cUnit(rhs.cUnit)					{		
}
/*********************************************************************/
geometry& geometry::operator=(const geometry& rhs){
	if ( this != &rhs ){
		type	= rhs.type;
		xyz		= rhs.xyz;
		pdb		= rhs.pdb;
		mol2	= rhs.mol2;
		Molecule= rhs.Molecule;
		cUnit	= rhs.cUnit;
	}
	return *this;
}
/*********************************************************************/
geometry::geometry(geometry&& rhs) noexcept:
	type(rhs.type)						,
	xyz( move(rhs.xyz) )				,
	pdb( move(rhs.pdb) )				,
	mol2( move(rhs.mol2) )				,
	Molecule( move(rhs.Molecule) )		,
	cUnit( move(rhs.cUnit) )			{
	
}
/*********************************************************************/
geometry& geometry::operator=(geometry&& rhs) noexcept{
	if ( this != &rhs ){
		type	= rhs.type;
		xyz		= move(rhs.xyz);
		pdb		= move(rhs.pdb);
		mol2	= move(rhs.mol2);
		Molecule= move(rhs.Molecule);
		cUnit	= move(rhs.cUnit);
	}
	return *this;	
}
/*********************************************************************/
void geometry::read_QCPinput(const char* file_name, std::string program){
	
}
/*********************************************************************/
void geometry::read_QCPoutput(const char* file_name, std::string program,bool last){
	
}
/*********************************************************************/
void geometry::convert_to_ang(){
	const double bohrtoang = 0.52917726;
	if ( cUnit == Bohr ){
		for(unsigned int i=0;i<Molecule.nAtoms;i++){
			Molecule.atoms[i].xc *= bohrtoang;
			Molecule.atoms[i].yc *= bohrtoang;
			Molecule.atoms[i].zc *= bohrtoang;
		}
	}

}
/*********************************************************************/
void geometry::convert_to_bohr(){
	const double angtobohr = 1.0/0.52917726;
	if ( cUnit == Ang ){
		for(unsigned int i=0;i<Molecule.nAtoms;i++){
			Molecule.atoms[i].xc *= angtobohr;
			Molecule.atoms[i].yc *= angtobohr;
			Molecule.atoms[i].zc *= angtobohr;
		}
	}
	
}
/*********************************************************************/
void geometry::write_to_file(std::string out_name,std::string format){
	if ( format == "xyz" ){
		if ( type == xyz_){
			xyz.write_xyz(out_name);
		}else{
			xyz = XYZ(Molecule);
			xyz.write_xyz(out_name);
		}		
	}else if ( format == "pdb"){
		if ( type == pdb_){
			pdb.write_pdb(out_name);
		}else{
			pdb.init_from_system(Molecule);
			pdb.write_pdb(out_name);
		}
	}
}
/*********************************************************************/
void geometry::center_coord(){
	double center_x = Molecule.atoms[0].xc;
	double center_y = Molecule.atoms[0].yc;
	double center_z = Molecule.atoms[0].zc;
	for(int i=0;i<Molecule.nAtoms;i++){
		Molecule.atoms[i].xc -= center_x;
		Molecule.atoms[i].yc -= center_y;
		Molecule.atoms[i].zc -= center_z;
	}
}
///////////////////////////////////////////////////////////////////////