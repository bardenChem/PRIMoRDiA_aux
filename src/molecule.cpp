//molecule.cpp

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
#include <cmath>
#include <fstream>
#include <iostream>

#include "../include/global.h"
#include "../include/atom.h"
#include "../include/molecule.h"


using std::vector;
using std::string;
using std::move;

/***********************************************************/
molecule::molecule():
	nAtoms(0)		,
	nElectrons(0)	,
	name("nonamed")	,
	type("none")	,
	fCharge(0.0)	,
	molar_mass(0.0)	{

	for(unsigned int i=0;i<3;i++){
		ver_inf[i] = 0.000;
		ver_sup[i] = 0.000;
	}
}
/**********************************************************/
molecule::molecule(std::vector<atom> ats,
				std::string nme		):
				
	nAtoms( ats.size() ),
	nElectrons(0)		,
	name(nme)			,
	type("none")		,
	fCharge(0.0)		,
	molar_mass(0.0)		,
	atoms(ats)			{
		
	for(unsigned int i=0;i<nAtoms;i++){
		nElectrons	+= atoms[i].aNmb;
		molar_mass	+= atoms[i].aMass;
		fCharge		+= atoms[i].pCharge;
		if( i == 0 ){
			ver_inf[0] = atoms[i].xc;
			ver_inf[0] = atoms[i].yc;
			ver_inf[0] = atoms[i].zc;
			ver_sup[0] = atoms[i].xc;
			ver_sup[0] = atoms[i].yc;
			ver_sup[0] = atoms[i].zc;
		}
		else{
			if ( ver_inf[0] > atoms[i].xc ) ver_inf[0] = atoms[i].xc;
			if ( ver_inf[1] > atoms[i].yc ) ver_inf[1] = atoms[i].yc;
			if ( ver_inf[2] > atoms[i].zc ) ver_inf[2] = atoms[i].zc;
			if ( ver_sup[0] < atoms[i].xc ) ver_sup[0] = atoms[i].xc;
			if ( ver_sup[1] < atoms[i].yc ) ver_sup[1] = atoms[i].yc;
			if ( ver_sup[2] < atoms[i].zc ) ver_sup[2] = atoms[i].zc;
		}		
	}	
}
/**********************************************************/
molecule::~molecule(){}
/**********************************************************/
molecule::molecule(const molecule& rhs):
	nAtoms(rhs.nAtoms)			,
	nElectrons(rhs.nElectrons)	,
	name(rhs.name)				,
	type(rhs.type)				,
	atoms(rhs.atoms)			,
	molar_mass(rhs.molar_mass)  {
		
	for(int i=0;i<3;i++){
		ver_inf[i] = rhs.ver_inf[i];
		ver_sup[i] = rhs.ver_sup[i];
	}	
}
/*********************************************************/
molecule::molecule(molecule&& rhs) noexcept:
	nAtoms(rhs.nAtoms)				,
	nElectrons(rhs.nElectrons)		,
	name( move(rhs.name) )			,
	type( move(rhs.type) )			,
	atoms( move(rhs.atoms) ) 		,
	molar_mass(rhs.molar_mass)		{
	
	for(int i=0;i<3;i++){
		ver_inf[i] = rhs.ver_inf[i];
		ver_sup[i] = rhs.ver_sup[i];
	}
}
/**********************************************************/
molecule& molecule::operator=(const molecule& rhs){
	if( this != &rhs ){
		nAtoms 		= rhs.nAtoms;
		nElectrons 	= rhs.nElectrons;
		name 		= rhs.name;		
		type 		= rhs.type;
		atoms 		= rhs.atoms;
		molar_mass 	= rhs.molar_mass;
		
		for(int i=0;i<3;i++){
			ver_inf[i] = rhs.ver_inf[i];
			ver_sup[i] = rhs.ver_sup[i];
		}
	}
	return *this;
}
/*********************************************************/
molecule& molecule::operator=(molecule&& rhs) noexcept{
	if( this != &rhs ){
		nAtoms 		= rhs.nAtoms;
		nElectrons 	= rhs.nElectrons;
		name 		= move(rhs.name);
		type 		= move(rhs.type);
		atoms 		= move(rhs.atoms);
		molar_mass 	= rhs.molar_mass;
		
		for(int i=0;i<3;i++){
			ver_inf[i] = rhs.ver_inf[i];
			ver_sup[i] = rhs.ver_sup[i];
		}
	}
	return *this;
}
/*********************************************************/
void molecule::add_atom(atom a){
	atoms.emplace_back( move(a) );
	nAtoms++;
}
/*********************************************************/
void molecule::add_atom(double x		,
					  double y		, 
					  double z		, 
					  std::string el){
						  
	atoms.emplace_back(x,y,z,el);
	nAtoms++;
}
/*********************************************************/
void molecule::remove_atom(unsigned int i){
	atoms.erase( atoms.begin()+i );
	nAtoms--;
}
/*********************************************************/
std::ostream& operator<<(std::ostream& out, const molecule& obj){
	out <<"Outputting information about molecule object"
		<<"\nmolecule name " 			<< obj.name
		<<"\n\tNumber of atoms:"		<< obj.nAtoms
		<<"\n\tNumber of electrons:"	<< obj.nElectrons
		<<"\n\tFormal charge: "			<< obj.fCharge
		<<"\n\tMolar Mass: "			<< obj.molar_mass
		<<"\n\tSize of the atom container" << obj.atoms.size();
		
	return out;
}
/*********************************************************/
void molecule::print(){
	std::cout << *this << std::endl;
}
/*********************************************************/
void UnitTest_molecule(){
	
}
////////////////////////////////////////////////////////////