//residue.cpp

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

#include <string>
#include <vector>
#include <map>

#include "../include/global.h"
#include "../include/pdbAtom.h"
#include "../include/residue.h"

using std::string;
using std::vector;

////////////////////////////////////////////////////////////////////////////
std::map<string, int> AAnHydrogens = {
	{"ALA", 5}, {"ARG",13}, {"ASN", 6}, {"ASP", 4}, 
	{"CYS", 5}, {"GLN", 8},	{"GLU", 6}, {"GLY", 3}, 
	{"HIS", 8}, {"ILE",11}, {"LEU",11},	{"LYS",13},
	{"MET", 9}, {"PHE", 9}, {"PRO", 7}, {"SER", 5}, 
	{"THR", 7}, {"TRP",10}, {"TYR", 9}, {"VAL", 9}
}
/*********************************************************/
std::map<string, AAres> R3name = {
	{"ALA", ALA}, {"ARG",ARG}, {"ASN", ASN}, {"ASP", ASP}, 
	{"CYS", CYS}, {"GLN",GLN}, {"GLU", GLU}, {"GLY", GLY}, 
	{"HIS", HIS}, {"ILE",ILE}, {"LEU", LEU}, {"LYS", LYS},
	{"MET", MET}, {"PHE",PHE}, {"PRO", PRO}, {"SER", SER}, 
	{"THR", THR}, {"TRP",TRP}, {"TYR", TYR}, {"VAL", VAL},
	{"HID", HIS}, {"HIE",HIS}, {"HIP", HIS}
}
/*********************************************************/
std::map<char, AAres>  R1name = {
	{'A', ALA}, {'R',ARG}, {'N', ASN}, {'D', ASP}, 
	{'C', CYS}, {'Q',GLN}, {'E', GLU}, {'G', GLY}, 
	{'H', HIS}, {'I',ILE}, {'L', LEU}, {'K', LYS},
	{'M', MET}, {'F',PHE}, {'P', PRO}, {'S', SER}, 
	{'T', THR}, {'W',TRP}, {'Y', TYR}, {'V', VAL}
}

///////////////////////////////////////////////////////////
residue::residue()	:
	type(UNK)		,
	AAname(OTH)     ,
	ligand(false)	,
	terminal(false)	,
	first(false)    ,
	nHydrogens(0)	,
	fCharge(0)		,
	pdb_index(0)	,
	nAtoms(0)		{
}
/*********************************************************/
residue::residue( vector<pdbAtom> resAtoms ):				
	r_atoms(resAtoms)						,
	type(resType)							,
	ligand(false)							,
	terminal(false)							,
	first(false)        					,
	pdb_index(0)							,
	fCharge(0.0)							{
		

	for(int i=0;i<resAtoms.size();i++){
		if ( resAtoms[i].is_hydrogen() ) nHydrogens++;
	}
	 
}
/*********************************************************/
residue::~residue(){}
/*********************************************************/
residue::residue(const residue& rhs):
	type(rhs.type)					,
	AAname(rhs.AAname)				,
	ligand(rhs.ligand)				,
	terminal(rhs.terminal)			,
	nHydrogens(rhs.nHydrogens)		,
	fCharge(rhs.fCharge)			,
	pdb_index(rhs.pdb_index)        ,
	nAtoms(rhs.nAtoms)				,
	r_atoms(rhs.r_atoms)			{	
}
/*********************************************************/
residue& residue::operator=(const residue& rhs){
	if ( this != &rhs ){
		res1n 		= rhs.res1n;				
		res3n 		= rhs.res3n;			
		type 		= rhs.type;
		AAname		= rhs.AAname;
		ligand		= rhs.ligand;				
		terminal	= rhs.terminal;	
		first		= rhs.first;
		pdb_index	= rhs.pdb_index;
		nHydrogens	= rhs.nHydrogens;		
		fCharge 	= rhs.fCharge;			
		nAtoms 		= rhs.nAtoms;				
		r_atoms 	= rhs.r_atoms;
	}
	return *this;
}
/*********************************************************/
residue::residue( residue&& rhs) noexcept:
	type( rhs.type )					,
	AAname(rhs.AAname)					,
	ligand(rhs.ligand)					,
	terminal(rhs.terminal)				,
	first(rhs.first)					,
	pdb_index(rhs.pdb_index )			,
	nHydrogens(rhs.nHydrogens)			,
	fCharge(rhs.fCharge)				,
	nAtoms( rhs.nAtoms) 				,
	r_atoms( move(rhs.r_atoms) )		{
	
}
/*********************************************************/
residue& residue::operator=( residue&& rhs) noexcept{
	if ( this != &rhs ){
		type 		= rhs.type;
		AAname		= rhs.AAname;
		ligand		= rhs.ligand;				
		terminal	= rhs.terminal;	
		first		= rhs.first;
		pdb_index	= rhs.pdb_index;
		nHydrogens	= rhs.nHydrogens;		
		fCharge 	= rhs.fCharge;			
		nAtoms 		= rhs.nAtoms;				
		r_atoms 	= move(rhs.r_atoms);
	}
	return *this;	
}
/*********************************************************/
bool residue::is_ion(){
	if ( r_atoms[0].res_name == "Cl-" ||
		 r_atoms[0].res_name == "Na+" ||
		 r_atoms[0].res_name == "K+"  ||
		 r_atoms[0].res_name == "Mg+" ||
		 r_atoms[0].res_name == "Zn+" ||
		 r_atoms[0].res_name == "SO4" ) {
	
		type = ION;
		return true;
	}
}
/*********************************************************/
void residue::set_charge(){
	
	int base_HN = get_AAnHy(AAname);
	
	if ( type == AA ){
		switch ( AAname ){
			case ASP:
			case GLU:
				if ( first || terminal ){
					fCharge = nHydrogens - base_HN - 2;				
				}else{
					fCharge = nHydrogens - base_HN  -1;
				} 
			break;
			//-------
			case ARG:
			case LYS:
				if ( first || terminal ) {
					fCharge = nHydrogens - base_HN;					
				}
				else {
					fCharge = nHydrogens - base_HN +1;
				}
			break;
			//-------
			case HIS:
				if ( first || terminal ) {
					fCharge = nHydrogens - base_HN;
				}
				else {
					fCharge = nHydrogens - base_HN +1;
				}
			break;
			//-------
			case CYS:
				fCharge = nHydrogens - base_HN - 1;
			break;
			//-------
			default:
				if ( first || terminal ){
					fCharge = nHydrogens - base_HN -1;
				}
				else {
					fCharge = 0;
				}
			break;			
		}		
	}
	else if ( type == ION ){
		if ( r_atoms[0].res_name == "Cl-" )
			fCharge = -1;
		else if  ( 	r_atoms[0].res_name == "Na+" ||
					r_atoms[0].res_name == "K+"  )
			fCharge = +1;
		else if (   r_atoms[0].res_name == "Mg+" ||
					r_atoms[0].res_name == "Zn+" )
			fCharge = +2;
		else if (  r_atoms[0].res_name == "SO4" ) 
			fCharge = -2;
		else fCharge = 0;
		
	}	
}
////////////////////////////////////////////////////////////////////