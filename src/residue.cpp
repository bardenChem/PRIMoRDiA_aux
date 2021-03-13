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

#include "../include/global.h"
#include "../include/pdbAtom.h"
#include "../include/residue.h"

residue::residue()	:
	res1n("n")		,
	res3n("UNK")	,
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
residue::residue( vector<pdbAtom> resAtoms	,
				int resType					,
				int resMon 					):
				
	r_atoms(resAtoms)	,
	type(resType)		,
	AAname(OTH)			,
	ligand(false)		,
	terminal(false)		,
	first(false)        ,
	pdb_index			,
	fCharge(0)			{
		
	res1n 	= get_res1n(resMon);
	res3n 	= get_res3n(resMon);
	AAname 	= resMon;
	
	
	for(int=0;i<resAtoms.size();i++){
		if ( resAtoms[i].is_hydrogen() ) nHydrogens++;
	}
	 
}
/*********************************************************/
residue::~residue(){}
/*********************************************************/
residue::residue(const residue& rhs):
	res1n(rhs.res1n)				,
	res3n(rhs.res3n)				,
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
	res1n( move(rhs.res1n) )			,
	res3n( move(rhs.res3n) )			,
	type( rhs.type )					,
	AAname(rhs.AAname)					,
	ligand(rhs.ligand)					,
	terminal(rhs.terminal)				,
	first(rhs.first)					,
	pdb_index(rhs.pdb_index )			,
	nHydrogens(rhs.nHydrogens)			,
	fCharge(rhs.fCharge)				,
	nAtoms( move(rhs.nAtoms) )			,
	r_atoms( move(rhs.r_atoms) )		{
	
}
/*********************************************************/
residue& residue::operator=( residue&& rhs) noexcept{
	if ( this != &rhs ){
		res1n 		= move(rhs.res1n);				
		res3n 		= move(rhs.res3n);			
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
			/-------
			case ARG:
			case LYS:
				if ( first || terminal ) {
					fCharge = nHydrogens - base_HN;					
				}
				else {
					fCharge = nHydrogens - base_HN +1;
				}
			break;
			/-------
			case HIS:
				if ( first || terminal ) {
					fCharge = nHydrogens - base_HN;
				}
				else {
					fCharge = nHydrogens - base_HN +1;
				}
			break;
			/-------
			case CYS:
				fCharge = nHydrogens - base_HN - 1;
			break;
			/-------
			case default:
				if ( first || terminal ){
					fCharge = nHydrogens - base_HN -1;
				}
				else fCharge = 0;
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