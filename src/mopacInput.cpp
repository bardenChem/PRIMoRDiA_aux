//gamessInput.cpp

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
#include <experimental/filesystem>

#include "../include/global.h"
#include "../include/atom.h"
#include "../include/molecule.h" 
#include "../include/geometry.h"
#include "../include/mopacInput.h"
#include "../include/pdbModel.h"
#include "../include/residue.h"
#include "../include/pdbAtom.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::move;

namespace fs = std::experimental::filesystem;

/*************************************************************/
mopac_input::mopac_input()	:
	charge(0)				,
	COSMO(false)			,
	MOZYME(false)			,
	QMMM(false)				,
	multiplicity(SINGLET)	,
	method(AM1)				,
	rtype(mSCF)				{
}
/*************************************************************/
mopac_input::~mopac_input(){}
/*************************************************************/
void mopac_input::init( int chg				,
						unsigned int mpcty	,
						string solvent		,
						string lmo			,
						string Method)		{
							
	if ( solvent== "h2o"	) COSMO  = true;
	if ( lmo 	== "mozyme"	) MOZYME = true;
	
	
	keywords.emplace_back(Method);
	keywords.emplace_back(" 1SCF ALLVECS VECTOR AUX LARGE LET");
		
	if ( MOZYME ){
		keywords.emplace_back(" MOZYME");
		//keywords.emplace_back(" charge=n");
	}
	else{
		keywords.emplace_back(" charge=");
		keywords.emplace_back( std::to_string(charge) );
		switch ( mpcty ){
		case 2:
			multiplicity = DOUBLET;
			keywords.emplace_back(" Doublet");
		break;
		case 3:
			multiplicity = TRIPLET;
			keywords.emplace_back(" Triplet");
		break;
		case 4:
			multiplicity = QUARTET;
			keywords.emplace_back(" Quartet");
		break;
		case 5:
			multiplicity = QUINTET;
			keywords.emplace_back(" Quintet");	
		break;
		default:
			keywords.emplace_back(" Singlet");
		break;
	}
	}
	if ( COSMO  ) keywords.emplace_back(" eps=78.4");
	
	charge = chg;	
	
	for (int i=0;i<keywords.size();i++){
		out_file << keywords[i];
	}
	out_file << endl;
	out_file << endl;
}
/*************************************************************/
void mopac_input::molin_init(pdbModel& qc_region, pdbModel& mm_region, std::string Method){
	
}
/*************************************************************/
void mopac_input::mark_charge(molecule&mol, std::string out_name, pdbModel& topol, std::string _residue, int charge){
	
	std::vector<int> marked_atoms; 	
	std::string mark;
	
	if 		( charge > 0 ) mark = "(+)";
	else if ( charge < 0 ) mark = "(-)";
	unsigned int cnt = 0;
	
	for (unsigned int i=0;i<topol.nResidues;i++){
		if ( topol.monomers[i].name == _residue ){
			marked_atoms.push_back(cnt);
			cnt += topol.monomers[i].nAtoms;
		}
	}
	
	for (unsigned int i=0;i<keywords.size();i++){
		out_file << keywords[i];
	}
	out_file <<"\n" << endl;
	out_file << endl;
	out_file << std::fixed;
	out_file.precision(3);
	
	cnt = 0;	
	for(int i=0;i<mol.nAtoms;i++){
		std::string marking = "";
		if( i == marked_atoms[cnt] ) { 
			marking = mark;
			cnt++;
		}
		out_file << mol.atoms[i].element
				 << marking
				 << " "
				 << mol.atoms[i].xc
				 << " 1 "
				 << mol.atoms[i].yc
				 << " 1 "
				 << mol.atoms[i].zc
				 << "\n";
		marking = "";
	}
	
	out_file.close();
} 
/*************************************************************/
void mopac_input::write_file(molecule& mol, std::string out_name ){
	
	out_name +=".mop";
	
	out_file.open( out_name.c_str() );
	
	for (int i=0;i<keywords.size();i++){
		out_file << keywords[i];
	}
	out_file <<"\n" << endl;
	out_file << endl;
	out_file << std::fixed;
	out_file.precision(3);
	
	for(int i=0;i<mol.nAtoms;i++){
		out_file << mol.atoms[i].element 
				 << " "
				 << mol.atoms[i].xc
				 << " 1 "
				 << mol.atoms[i].yc
				 << " 1 "
				 << mol.atoms[i].zc
				 << "\n";
	}
	
	out_file.close();
	
	/*
	 * Fazer um scritp em python para rodar os cálculos 
	 * em paralelo e comprimir os resultados
	 */
}
/*************************************************************/
void mopac_input::read_from_input(	const char* inp_file,
									std::string out_name){
	
}
/*************************************************************/
std::ostream& operator<<(std::ostream& out, const mopac_input& obj){
	
}
/*************************************************************/
void mopac_input::print(){
	
}
/*************************************************************/
void UnitTest_mopac_input(){
	
}
////////////////////////////////////////////////////////////////