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
	keywords.emplace_back(" 1SCF ALLVECS VECTOR AUX LARGE");
	keywords.emplace_back(" charge=");
	keywords.emplace_back( std::to_string(charge) );
	
	if ( MOZYME ) keywords.emplace_back(" MOZYME");
	if ( COSMO  ) keywords.emplace_back(" eps=78.4");
	
	charge = chg;
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
	
	for (int i=0;i<keywords.size();i++){
		out_file << keywords[i];
	}
	out_file << endl;
	out_file << endl;
}
/*************************************************************/
void mopac_input::write_file(molecule& mol, std::string out_name ){
	
	out_name +=".mop";
	
	out_file.open( out_name.c_str() );
	
	for (int i=0;i<keywords.size();i++){
		out_file << keywords[i];
	}
	out_file << endl;
	out_file << endl;
	
	for(int i=0;i<mol.nAtoms;i++){
		out_file << mol.atoms[i].element 
				 << " "
				 << mol.atoms[i].xc
				 << " 1 "
				 << mol.atoms[i].zc
				 << " 1 "
				 << mol.atoms[i].zc
				 << "\n";
	}
	
	out_file.close();
}
/*************************************************************/
void mopac_input::read_from_input(	const char* inp_file,
									std::string out_name){
	
}
////////////////////////////////////////////////////////////////