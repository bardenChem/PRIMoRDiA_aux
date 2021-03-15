//orcaInput.cpp

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
#include <experimental/filesystem>

#include "../include/global.h"
#include "../include/molecule.h" 
#include "../include/geometry.h"
#include "../include/mopacInput.h"

using std::string;
using std::vector;

/***************************************/
orcaInput::orcaInput():
	multi(1),
	charge(0),
	nprocs(4),
	rtyp("SinglePoint"),
	basis("3-21G"){
	
}	
/*****************************************************/
orcaInput::~orcaInput(){
	
}
/*****************************************************/
void orcaInput::init(	const system& molecule	,
						int chg					, 
						int mlt					,	 
						std::string& mth		,
						std::string& _rtyp		,
						int _nprocs				,
						std::string _basis)		{
	
	multi	= mlt;
	charge	= chg;
	method	= mth;
	rtyp	= _rtyp;
	nprocs	= _nprocs;
	basis	= _basis;
						
	vector<string> keywords;
	keywords.emplace_back("!PAL");
	keywords.emplace_back( std::to_string(_nprocs) );
	keywords.emplace_back("\n");
	keywords.emplace_back("!"+method);
	keywords.emplace_back(" "+rtyp);
	keywords.emplace_back(" "+basis);
	keywords.emplace_back("\n");	
	
}
/******************************************************/
void orcaInput::write_input_file(std::string out_file){
	
}
///////////////////////////////////////////////////////