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
#include "../include/atom.h"
#include "../include/molecule.h" 
#include "../include/geometry.h"
#include "../include/orcaInput.h"

using std::string;
using std::vector;

/***************************************/
orcaInput::orcaInput()	:
	multi(1)			,
	charge(0)			,
	nprocs(4)			,
	qm_method("HF")		,
	rtyp("SinglePoint")	,
	basis("3-21G")		{
}	
/*****************************************************/
orcaInput::~orcaInput(){}
/*****************************************************/
orcaInput::orcaInput(int chg			,
					int mlt				,
					int _nprocs			,
					std::string _rtyp)	:
					multi(mlt)			,
					charge(chg)			,
					nprocs(_nprocs)		,
					qm_method("HF")		,
					rtyp(_rtyp)			,
					basis("3-21G")		{
	
	
}
/*****************************************************/
void orcaInput::write_inp(	const molecule& mol		,
							std::string mth			,
							std::string _basis		,
							std::string out_file)	{
	
	qm_method	= mth;
	basis		= _basis;
	out_file	+= ".inp";
	string PAL	= "!PAL";
	PAL			+= std::to_string(nprocs);
	PAL			+= "\n";

	if (nprocs == 1) {
		PAL = "";
	}
	
	
	if ( mol.nElectrons % 2 != 0 && std::abs(charge) %2 == 0 ){
		multi = 2;		
	}else if ( mol.nElectrons % 2 != 0 && std::abs(charge) %2 != 0 ){
		multi = 1;
	}
	
	if ( rtyp == "SinglePoint" ){ rtyp = ""; }
	
	bool zora = false;
	for(int i=0;i<mol.nAtoms;i++){
		if ( mol.atoms[i].aNmb > 36 && !zora ){
			basis +=" ZORA"; 
			zora = true;
		}
	}
	
	
	out_fl.open( out_file.c_str() );
	out_fl	<< PAL
			<< "! " 
			<< qm_method
			<< " "
			<< rtyp
			<< " "
			<< basis
			<< "\n"
			<< "%output \n"
			<< "print[p_mos] 1\n"
			<< "print[p_overlap] 5\n"
			<< "end #output\n"
			<< "* xyz "
			<< charge 
			<< " "
			<< multi
			<< " \n";
			
	
	for(int i=0;i<mol.nAtoms;i++){
		out_fl	<< mol.atoms[i].element
				<< " "
				<< mol.atoms[i].xc
				<< " "
				<< mol.atoms[i].yc
				<< " "
				<< mol.atoms[i].zc
				<< std::endl;
	}
	out_fl << "*";
	out_fl.close();
}
///////////////////////////////////////////////////////