//orcaInput.h
//Declarations for functions and source file with constant values. 

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

#ifndef ORCA_INP_H
#define ORCA_INP_H

#include <string>
#include <fstream>

class molecule;

//=======================================
enum OrcaRunType{
	SP			,
	oOPT		,	
	NUmOfopts
};
/********************************************************************/
class orcaInput{
	public:
		unsigned int multi;
		int charge;
		unsigned int nprocs;
		std::string qm_method;
		std::string rtyp;
		std::string basis;
		std::ofstream out_fl;
		orcaInput();	
		~orcaInput();
		orcaInput(int chg, int mlt,int _nprocs,std::string _rtyp);
		orcaInput(const orcaInput& rhs) = delete;
		orcaInput& operator=(const orcaInput& rhs) = delete;
		void write_inp(const molecule& mol, std::string mth,std::string _basis,std::string out_file);
};

///////////////////////////////////////////////////////////////////////
#endif