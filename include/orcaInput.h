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

#ifndef ORCA_INP
#define ORCA_INP

#include <string>
#include "../include/global.h"

class system;

//=======================================
enum OrcaRunType{
	SP,
	OPT,	
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
		std::ostream out_fl;
		std::string basis;
		orcaInput();	
		~orcaInput();
		orcaInput(const orcaInput& rhs) = delete;
		orcaInput& operator=(const orcaInput& rhs) = delete;
		void init(const system& molecule,int chg, int mlt, std::string& mth, std::string& _rtyp, int _nprocs,std::string _basis);
		void write_input_file(std::string out_file);
};

///////////////////////////////////////////////////////////////////////
#endif