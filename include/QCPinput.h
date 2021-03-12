//QCPinput.h

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

#ifndef QCP_INPUT
#define QCP_INPUT

#include <string>

enum package{
	MOPAC,
	GAMESS,
	ORCA,
	
	NumOfPackages
}

/****************************************/

class QCPinput{
	public:
		package program;
		unsigned int base_multi;
		int base_charge;
		std::string QMmethod;
		std::string base_basis;
		std::string geo_ext;
		QCPinput();
		QCPinput(std::string _folder, std::string _geo_ext, std::string _basis, std::string _method)
		~QCPinput();
		QCPinput(const QCPinput& rhs) delete;
		QCPinput& operator=(const QCPinput& rhs) delete;
		void make_input_from_folder(int QCP,unsigned int bs_mlt,int bs_chg);
		// for finite differences calculations
		void make_input_from_folder_FD(int QCP,unsigned int bs_mlt,int bs_chg,int chg_diff); 		
};


#endif 