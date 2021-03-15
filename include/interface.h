//interface.h

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

#ifndef INTERFACE_H
#define INTERFACE_H

#include <string>
#include <vector>

class interface{
	public:
		int m_argc;
		std::vector<std::string> m_argv;
		interface();
		~interface();
		interface(int argc, char* argv[]);
		interface(const interface& rhs) = delete;
		interface& operator=(const interface& rhs) = delete;		
		void run();
		void print_options();
		void help();
};

#endif