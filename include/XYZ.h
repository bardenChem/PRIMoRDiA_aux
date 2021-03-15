//XYZ.h

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

#ifndef XYZ_H
#define XYZ_H

#include <vector>
#include <string>

class atom;
class molecule;

class XYZ{
	public:
		unsigned int nAtoms;
		std::string name;
		std::vector<atom> atoms;
		std::string commentary;
		XYZ();
		~XYZ();
		XYZ(const char* xyz_file);
		XYZ(const molecule& mol );
		XYZ(const XYZ& rhs);
		XYZ& operator=(const XYZ& rhs);
		XYZ(XYZ&& rhs) noexcept;
		XYZ& operator=(XYZ&& rhs) noexcept;
		void write_xyz(std::string out_name);
		molecule get_molecule();
};

#endif 