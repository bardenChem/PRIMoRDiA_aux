//pdbAtom.h

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

#ifndef PDB_ATOM_H
#define PDB_ATOM_H

#include <string>

class pdbAtom{
	public:
		std::string atom_name;
		unsigned int indx;
		std::string res_name;
		int res_indx;
		std::string chain_name;
		float occupancy;
		double b_factor;
		bool sideC;
		double xc;
		double yc;
		double zc;
		pdbAtom();
		~pdbAtom();
		pdbAtom(const pdbAtom& rhs);
		pdbAtom& operator=(const pdbAtom& rhs);
		pdbAtom(pdbAtom&& rhs) noexcept;
		pdbAtom& operator=(pdbAtom&& rhs) noexcept;
		pdbAtom(std::string& pdb_line);
		bool operator==(const pdbAtom& rhs);
		bool is_hydrogen();
		
};

#endif 