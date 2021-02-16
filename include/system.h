//system.h

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

#ifndef SYSTEM
#define SYSTEM

#include <vector>
#include <string>

/*********************************************************************/
class atom{
	public:
		double xc;
		double yc;
		double zc;
		std::string element;
		unsigned int nmb;
		atom();
		~atom();
		atom(double xc, double yc, double zc, std::string type);
		atom(const atom& rhs);
		atom& operator=(const atom& rhs);
		atom(atom&& rhs) noexcept;
		atom&& operator=(atom&& rhs) noexcept;
};
/*********************************************************************/
class pdbAtom{
	public:
		std::string atom_name;
		unsigned int indx;
		std::string res_name;
		std::string chain_name;
		bool sideC;
		atom atom_info;
		pdbAtom();
		~pdbAtom();
		pdbAtom(const pdbAtom& rhs);
		pdbAtom& operator=(const pdbAtom& rhs);
		pdbAtom(pdbAtom&& rhs) noexcept;
		pdbAtom&& operator=(pdbAtom&& rhs) noexcept;
};
/*********************************************************************/
class residue{
	public:
		std::string type;
		residue();
		~residue();
};
/*********************************************************************/
class chain{
	public:
		std::string index;
		std::vector<residue> monomers;
		chain();
		~chain();
};
/*********************************************************************/
class system{
	public:
		std::string name;
		std::string type;
		
};


#endif